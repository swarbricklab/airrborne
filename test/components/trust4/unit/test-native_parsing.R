# Check native TRUST4 imputation and original record support.
# Run: python test/run_tests.py --component trust4 --level unit

load_test_modules(c("trust4_parsing"))

WORK_DIR <- tempfile("trust4-parsing-unit-")
dir.create(WORK_DIR)

# A headerless 13-column cdr3.out. Only columns 1-2 (consensus_id, index) and 10
# (the CDR3 score) are read; the rest are filled so the width is realistic.
write_cdr3_out <- function(path, consensus_id, index, score, n_columns = 13L) {
  df <- data.frame(
    cid = consensus_id, idx = index, v = "IGHV1-8*01", d = ".", j = "IGHJ5*02",
    c = "IGHM", cdr1 = "AAA", cdr2 = "CCC", cdr3 = "TGTGCC", score = score,
    reads = "1", germline_similarity = "95.0", complete_vdj = "1",
    stringsAsFactors = FALSE
  )
  write.table(df[, seq_len(n_columns), drop = FALSE], path, sep = "\t",
              row.names = FALSE, col.names = FALSE, quote = FALSE)
}

# --- The imputed-CDR3 sentinel ------------------------------------------------

test_that("the imputed sentinel and its tolerance are TRUST4's own values", {
  expect_identical(TRUST4_IMPUTED_CDR3_SCORE, 0.01)
  expect_identical(TRUST4_IMPUTED_CDR3_TOLERANCE, 1e-6)
})

# --- read_trust4_cdr3_scores() ------------------------------------------------

test_that("an empty CDR3 file is valid only for an empty AIRR population", {
  path <- file.path(WORK_DIR, "empty_cdr3.out")
  file.create(path)
  scores <- read_trust4_cdr3_scores(path)
  expect_identical(scores, data.frame(consensus_id = character(),
    cdr3_key = character(), cdr3_score = numeric()))
  expect_identical(flag_imputed_cdr3_records(character(), scores$cdr3_key,
    scores$cdr3_score), logical())
  expect_error(flag_imputed_cdr3_records("assemble0_0", scores$cdr3_key,
    scores$cdr3_score), "missing a CDR3 score")
  expect_error(suppressWarnings(read_trust4_cdr3_scores(paste0(path, ".missing"))),
    "cannot open")
})

describe("a three-call cdr3.out", {
  path <- file.path(WORK_DIR, "cdr3.out")
  write_cdr3_out(path, c("assemble0", "assemble0", "assemble9"),
                 c("0", "1", "3"), c("1.00", "1.00", "0.83"))
  scores <- read_trust4_cdr3_scores(path)

  it("joins columns 1 and 2 into a score lookup key", {
    expect_identical(scores$cdr3_key,
                     c("assemble0_0", "assemble0_1", "assemble9_3"))
  })
  it("reads column 10 as the numeric score", {
    expect_identical(scores$cdr3_score, c(1.00, 1.00, 0.83))
  })
  it("keeps the assembly each call came from", {
    expect_identical(scores$consensus_id,
                     c("assemble0", "assemble0", "assemble9"))
    # Two contigs, three calls: this is the QC audit's assembly-level denominator.
    expect_identical(length(unique(scores$consensus_id)), 2L)
    expect_identical(nrow(scores), 3L)
  })
})

test_that("a single-call cdr3.out reads as one row", {
  path <- file.path(WORK_DIR, "cdr3_single.out")
  write_cdr3_out(path, "assemble0", "0", "1.00")
  scores <- read_trust4_cdr3_scores(path)
  expect_identical(scores$cdr3_key, "assemble0_0")
  expect_identical(scores$cdr3_score, 1.00)
})

test_that("an unparseable score comes back NA rather than raising", {
  path <- file.path(WORK_DIR, "cdr3_bogus.out")
  write_cdr3_out(path, "assemble0", "0", "bogus")
  expect_identical(read_trust4_cdr3_scores(path)$cdr3_score, NA_real_)
})

test_that("a narrower cdr3.out stops the stage", {
  path <- file.path(WORK_DIR, "cdr3_short.out")
  write_cdr3_out(path, "assemble0", "0", "1.00", n_columns = 12L)
  expect_error(read_trust4_cdr3_scores(path),
               "TRUST4 cdr3.out has fewer than 13 columns", fixed = TRUE)
})

test_that("two calls sharing one key stop the stage", {
  path <- file.path(WORK_DIR, "cdr3_dup.out")
  write_cdr3_out(path, c("assemble0", "assemble0"), c("0", "0"),
                 c("1.00", "1.00"))
  expect_error(read_trust4_cdr3_scores(path),
               "duplicate <consensus_id>_<index> keys", fixed = TRUE)
})

# --- flag_imputed_cdr3_records() ----------------------------------------------

describe("classifying records against their CDR3 scores", {
  keys <- c("assemble0_0", "assemble5_0", "assemble9_3")
  values <- c(1.00, 0.01, 0.83)

  it("flags only the sentinel score", {
    expect_identical(flag_imputed_cdr3_records(keys, keys, values),
                     c(FALSE, TRUE, FALSE))
  })
  it("returns nothing for no records", {
    expect_identical(flag_imputed_cdr3_records(character(0), keys, values),
                     logical(0))
  })
  it("classifies one record on its own", {
    expect_true(flag_imputed_cdr3_records("assemble5_0", keys, values))
  })
  it("follows the score, not the row order", {
    expect_identical(
      flag_imputed_cdr3_records(keys[c(2, 3, 1)], keys, values),
      c(TRUE, FALSE, FALSE)
    )
  })
})

test_that("variant rank changes preserve the repeated contig imputation status", {
  keys <- c("assemble0_0", "assemble0_1", "assemble5_0", "assemble5_1")
  expect_identical(flag_imputed_cdr3_records(keys[c(2, 1, 4, 3)], keys,
                                            c(1, 1, 0.01, 0.01)),
                   c(FALSE, FALSE, TRUE, TRUE))
})

# The tolerance exists to absorb the decimal text round-trip, not to widen the
# rule: 0.00 and the smallest nonzero motif score 1/6 stay outside it.
test_that("the tolerance admits a rounding wobble and nothing else", {
  ids <- c("a_1", "a_2", "a_3", "a_4", "a_5")
  expect_identical(
    flag_imputed_cdr3_records(
      ids, ids, c(0.01, 0.01 + 5e-7, 0.01 + 2e-6, 0.00, 1 / 6)),
    c(TRUE, TRUE, FALSE, FALSE, FALSE)
  )
})

test_that("a record with no score in cdr3.out stops the stage", {
  expect_error(
    flag_imputed_cdr3_records(c("nope_1", "nope_2", "nope_3", "nope_4"),
                              "assemble0_0", 1.00),
    "missing a CDR3 score for 4 sequence_id(s), e.g. nope_1, nope_2, nope_3",
    fixed = TRUE
  )
})

test_that("a matched but non-numeric score stops the stage", {
  expect_error(
    flag_imputed_cdr3_records("assemble0_0", "assemble0_0", NA_real_),
    "non-numeric CDR3 score for a matched sequence_id", fixed = TRUE
  )
})
