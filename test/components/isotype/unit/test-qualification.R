# Check C-call qualification boundaries and preserved common-QC observations.
# Run: python test/run_tests.py --component isotype --level unit

load_test_modules(c("table_io", "gene_calls", "isotype_calls", "isotype_resolution"))

# ------------------------------------------------------------ C callability ---

# The configured boundary masks a field, not a receptor: a C call without enough
# aligned evidence is unusable as an isotype, while the record itself remains a
# perfectly good receptor.
describe("c_call_qc_flags at a configured floor of 40 informative positions", {
  flags <- c_call_qc_flags(
    c("IGHG1*01", "IGHG1*01", "IGHG1*01", ""),
    c("IGH", "IGH", "IGK", "IGH"),
    c(39, 40, 80, NA),
    40
  )

  it("rejects a call one position below the floor",
     expect_identical(flags$callable[[1]], FALSE))
  it("accepts a call sitting exactly on the floor",
     expect_identical(flags$callable[[2]], TRUE))
  it("rejects a heavy-chain call on a kappa record however deep the evidence",
     expect_identical(flags$callable[[3]], FALSE))
  it("rejects an absent call",
     expect_identical(flags$callable[[4]], FALSE))
  # Absent and failed are different readings, so they are reported separately.
  it("counts every non-empty call as present",
     expect_identical(flags$present[1:3], rep(TRUE, 3)))
  it("counts an empty call as absent",
     expect_identical(flags$present[[4]], FALSE))
})


test_that("zero C coverage removes only the length floor", {
  flags <- c_call_qc_flags(
    c("IGHM", "IGHM", "", "TRBC1"), rep("IGH", 4), c(0, NA, 0, 10), 0)
  expect_identical(flags$callable, c(TRUE, FALSE, FALSE, FALSE))
})

test_that("the shared call table retains failed IGH records and original calls", {
  records <- data.frame(
    library_id = "0460", sequence_id = c("s3", "s1", "s2", "s4", "s5", "t1"),
    locus = c(rep("IGH", 5), "TRB"),
    c_call = c("IGHG1*01", "IGHM*01", "IGHA1,IGHA2", "", "TRBC1", "TRBC1"),
    c_informative_positions = c(39, 40, 80, NA, 80, 80),
    consensus_count = c(3, 5, 7, 11, 13, 17),
    stringsAsFactors = FALSE
  )
  calls <- build_qualified_isotype_calls(records, "0460", 40)
  expect_identical(names(calls), ISOTYPE_CALL_COLUMNS)
  expect_identical(calls$library_id, rep("0460", 5))
  expect_identical(calls$sequence_id, paste0("s", 1:5))
  expect_identical(calls$c_call,
                   c("IGHM*01", "IGHA1,IGHA2", "IGHG1*01", "", "TRBC1"))
  expect_identical(calls$qualified_c_call,
                   c("IGHM*01", "IGHA1,IGHA2", rep(NA_character_, 3)))
  expect_identical(calls$c_call_present, c(TRUE, TRUE, TRUE, FALSE, TRUE))
  expect_identical(calls$c_call_enough_evidence, c(TRUE, TRUE, FALSE, FALSE, TRUE))
  expect_identical(calls$c_call_locus_match, c(TRUE, TRUE, TRUE, TRUE, FALSE))
  expect_identical(calls$c_call_callable, c(TRUE, TRUE, FALSE, FALSE, FALSE))
  igh <- records[records$locus == "IGH", ]
  joined <- join_qualified_isotype_calls(igh, calls)
  expect_identical(joined$sequence_id, igh$sequence_id)
  expect_identical(joined$consensus_count, igh$consensus_count)
  expect_identical(joined$qualified_c_call, c(NA_character_, "IGHM*01", "IGHA1,IGHA2",
                                   NA_character_, NA_character_))
})

test_that("call qualification round trips with empty calls and text identifiers", {
  records <- data.frame(
    library_id = "0460", sequence_id = c("0001", "0002"), locus = "IGH",
    c_call = c("IGHM*01", NA_character_), c_informative_positions = c(40, NA)
  )
  calls <- build_qualified_isotype_calls(records, "0460", 40)
  path <- tempfile(fileext = ".tsv")
  on.exit(unlink(path), add = TRUE)
  write_pipeline_table(calls, path)
  restored <- read_pipeline_table(path)
  expect_identical(restored$library_id, rep("0460", 2))
  expect_identical(restored$sequence_id, c("0001", "0002"))
  expect_identical(restored$qualified_c_call, c("IGHM*01", ""))
  expect_identical(restored$c_call_callable, c("TRUE", "FALSE"))
  joined <- join_qualified_isotype_calls(records, restored)
  expect_identical(resolve_isotype(joined$qualified_c_call), c("IgM", NA_character_))
})

test_that("an empty IGH population still has the complete sidecar header", {
  records <- data.frame(library_id = "A", sequence_id = "t1", locus = "TRB",
                         c_call = "TRBC1", c_informative_positions = 80)
  empty <- build_qualified_isotype_calls(records, "A", 40)
  expect_identical(names(empty), ISOTYPE_CALL_COLUMNS)
  expect_equal(nrow(empty), 0L)
  expect_identical(empty$qualified_c_call, character())
  expect_identical(empty$c_call_callable, logical())
  expect_identical(build_qualified_isotype_calls(records[FALSE, ], "A", 40), empty)
})

test_that("malformed library keys cannot produce a qualification sidecar", {
  records <- data.frame(library_id = "A", sequence_id = c("s1", "s2"), locus = "IGH",
                         c_call = "IGHM*01", c_informative_positions = 40)
  expect_error(build_qualified_isotype_calls(records, "B", 40), "library_id other than")
  records$sequence_id <- c("s1", "s1")
  expect_error(build_qualified_isotype_calls(records, "A", 40), "unique")
  records$sequence_id <- c("s1", "")
  expect_error(build_qualified_isotype_calls(records, "A", 40), "non-missing")
})

test_that("joining qualification preserves the original annotation", {
  records <- data.frame(library_id = "NA", sequence_id = c("a", "b"), locus = "IGH",
    c_call = c("IGHM", "IGHA1"), c_informative_positions = c(40, 1))
  calls <- build_qualified_isotype_calls(records, "NA", 40)
  joined <- join_qualified_isotype_calls(records, calls)
  expect_identical(joined$c_call, records$c_call)
  expect_identical(joined$qualified_c_call, c("IGHM", NA_character_))
})

igh_records <- function(sequence_id = c("s1", "s2", "s3")) {
  data.frame(library_id = "WTS_1_S1", sequence_id = sequence_id,
             consensus_count = c(3, 5, 2)[seq_along(sequence_id)],
             stringsAsFactors = FALSE)
}
qualified_igh_calls <- function(sequence_id = c("s1", "s2", "s3"),
                         qualified_c_call = c("IGHM*01", "IGHG1,IGHG3", "IGHA1")) {
  data.frame(library_id = "WTS_1_S1", sequence_id = sequence_id,
             qualified_c_call = qualified_c_call[seq_along(sequence_id)], stringsAsFactors = FALSE)
}

describe("join_qualified_isotype_calls", {
  it("attaches each record's constant-region call", {
    joined <- join_qualified_isotype_calls(igh_records(), qualified_igh_calls())
    expect_identical(joined$qualified_c_call, c("IGHM*01", "IGHG1,IGHG3", "IGHA1"))
  })
  it("matches on the key, not on row order", {
    shuffled <- qualified_igh_calls(c("s3", "s1", "s2"),
                             c("IGHA1", "IGHM*01", "IGHG1,IGHG3"))
    expect_identical(
      join_qualified_isotype_calls(igh_records(), shuffled)$qualified_c_call,
      c("IGHM*01", "IGHG1,IGHG3", "IGHA1")
    )
  })
  it("leaves every other column alone", {
    joined <- join_qualified_isotype_calls(igh_records(), qualified_igh_calls())
    expect_identical(names(joined),
                     c("library_id", "sequence_id", "consensus_count", "qualified_c_call"))
  })
  it("gives a library with no IGH records an empty call column", {
    empty <- igh_records()[0, , drop = FALSE]
    joined <- join_qualified_isotype_calls(empty, qualified_igh_calls()[0, , drop = FALSE])
    expect_identical(joined$qualified_c_call, character(0))
  })
  it("refuses a duplicate key on the AIRR side", {
    expect_error(
      join_qualified_isotype_calls(igh_records(c("s1", "s1", "s3")),
                                        qualified_igh_calls()),
      "duplicate sequence_id"
    )
  })
  it("refuses a duplicate key on the sidecar side", {
    expect_error(
      join_qualified_isotype_calls(
        igh_records(), qualified_igh_calls(c("s1", "s1", "s3"))),
      "duplicate sequence_id"
    )
  })
  # A record the sidecar does not describe would silently get no isotype, and a
  # sidecar row with no record would mean the sidecar describes something else.
  it("refuses a record the sidecar does not describe", {
    expect_error(
      join_qualified_isotype_calls(igh_records(),
                                        qualified_igh_calls(c("s1", "s2"))),
      "missing=1, extra=0"
    )
  })
  it("refuses a sidecar row with no matching record", {
    expect_error(
      join_qualified_isotype_calls(igh_records(c("s1", "s2")),
                                        qualified_igh_calls()),
      "missing=0, extra=1"
    )
  })
  # Both tables describe one library, so a sequence_id is the whole key. The stage
  # is what establishes that, by refusing any row from another library.
  it("keys on the sequence_id alone", {
    records <- data.frame(library_id = "A", sequence_id = c("x|y", "y"),
                          consensus_count = c(1, 2), stringsAsFactors = FALSE)
    features <- data.frame(library_id = "A", sequence_id = c("y", "x|y"),
                           qualified_c_call = c("IGHA1", "IGHM"), stringsAsFactors = FALSE)
    expect_identical(join_qualified_isotype_calls(records, features)$qualified_c_call,
                     c("IGHM", "IGHA1"))
  })
})

test_that("C-call completeness follows scalar-or-record-length rules", {
  flags <- c_call_qc_flags(rep("IGHM", 3), "IGH", 40, 40, TRUE)
  expect_identical(flags$callable, rep(TRUE, 3))
  flags <- c_call_qc_flags(rep("IGHM", 3), "IGH", 40, 40, c(TRUE, FALSE, TRUE))
  expect_identical(flags$callable, c(TRUE, FALSE, TRUE))
  for (complete in list(logical(), c(TRUE, FALSE), rep(TRUE, 4))) {
    expect_error(c_call_qc_flags(rep("IGHM", 3), "IGH", 40, 40, complete),
                 "equal lengths or length one")
  }
  for (complete in list(TRUE, logical())) {
    empty <- c_call_qc_flags(character(), character(), numeric(), 40, complete)
    expect_identical(nrow(empty), 0L)
    expect_identical(empty$callable, logical())
  }
  expect_error(c_call_qc_flags(character(), character(), numeric(), 40,
                               c(TRUE, FALSE)), "equal lengths or length one")
})
