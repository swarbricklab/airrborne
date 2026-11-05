# Check common record QC, including exact boundaries and overlapping failures.
# Run: python test/run_tests.py --component record_qc --level unit

load_test_modules(c("table_io", "gene_calls", "record_qc"))

# ------------------------------------------- one common post-IgBLAST QC rule ---

# Nine records against the shared gate: three are retained and each of the other
# six is engineered to fail exactly one criterion. Record 4 has a valid short
# alignment, and record 9 carries the V gene TRA and TRD share; both are kept. d_call is absent
# throughout, so these records exercise only the V/J, N-fraction, productivity,
# junction-length and support criteria; the D criterion has its own frame below.
qc_criterion_records <- function() {
  data.frame(
    locus = c(rep("IGH", 8), "TRD"),
    v_call = c("IGHV1-2*01", "IGKV1-5*01", rep("IGHV1-2*01", 6),
               "TRAV38-2/DV8*01"),
    j_call = c("IGHJ4*02", "IGHJ4*02", "IGKJ1*01", rep("IGHJ4*02", 5),
               "TRDJ1*01"),
    d_call = rep("", 9),
    sequence_alignment = c(
      paste0(rep("A", 200), collapse = ""),
      paste0(rep("A", 200), collapse = ""),
      paste0(rep("A", 200), collapse = ""),
      paste0(rep("A", 199), collapse = ""),
      paste0(paste0(rep("A", 200), collapse = ""),
             paste0(c(rep("N", 12), rep("n", 11)), collapse = "")),
      paste0(rep("A", 200), collapse = ""),
      paste0(rep("A", 200), collapse = ""),
      paste0(rep("A", 200), collapse = ""),
      paste0(rep("A", 200), collapse = "")
    ),
    productive = c(rep("T", 5), "F", "T", "T", "T"),
    junction_length = c(rep(48, 6), 47, 48, 48),
    consensus_count = c(rep(1, 7), 0, 1),
    stringsAsFactors = FALSE
  )
}

describe("airr_qc_flags on nine records with individual criterion failures", {
  flags <- airr_qc_flags(qc_criterion_records(), 0.1)
  passes <- airr_qc_pass(flags)

  it("keeps a record satisfying every criterion",
     expect_identical(passes[[1]], TRUE))
  it("drops a record whose V call belongs to another locus",
     expect_identical(passes[[2]], FALSE))
  it("drops a record whose J call belongs to another locus",
     expect_identical(passes[[3]], FALSE))
  # Common QC has no alignment-length floor. Regional SHM eligibility owns
  # coverage requirements, so this valid 199-position record is retained.
  it("keeps a short alignment that satisfies common record QC",
     expect_identical(passes[[4]], TRUE))
  it("drops a record whose alignment is more than a tenth ambiguous",
     expect_identical(passes[[5]], FALSE))
  it("drops a non-productive record",
     expect_identical(passes[[6]], FALSE))
  it("drops a record whose junction length is not a multiple of three",
     expect_identical(passes[[7]], FALSE))
  it("drops a record with no read support",
     expect_identical(passes[[8]], FALSE))
  it("keeps a TRD record called with the V gene TRA and TRD share",
     expect_identical(passes[[9]], TRUE))

  # Record 5 carries 12 uppercase and 11 lowercase N over 223 aligned bases,
  # which clears the 0.1 ceiling only if both cases are counted.
  it("counts uppercase and lowercase N alike",
     expect_identical(flags$low_n[[5]], FALSE))
  it("finds every other record's ambiguity acceptable",
     expect_identical(flags$low_n[-5], rep(TRUE, 8)))

  it("has no whole-alignment length criterion",
     expect_false("enough_alignment" %in% names(flags)))
})

test_that("the aligned N ceiling includes its boundary and ignores alignment gaps", {
  records <- qc_criterion_records()[rep(1, 5), , drop = FALSE]
  records$sequence_alignment <- c(
    paste0(strrep("A", 91), strrep("N", 9), ".-"),
    paste0(strrep("A", 90), "NNNNNnnnnn", ".-"),
    paste0(strrep("A", 89), strrep("N", 11), ".-"), "", ".--."
  )
  expect_identical(airr_qc_flags(records, 0.1)$low_n,
                   c(TRUE, TRUE, FALSE, FALSE, FALSE))
  expect_identical(airr_qc_flags(records, 1)$low_n,
                   c(TRUE, TRUE, TRUE, FALSE, FALSE))
})

test_that("common QC measures junction length separately from grouping eligibility", {
  records <- qc_criterion_records()[rep(1, 7), , drop = FALSE]
  records$junction_length <- c(0, 3, 4, -3, Inf, NA, 48)
  expect_identical(airr_qc_flags(records)$junction_mod3,
                   c(TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, TRUE))
})

# Two callers need the conjunction of the criteria: the gate that keeps records,
# and the population audit that reports what the gate cost. Both read it from
# airr_qc_fail_counts(), so there is one evaluation rather than two that could
# drift. A fail count is only meaningful while every flag is logical and non-NA
# -- an NA would otherwise make the count NA, the pass mask NA, and the filtered
# table silently gain a row of NAs -- so that precondition is asserted, and the
# last two cases below are what hold it.

describe("the one evaluation of the common QC rule", {
  flags <- airr_qc_flags(qc_criterion_records(), 0.1)

  it("counts the criteria each record fails",
     expect_identical(airr_qc_fail_counts(flags)[1:3], c(0L, 1L, 1L)))
  it("passes exactly the records that fail nothing",
     expect_identical(airr_qc_pass(flags), airr_qc_fail_counts(flags) == 0L))
  it("agrees with the conjunction read criterion by criterion",
     expect_identical(airr_qc_pass(flags),
                      Reduce(`&`, as.list(flags))))

  it("returns nothing for a table with no records", {
    empty <- flags[0L, , drop = FALSE]
    expect_identical(airr_qc_fail_counts(empty), integer(0))
    expect_identical(airr_qc_pass(empty), logical(0))
  })

  it("refuses a flag column that is not logical", {
    text_flags <- flags
    text_flags$productive <- as.character(text_flags$productive)
    expect_error(airr_qc_pass(text_flags), regexp = "productive")
  })

  it("refuses an NA flag, which would otherwise pass NA rows into the output", {
    na_flags <- flags
    na_flags$low_n[2] <- NA
    expect_error(airr_qc_pass(na_flags), regexp = "low_n")
  })
})

# A D call is optional; an incoherent one is not. The TRA row carrying a TRDD
# call is the case this criterion exists for: IgBLAST does report it, and it is
# the classic TRA/TRD ambiguity rather than a clean TRA rearrangement.
d_call_records <- function() {
  data.frame(
    locus = c("IGH", "IGH", "TRB", "TRB", "TRA", "TRA", "IGK", "TRD"),
    v_call = c(rep("IGHV1-2*01", 2), rep("TRBV5-1*01", 2),
               rep("TRAV38-2/DV8*01", 2), "IGKV1-5*01", "TRDV1*01"),
    j_call = c(rep("IGHJ4*02", 2), rep("TRBJ2-7*01", 2),
               rep("TRAJ1*01", 2), "IGKJ1*01", "TRDJ1*01"),
    d_call = c("IGHD3-10*01", "", "TRBD1*01", "TRDD2*01",
               "", "TRDD2*01", "IGHD3-10*01", "TRDD2*01"),
    sequence_alignment = rep(paste0(rep("A", 200), collapse = ""), 8),
    productive = rep("T", 8),
    junction_length = rep(48, 8),
    consensus_count = rep(1, 8),
    stringsAsFactors = FALSE
  )
}

describe("the D criterion of airr_qc_flags", {
  d_locus_match <- airr_qc_flags(d_call_records(), 0.1)$d_locus_match

  it("accepts a heavy-chain D segment on an IGH record",
     expect_identical(d_locus_match[[1]], TRUE))
  it("accepts an absent D call, which IgBLAST often cannot place",
     expect_identical(d_locus_match[[2]], TRUE))
  it("accepts a TRB D segment on a TRB record",
     expect_identical(d_locus_match[[3]], TRUE))
  it("rejects a TRD D segment on a TRB record",
     expect_identical(d_locus_match[[4]], FALSE))
  it("accepts a TRA record with no D call, since TRA joins V straight to J",
     expect_identical(d_locus_match[[5]], TRUE))
  it("rejects the TRA record carrying a TRDD call, the classic TRA/TRD ambiguity",
     expect_identical(d_locus_match[[6]], FALSE))
  it("rejects a D call on IGK, which rearranges no D segment",
     expect_identical(d_locus_match[[7]], FALSE))
  it("accepts a TRD D segment on a TRD record",
     expect_identical(d_locus_match[[8]], TRUE))

  it("leaves the V/J criterion unmoved by an incoherent D call",
     expect_identical(airr_qc_flags(d_call_records(), 0.1)$locus_match,
                      rep(TRUE, 8)))
})

# A short alignment is a normal outcome of 5' V truncation in bulk RNA-seq, so a
# short but otherwise valid record is retained at every locus.
short_alignment_records <- function() {
  data.frame(
    locus = c("TRB", "IGH"),
    v_call = c("TRBV5-1*01", "IGHV1-2*01"),
    j_call = c("TRBJ2-7*01", "IGHJ4*02"),
    d_call = c("TRBD1*01", ""),
    sequence_alignment = rep(paste0(rep("A", 150), collapse = ""), 2),
    productive = c("T", "T"),
    junction_length = c(48, 48),
    consensus_count = c(1, 1),
    stringsAsFactors = FALSE
  )
}

describe("a 150-position alignment against the common QC gate", {
  records <- short_alignment_records()
  passes <- airr_qc_pass(airr_qc_flags(records, 0.1))

  it("is retained at TRB", expect_identical(passes[[1]], TRUE))
  it("is retained at IGH", expect_identical(passes[[2]], TRUE))
})
