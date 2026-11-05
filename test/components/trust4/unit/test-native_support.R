# Check TRUST4-reported read support and yield using input read-pair denominators.
# Run: python test/run_tests.py --component trust4 --level unit

load_test_modules(c("table_io", "repertoire_support"))

SUPPORT_COLUMNS <- c(
  "library_id", "trust4_consensus_count", "bcr_trust4_consensus_count",
  "tcr_trust4_consensus_count",
  "igh_trust4_consensus_count", "igk_trust4_consensus_count",
  "igl_trust4_consensus_count", "tra_trust4_consensus_count",
  "trb_trust4_consensus_count", "trd_trust4_consensus_count",
  "trg_trust4_consensus_count"
)

test_that("the loci are the seven canonical ones, B-cell first", {
  expect_identical(SUPPORT_BCR_LOCI, c("IGH", "IGK", "IGL"))
  expect_identical(SUPPORT_TCR_LOCI, c("TRA", "TRB", "TRD", "TRG"))
  expect_identical(SUPPORT_LOCI,
                   c("IGH", "IGK", "IGL", "TRA", "TRB", "TRD", "TRG"))
})

# Distinct per-locus values so a mis-ordered column would not go unnoticed.
BY_LOCUS <- stats::setNames(c(100, 20, 3, 7000, 800, 90, 5), SUPPORT_LOCI)

describe("a library with support on every locus", {
  row <- build_repertoire_support_row("WTS_1_S1", BY_LOCUS)

  it("is one row",
     expect_identical(nrow(row), 1L))
  it("carries the released columns in the released order",
     expect_identical(names(row), SUPPORT_COLUMNS))
  it("names the library",
     expect_identical(row$library_id, "WTS_1_S1"))
  it("sums the four T-cell loci",
     expect_equal(row$tcr_trust4_consensus_count, 7895))
  it("sums the three B-cell loci",
     expect_equal(row$bcr_trust4_consensus_count, 123))
  it("totals the two classes",
     expect_equal(row$trust4_consensus_count, 8018))
  it("reports each locus under its own column",
     expect_equal(unname(unlist(row[, paste0(tolower(SUPPORT_LOCI), "_trust4_consensus_count")])),
                  unname(BY_LOCUS)))
})

# calculate_isotype.R and this stage both promise the shared aggregate exactly
# one row per library, so an empty assembly is a zero row rather than no row.
describe("a library that assembled nothing", {
  row <- build_repertoire_support_row(
    "WTS_2_S1", stats::setNames(rep(0, length(SUPPORT_LOCI)), SUPPORT_LOCI)
  )

  it("still reports one row",
     expect_identical(nrow(row), 1L))
  it("still carries the released columns",
     expect_identical(names(row), SUPPORT_COLUMNS))
  it("reports zero rather than a missing total",
     expect_equal(row$trust4_consensus_count, 0))
})

test_that("support on a single locus reaches only that class", {
  by_locus <- stats::setNames(rep(0, length(SUPPORT_LOCI)), SUPPORT_LOCI)
  by_locus[["TRB"]] <- 42
  row <- build_repertoire_support_row("WTS_3_S1", by_locus)
  expect_equal(row$tcr_trust4_consensus_count, 42)
  expect_equal(row$bcr_trust4_consensus_count, 0)
  expect_equal(row$trb_trust4_consensus_count, 42)
  expect_equal(row$tra_trust4_consensus_count, 0)
})

# The vector is read by name, not by position, because the caller builds it by
# scattering a group-by result into a zero-filled vector.
test_that("the per-locus vector is read by name rather than by position", {
  expect_identical(
    build_repertoire_support_row("WTS_4_S1", BY_LOCUS),
    build_repertoire_support_row("WTS_4_S1", rev(BY_LOCUS))
  )
})

test_that("a locus missing from the vector is an error rather than a silent zero", {
  expect_error(
    build_repertoire_support_row("WTS_5_S1", BY_LOCUS[SUPPORT_LOCI != "TRD"])
  )
})

test_that("yield uses TRUST4-reported support and exact input-pair depth at all seven loci", {
  support <- build_repertoire_support_row("NA", setNames(0:6, SUPPORT_LOCI))
  depth <- function(pairs) data.frame(library_id = "NA", n_fastq_file_pairs = "2",
                                      input_read_pairs = pairs, total_bases = "600000000",
                                      stringsAsFactors = FALSE)
  yield <- add_repertoire_yield(support, depth("2000000"))
  columns <- paste0(tolower(SUPPORT_LOCI), "_yield_per_million_read_pairs")
  expect_identical(names(yield), c(
    "library_id", "n_fastq_file_pairs", "input_read_pairs", "input_bases", "yield_status",
    "trust4_consensus_count", "bcr_trust4_consensus_count", "tcr_trust4_consensus_count",
    as.vector(rbind(paste0(tolower(SUPPORT_LOCI), "_trust4_consensus_count"), columns))))
  expect_identical(c(yield$library_id, yield$n_fastq_file_pairs, yield$input_bases),
                   c("NA", "2", "600000000"))
  expect_equal(unname(unlist(yield[columns])), (0:6) / 2)
  expect_identical(yield$yield_status, "available")
  expect_true(all(is.na(add_repertoire_yield(support, depth("0"))[columns])))
  expect_error(add_repertoire_yield(support, depth("1.5")), "integer")
})
