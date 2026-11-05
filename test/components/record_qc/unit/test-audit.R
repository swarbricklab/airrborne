# Check observed QC populations, explicit units and conserved read support.
# Run: python test/run_tests.py --component record_qc --level unit

load_test_modules(c(
  "table_io", "gene_calls", "assembly_backbone", "rearrangement_keys", "igh_clones",
  "igh_grouping", "qc_audit", "rearrangement_groups", "record_qc"))

grouping_audit_fixture <- function(loci = QC_AUDIT_LOCI) {
  do.call(rbind, lapply(seq_along(loci), function(index) {
    locus <- loci[index]
    rows <- data.frame(
      library_id = "0460", locus = locus,
      sequence_id = paste0("assemble", rep(c(2 * index - 1, 2 * index), each = 2),
                           "_", rep(c(0, 1), 2)),
      consensus_count = c(100, 3, 5, 7), junction = rep("TGTGCGAAA", 4),
      junction_aa = "CAK", v_call = paste0(locus, "V1*01"),
      j_call = paste0(locus, "J1*01"), v_gene = paste0(locus, "V1"),
      j_gene = paste0(locus, "J1"), v_gene_status = "resolved", j_gene_status = "resolved",
      v_identity = 0.99, sequence_alignment = paste0(c("A", "C", "G", "T"), "ACGTACGTACGT"),
      sequence = paste0(c("A", "C", "G", "T"), "ACGTACGTACGT"),
      stringsAsFactors = FALSE
    )
    if (locus == "IGH") {
      rows$v_gene_status[1] <- "ambiguous"
      rows$v_gene[1] <- ""
    } else {
      rows$junction[1] <- "TGTNCGAAA"
    }
    rows
  }))
}

grouping_counts <- function(rows, stage, locus) {
  as.numeric(rows[rows$stage == stage & rows$locus == locus,
                  c("n_units_in", "n_units_out", "consensus_count_in", "consensus_count_out")])
}


test_that("grouping audits preserve BCR support and retain eligible TCR variants", {
  records <- grouping_audit_fixture()
  mapping <- prepare_grouping_fixture(records, "0460", 0)
  rows <- qc_audit_grouping_rows("0460", records, mapping, QC_AUDIT_LOCI)
  for (locus in QC_AUDIT_LOCI) {
    expect_identical(grouping_counts(rows, "grouping_eligible", locus),
      if (locus %in% c("IGH", "IGK", "IGL")) c(4, 2, 115, 12) else c(4, 3, 115, 15))
    if (locus %in% c("IGH", "IGK", "IGL")) {
      expect_identical(grouping_counts(rows, "assembly_backbone_collapse", locus), c(2, 1, 12, 12))
      expect_identical(grouping_counts(rows, "inference_representatives", locus), c(1, 1, 12, 12))
    } else {
      expect_equal(sum(rows$stage == "assembly_backbone_collapse" & rows$locus == locus), 0)
      expect_identical(grouping_counts(rows, "inference_representatives", locus), c(3, 3, 15, 15))
    }
  }
  expect_false("all" %in% rows$locus)
  path <- tempfile(fileext = ".tsv")
  on.exit(unlink(path))
  write_qc_audit(rows, path)
  back <- read_qc_audit(path)
  expect_identical(back[, QC_AUDIT_NUMERIC_COLUMNS], rows[, QC_AUDIT_NUMERIC_COLUMNS])
  expect_identical(back$library_id, rows$library_id)
})

test_that("an empty applicable branch records zero counts with explicit units", {
  records <- grouping_audit_fixture("IGH")[0, ]
  rows <- qc_audit_grouping_rows("0460", records, make_empty_grouping(), "IGH")
  expect_identical(nrow(rows), 4L)
  expect_true(all(rows$n_units_in == 0 & rows$n_units_out == 0))
  expect_true(all(rows$consensus_count_in == 0 & rows$consensus_count_out == 0))
  expect_identical(rows$output_unit, c("rearrangement", "backbone_representative", "inference_representative", "inference_representative"))
})

test_that("all-ineligible input retains its measured source population", {
  records <- grouping_audit_fixture("IGH")
  records$junction <- ""
  mapping <- prepare_igh_grouping(records, "0460", 0)
  rows <- qc_audit_grouping_rows("0460", records, mapping, "IGH")
  expect_identical(grouping_counts(rows, "grouping_eligible", "IGH"), c(4, 0, 115, 0))
  expect_identical(grouping_counts(rows, "assembly_backbone_collapse", "IGH"), c(0, 0, 0, 0))
})

test_that("missing links, wrong contexts and corrupt transferred support fail", {
  records <- grouping_audit_fixture("IGH")
  mapping <- prepare_igh_grouping(records, "0460", 0)
  expect_error(qc_audit_grouping_rows("0460", records, mapping[-1, ], "IGH"), "ordered source keys")
  bad <- mapping
  bad$sequence_id[2] <- bad$sequence_id[1]
  expect_error(qc_audit_grouping_rows("0460", records, bad, "IGH"), "ordered source keys")
  bad <- mapping
  bad$locus[1] <- "TRB"
  expect_error(qc_audit_grouping_rows("0460", records, bad, "IGH"), "disagrees")
  bad <- mapping
  bad$inference_consensus_count[4] <- 103
  expect_error(qc_audit_grouping_rows("0460", records, bad, "IGH"), "preserve eligible support")
  bad <- mapping
  bad$inference_consensus_count[1] <- 100
  expect_error(qc_audit_grouping_rows("0460", records, bad, "IGH"), "links or support")
})

test_that("canonical empty loci and native totals retain defined denominators", {
  rows <- qc_audit_gate_rows("NA", "retained_native", c(TRUE, FALSE, TRUE),
    c(4, 5, 7), c("IGH", "IGH", "TRB"), analysis = "common", stage_index = 1,
    locus_basis = "native", kind = "eligibility", support_conserved = FALSE)
  expect_identical(grouping_counts(rows, "retained_native", "IGH"), c(2, 1, 9, 4))
  expect_identical(grouping_counts(rows, "retained_native", "TRD"), c(0, 0, 0, 0))
  expect_identical(grouping_counts(rows, "retained_native", "all"), c(3, 2, 16, 11))
  expect_error(qc_audit_gate_rows("NA", "gate", NA, 1, "IGH",
    analysis = "common", stage_index = 1), "decision mask")
  expect_error(qc_audit_locus_rows("NA", "retained", 0, loci = "IGH",
    analysis = "common", stage_index = 1), "all canonical loci")
})

test_that("common-QC criteria count overlapping applied subsets", {
  flags <- data.frame(
    locus_match = c(TRUE, FALSE, TRUE, TRUE, TRUE, TRUE),
    d_locus_match = rep(TRUE, 6),
    low_n = c(TRUE, TRUE, FALSE, TRUE, TRUE, FALSE),
    productive = c(TRUE, TRUE, TRUE, FALSE, TRUE, FALSE),
    junction_mod3 = c(TRUE, TRUE, TRUE, TRUE, FALSE, TRUE),
    positive_support = rep(TRUE, 6)
  )
  rows <- qc_audit_qc_rows("G1", flags, airr_qc_pass(flags),
                          rep(1, 6), rep("IGH", 6), 0.1)
  expect_equal(nrow(rows), 7L * 8L)
  expect_identical(grouping_counts(rows, "qc_pass_total", "IGH"), c(6, 1, 6, 1))
  expect_identical(grouping_counts(rows, "qc_locus_match", "IGH"), c(6, 5, 6, 5))
  expect_identical(grouping_counts(rows, "qc_productive", "IGH"), c(6, 4, 6, 4))
  expect_true(all(rows$n_units_out[rows$locus == "TRD"] == 0))
  disabled <- qc_audit_qc_rows("G1", flags[1, ], TRUE, 1, "IGH", 1)
  expect_true(all(disabled$decision_status[disabled$stage == "qc_low_n"] == "disabled"))
  path <- tempfile()
  on.exit(unlink(path))
  write_qc_audit(rows, path)
  back <- read_qc_audit(path)
  expect_identical(names(back), QC_AUDIT_COLUMNS)
  expect_identical(back[, QC_AUDIT_NUMERIC_COLUMNS], rows[, QC_AUDIT_NUMERIC_COLUMNS])
})

test_that("audit validation rejects invalid counts and preserves missing support", {
  expect_error(qc_audit_rows("G1", "stage", 2, 1, analysis = "common", stage_index = 1),
               "checkpoint")
  expect_error(qc_audit_rows("G1", "stage", 1, 2, analysis = "common", stage_index = 1,
                            kind = "eligibility"), "subset")
  expect_error(qc_audit_rows("G1", "stage", 1, analysis = "common", stage_index = 1,
                            locus = "XXX"), "locus")
  expect_error(qc_audit_rows("G1", "stage", 0.5, analysis = "common", stage_index = 1),
               "Non-integer")
  path <- tempfile()
  on.exit(unlink(path))
  big <- qc_audit_rows("G1", "stage", 2e6, 1e5, 5613344, 1e6,
    analysis = "common", stage_index = 1, kind = "eligibility", support_conserved = FALSE)
  write_qc_audit(big, path)
  lines <- readLines(path)
  expect_false(any(grepl("e[+-]", lines)))
  expect_true(any(grepl("\t2000000\t", lines)))
  expect_true(any(grepl("\t5613344\t", lines)))
  null <- qc_audit_rows("G1", "native_calls", 42, analysis = "assembly",
                        stage_index = 1, input_unit = "cdr3_call", locus_basis = "native")
  write_qc_audit(null, path)
  expect_true(is.na(read_qc_audit(path)$consensus_count_in))
})
