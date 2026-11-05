# Run native parsing while preserving source identifiers and observed CDR3s.
# Run: python test/run_tests.py --component trust4 --level integration

load_test_modules(c())

work_dir <- tempfile("trust4-parsing-")
dir.create(work_dir)

# Write a headerless 13-column TRUST4 cdr3.out fixture. Only columns 1-2
# (consensus_id, index -> sequence_id) and 10 (CDR3 score) are read.
write_cdr3 <- function(path, ids, scores) {
  cid <- sub("_[^_]*$", "", ids)
  idx <- sub("^.*_", "", ids)
  df <- data.frame(
    cid = cid, idx = idx, v = "IGHV1-8*01", d = ".", j = "IGHJ5*02", c = "IGHM",
    cdr1 = "AAA", cdr2 = "CCC", cdr3 = "TGTGCC", score = sprintf("%.2f", scores),
    reads = "1", germline_similarity = "95.0", complete_vdj = "1",
    stringsAsFactors = FALSE
  )
  write.table(df, path, sep = "\t", row.names = FALSE, col.names = FALSE,
              quote = FALSE)
}

run_parser <- function(input, output, cdr3, audit = tempfile("audit-", work_dir)) {
  run <- run_pipeline_script("parse_trust4_airr.R", c(
    paste0("--input_file=", input),
    paste0("--output_file=", output),
    paste0("--cdr3_file=", cdr3),
    paste0("--qc_audit_file=", audit),
    "--library_id=GUSS_TEST"
  ))
  c(run, list(audit = audit))
}

# One (stage, locus) cell of a written QC audit shard.
audit_cell <- function(path, stage, locus, column) {
  rows <- read.delim(path, stringsAsFactors = FALSE, check.names = FALSE)
  as.numeric(rows[[column]][rows$stage == stage & rows$locus == locus])
}

# --- Preserve exact sequences and source identities --------------------------
input_file <- file.path(work_dir, "native.tsv")
output_file <- file.path(work_dir, "parsed.tsv")
cdr3_file <- file.path(work_dir, "cdr3.out")
native <- data.frame(
  sequence_id = c("assemble0_0", "assemble0_1", "assemble9_3"),
  sequence = c("AAACCC", "AAAGCC", "AAACCC"),
  locus = "IGH",
  consensus_count = c(10, 2, 5),
  stringsAsFactors = FALSE
)
write.table(native, input_file, sep = "\t", row.names = FALSE, quote = FALSE)
write_cdr3(cdr3_file, native$sequence_id, rep(1.00, nrow(native)))

res <- run_parser(input_file, output_file, cdr3_file)
expect_script_ok(res)

test_that("identical sequences retain all original identifiers and support", {
  parsed <- read.delim(output_file, stringsAsFactors = FALSE, check.names = FALSE)
  expect_identical(parsed$sequence_id, native$sequence_id)
  expect_identical(parsed$sequence, native$sequence)
  expect_identical(as.numeric(parsed$consensus_count), c(10, 2, 5))
})

test_that("the parser audit preserves sources and leaves exact merging to grouping", {
  audit <- read.delim(res$audit, stringsAsFactors = FALSE)
  expect_false("exact_sequence_collapse" %in% audit$stage)
  expect_identical(audit_cell(res$audit, "igblast_locus_routing", "all", "n_units_out"), 3)
  expect_identical(audit_cell(res$audit, "igblast_locus_routing", "all", "consensus_count_out"), 17)
  expect_identical(
    audit_cell(res$audit, "trust4_reported_airr", "IGH", "n_units_in"), 3)
  expect_identical(
    audit_cell(res$audit, "trust4_cdr3_calls", "all", "n_units_in"), 3)
  expect_identical(
    audit_cell(res$audit, "trust4_contigs_with_cdr3_call", "all", "n_units_in"), 2)
  expect_identical(
    audit_cell(res$audit, "trust4_reported_airr", "TRB", "n_units_in"), 0)
})

# --- Imputed CDR3 (score 0.01) is dropped before annotation --------------------
# assemble5_0 is imputed and must be removed; a partial-motif score (0.83) and a
# fully observed score (1.00) must be retained.
imp_input <- file.path(work_dir, "native_imp.tsv")
imp_out <- file.path(work_dir, "parsed_imp.tsv")
imp_cdr3 <- file.path(work_dir, "cdr3_imp.out")
native_imp <- data.frame(
  sequence_id = c("assemble0_0", "assemble5_0", "assemble5_1"),
  sequence = c("AAACCC", "TTTGGG", "TTTAGG"),
  locus = "IGH",
  consensus_count = c(10, 7, 3),
  stringsAsFactors = FALSE
)
write.table(native_imp, imp_input, sep = "\t", row.names = FALSE, quote = FALSE)
write_cdr3(imp_cdr3, native_imp$sequence_id, c(1.00, 0.01, 0.83))

imp <- run_parser(imp_input, imp_out, imp_cdr3)
expect_script_ok(imp)

test_that("an imputed CDR3 is dropped before annotation", {
  parsed_imp <- read.delim(imp_out, stringsAsFactors = FALSE, check.names = FALSE)
  # imputed assemble5_0 dropped; the two non-imputed distinct sequences remain
  expect_identical(parsed_imp$sequence_id, c("assemble0_0", "assemble5_1"))
  expect_identical(as.numeric(parsed_imp$consensus_count), c(10, 3))
  expect_false("assemble5_0" %in% parsed_imp$sequence_id)
  # The audit carries the same drop as a released number, with its denominator.
  expect_identical(
    audit_cell(imp$audit, "imputed_cdr3_filter", "IGH", "n_units_in"), 3)
  expect_identical(
    audit_cell(imp$audit, "imputed_cdr3_filter", "IGH", "n_units_out"), 2)
  expect_identical(
    audit_cell(imp$audit, "imputed_cdr3_filter", "IGH", "consensus_count_in"), 20)
  expect_identical(
    audit_cell(imp$audit, "imputed_cdr3_filter", "IGH", "consensus_count_out"), 13)
})

# --- Conflicting locus for one exact sequence must fail ----------------------
conflict_file <- file.path(work_dir, "conflict.tsv")
conflict_out <- file.path(work_dir, "conflict_out.tsv")
conflict <- native[c(1, 3), ]
conflict$locus <- c("IGH", "TRB")
write.table(conflict, conflict_file, sep = "\t", row.names = FALSE, quote = FALSE)
conflict <- run_parser(conflict_file, conflict_out, cdr3_file)

test_that("a conflicting locus for one exact sequence stops the stage", {
  expect_true(conflict$status != 0L)
  # The stage refuses rather than picking a locus, so nothing is released: the
  # parsed table the next rule would read is never written.
  expect_false(file.exists(conflict_out))
})
