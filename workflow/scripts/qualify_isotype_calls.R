#!/usr/bin/env Rscript
# Qualify constant calls for every common-QC IGH record without deleting rows.
# The keyed sidecar feeds isotype, SHM strata, selection and library AIRR export.
# build_qualified_isotype_calls() in workflow/lib/isotype_calls.R decides each
# call, and resolve_isotype() in workflow/lib/isotype_resolution.R names its
# class and subclass. This script writes the sidecar and one audit row per
# criterion.

suppressPackageStartupMessages(library(optparse))

source(file.path(dirname(sub("^--file=", "",
  grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1])),
  "..", "lib", "runtime.R"))
load_pipeline_modules(c(
  "table_io", "gene_calls", "isotype_calls", "isotype_resolution", "qc_audit"))
stage_log <- start_stage_log("qualify_isotype_calls")

script_options <- parse_script_arguments(list(
  make_option("--input_file", type = "character",
              help = "Common QC-pass AIRR table for one library. Required."),
  make_option("--library_id", type = "character",
              help = "Library the records belong to. Required."),
  make_option("--min_c_positions", type = "integer",
              help = paste("Fewest informative C alignment positions for a callable",
                           "constant call; zero disables the length floor. Required.")),
  make_option("--qc_audit_file", type = "character",
              help = "Path to write the qualification audit to. Required."),
  make_option("--output_file", type = "character",
              help = "Path to write the qualified call sidecar to. Required.")
))
input_file <- script_options$input_file
output_file <- script_options$output_file
sample_label <- script_options$library_id
min_c_positions <- script_options$min_c_positions

records <- read_pipeline_table(input_file)
calls <- build_qualified_isotype_calls(records, sample_label, min_c_positions)
# The isotype tables and clone members resolve the same qualified call with the
# same function, so the AIRR export names the isotype they count. A masked,
# missing or unresolved call leaves the column empty.
calls$isotype_class <- resolve_isotype(calls$qualified_c_call, "class")
calls$isotype_subclass <- resolve_isotype(calls$qualified_c_call, "subclass")
igh <- records[records$locus == "IGH", , drop = FALSE]
support <- as_read_support(igh$consensus_count[match(calls$sequence_id, igh$sequence_id)])
# One audit row per qualification criterion. The criteria share the same
# input population and can fail together, so their counts are not sequential
# losses. support_conserved is FALSE because a criterion selects a subset:
# the default would assert equal support only because support_in was supplied.
# A zero threshold is reported as disabled because the length floor is what the
# setting controls; presence, locus compatibility and a finite alignment
# measurement still apply, so records without a C alignment still fail.
audit <- lapply(c("present", "locus_match", "enough_evidence", "callable", "candidates_complete"), function(criterion) {
  qc_audit_rows(sample_label, paste0("c_call_", criterion), nrow(calls),
    sum(calls[[paste0("c_call_", criterion)]]), sum(support),
    sum(support[calls[[paste0("c_call_", criterion)]]]),
    analysis = "isotype", stage_index = 40 + match(criterion,
      c("present", "locus_match", "enough_evidence", "callable", "candidates_complete")),
    kind = if (criterion == "callable") "eligibility" else "criterion",
    decision_status = if (criterion == "enough_evidence" && min_c_positions == 0)
      "disabled" else "applied", support_conserved = FALSE)
})
write_qc_audit(do.call(rbind, audit), script_options$qc_audit_file)
log_info(
  "C-call QC present=%d callable=%d low_or_missing_evidence=%d locus_mismatch=%d",
  sum(calls$c_call_present), sum(calls$c_call_callable),
  sum(calls$c_call_present & !calls$c_call_enough_evidence),
  sum(calls$c_call_present & !calls$c_call_locus_match)
)
write_pipeline_table(calls, output_file)
log_info("Wrote %d qualified IGH call rows to %s", nrow(calls), output_file)
finish_stage_log(stage_log)
