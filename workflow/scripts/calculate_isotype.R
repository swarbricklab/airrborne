#!/usr/bin/env Rscript
# Calculate library IGH class and subclass composition from common-QC AIRR
# and qualified constant calls. Write complete ordered strata and their audit.
# build_isotype_metrics() in workflow/lib/isotype_metrics.R supplies the
# composition, the one-row summary and the callable populations the audit counts.
# Fractions use recovered read support callable at their own resolution, class or
# subclass. Unresolved calls keep their support in the total. With no IGH records,
# the ordered strata and summary retain zero counts and missing fractions.

suppressPackageStartupMessages(library(optparse))

source(file.path(dirname(sub("^--file=", "",
  grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1])),
  "..", "lib", "runtime.R"))
load_pipeline_modules(c(
  "table_io", "gene_calls", "isotype_calls", "isotype_resolution",
  "isotype_metrics", "qc_audit"))
stage_log <- start_stage_log("calculate_isotype")

script_options <- parse_script_arguments(list(
  make_option("--input_file", type = "character",
              help = "Common QC-pass AIRR table for one library. Required."),
  make_option("--isotype_calls_file", type = "character",
              help = "Qualified constant-call sidecar for the same library. Required."),
  make_option("--library_id", type = "character",
              help = "Library both tables belong to. Required."),
  make_option("--qc_audit_file", type = "character",
              help = "Path to write the callable-support audit to. Required."),
  make_option("--isotype_file", type = "character",
              help = "Path to write the class and subclass composition to. Required."),
  make_option("--summary_file", type = "character",
              help = "Path to write the one-row composition summary to. Required.")
))
input_file   <- script_options$input_file
calls_file   <- script_options$isotype_calls_file
output_file  <- script_options$isotype_file
sample_label <- script_options$library_id

df <- read_pipeline_table(input_file)
calls <- read_pipeline_table(calls_file)
log_info("Loaded %d rows from %s", nrow(df), input_file)
log_info("Loaded %d qualified IGH call rows from %s", nrow(calls), calls_file)

require_columns(df, c("sequence_id", "library_id", "locus",
                      "consensus_count"), "Common-QC AIRR")
require_columns(calls, c("library_id", "sequence_id", "qualified_c_call"),
                "Qualified isotype calls")
if (nrow(df) && any(as.character(df$library_id) != sample_label)) {
  stop("Common-QC AIRR contains a library_id other than ", sample_label)
}
if (nrow(calls) && any(as.character(calls$library_id) != sample_label)) {
  stop("Qualified isotype calls contain a library_id other than ", sample_label)
}
if (nrow(calls) && any(is_airr_null(calls$sequence_id))) {
  stop("Qualified isotype calls contain a missing sequence_id")
}

df <- df[as.character(df$locus) == "IGH", , drop = FALSE]
if (nrow(df) && any(is_airr_null(df$sequence_id))) {
  stop("Common-QC IGH AIRR contains a missing sequence_id")
}
df <- join_qualified_isotype_calls(df, calls)
df$consensus_count <- as_read_support(
  df$consensus_count, "consensus_count")
if (any(df$consensus_count <= 0)) {
  stop("common QC input contains non-positive read support")
}

metrics <- build_isotype_metrics(df)
audit <- lapply(names(metrics$resolved), function(resolution) {
  keep <- !is.na(metrics$resolved[[resolution]])
  qc_audit_rows(sample_label, paste0("isotype_", resolution, "_callable"),
    nrow(df), sum(keep), sum(df$consensus_count), sum(df$consensus_count[keep]),
    analysis = "isotype", stage_index = 45 + match(resolution, names(metrics$resolved)),
    kind = "eligibility", support_conserved = FALSE)
})
write_qc_audit(do.call(rbind, audit), script_options$qc_audit_file)
write_pipeline_table(data.frame(library_id = sample_label, metrics$composition),
                     output_file)
write_pipeline_table(data.frame(library_id = sample_label, metrics$summary),
                     script_options$summary_file)
log_info("Wrote class and subclass composition with their callable denominators")
finish_stage_log(stage_log)
