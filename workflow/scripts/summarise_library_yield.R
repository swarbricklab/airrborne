#!/usr/bin/env Rscript
# Sum TRUST4-reported CDR3 read support by canonical locus before imputed-CDR3
# removal, reannotation and common QC, then divide each locus sum by the
# library's input read pairs to report yield per million. The row also carries
# the library's input depth. Every library contributes one row, including
# all-zero sums for a header-only AIRR input. A library with no input read pairs
# has an unavailable yield rather than a zero one. build_repertoire_support_row()
# in repertoire_support.R constructs the totals; add_repertoire_yield() joins the
# depth and applies the input-read-pair denominator.

suppressPackageStartupMessages(library(optparse))

source(file.path(dirname(sub("^--file=", "",
  grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1])),
  "..", "lib", "runtime.R"))
load_pipeline_modules(c("table_io", "repertoire_support"))
stage_log <- start_stage_log("summarise_library_yield")

script_options <- parse_script_arguments(list(
  make_option("--airr_file", type = "character",
              help = "TRUST4 AIRR table for one library. Required."),
  make_option("--depth_file", type = "character",
              help = "One-row input depth table for the same library. Required."),
  make_option("--library_id", type = "character",
              help = "Library both tables belong to. Required."),
  make_option("--output_file", type = "character",
              help = "Path to write the one-row yield table to. Required.")
))
input_file   <- script_options$airr_file
output_file  <- script_options$output_file
sample_label <- script_options$library_id

df <- read_pipeline_table(
  input_file,
  required_columns = c("locus", "consensus_count"),
  context = "TRUST4 AIRR table")
log_info("Loaded %d rows, %d cols from %s", nrow(df), ncol(df), input_file)

consensus_count <- as_read_support(df$consensus_count, "consensus_count")

# A contig TRUST4 could not place on one of the seven canonical loci carries
# support that belongs to no class total, so it is reported and left out.
unexpected <- setdiff(unique(df$locus), SUPPORT_LOCI)
if (length(unexpected)) {
  log_warn("Excluded %.0f support on unexpected/blank loci from totals (values: %s)",
           sum(consensus_count[df$locus %in% unexpected]),
           paste(head(unexpected, 5), collapse = ", "))
}

support_by_locus <- vapply(SUPPORT_LOCI,
                           function(locus) sum(consensus_count[df$locus == locus]),
                           numeric(1))
log_info("Reported CDR3 support by locus: %s",
         paste(sprintf("%s=%.0f", SUPPORT_LOCI, support_by_locus),
               collapse = ", "))

result <- build_repertoire_support_row(sample_label, support_by_locus)
depth <- read_pipeline_table(script_options$depth_file,
                             required_columns = c("library_id", "n_fastq_file_pairs",
                                                  "input_read_pairs", "total_bases"))
if (nrow(depth) != 1L || depth$library_id != sample_label) {
  stop("Input depth must identify this library exactly once")
}
result <- add_repertoire_yield(result, depth)
write_pipeline_table(result, output_file)
log_info("Wrote %s: total=%.0f tcr=%.0f bcr=%.0f", output_file,
         result$trust4_consensus_count, result$tcr_trust4_consensus_count,
         result$bcr_trust4_consensus_count)

finish_stage_log(stage_log)
