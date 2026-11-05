#!/usr/bin/env Rscript
# Keep records that pass the shared post-annotation QC criteria at every locus.
# airr_qc_flags() in workflow/lib/record_qc.R tests locus agreement, productive
# rearrangement, junction length divisible by three, positive read support and
# aligned N fraction. params.airr_qc.max_alignment_n_fraction supplies the N
# ceiling. Each mask uses the same input rows, so failure counts can overlap.
#
# This script writes retained AIRR rows and the QC audit. A common-QC pass does
# not establish eligibility for every analysis: grouping resolves V/J and
# junction requirements, qualify_isotype_calls.R masks unsupported C calls and
# SHM applies coverage to each measured region. Those stages retain their own
# eligibility decisions rather than adding them to this shared filter.

source(file.path(dirname(sub("^--file=", "",
  grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1])),
  "..", "lib", "runtime.R"))
load_pipeline_modules(c(
  "table_io", "airr_io", "gene_calls", "qc_audit", "record_qc"))
stage_log <- start_stage_log("filter_airr")

suppressPackageStartupMessages({
  library(optparse)
})

script_options <- parse_script_arguments(list(
  make_option("--input_file", type = "character",
              help = "All-locus reannotated AIRR table for one library. Required."),
  make_option("--output_file", type = "character",
              help = "Path to write the QC-pass AIRR table to. Required."),
  make_option("--library_id", type = "character",
              help = "Library identifier every input row must carry. Required."),
  make_option("--max_n_fraction", type = "double",
              help = paste("Largest fraction of N bases among the non-gap",
                           "sequence_alignment positions, from 0 to 1; 1",
                           "removes the ceiling. Required.")),
  make_option("--qc_audit_file", type = "character",
              help = "Path to write the common-QC qc_audit shard to. Required.")
))
input_file <- script_options$input_file
output_file <- script_options$output_file
library_id <- script_options$library_id
max_n_fraction <- script_options$max_n_fraction
qc_audit_file <- script_options$qc_audit_file

db <- read_pipeline_table(input_file)
log_info("Loaded %d all-locus MakeDb-pass rows from %s", nrow(db), input_file)

require_columns(db, c("sequence_id", "library_id", "locus", "v_call", "j_call",
                      "d_call", "v_gene_status", "j_gene_status",
                      "junction", "sequence_alignment", "productive",
                      "junction_length", "consensus_count"),
                "reannotated AIRR")
if (anyDuplicated(db$sequence_id)) stop("sequence_id is not unique")
if (any(as.character(db$library_id) != library_id)) {
  stop("reannotated AIRR contains a library_id other than ", library_id)
}

support <- as_read_support(db$consensus_count, "consensus_count")
flags <- airr_qc_flags(db, max_n_fraction)
pass <- airr_qc_pass(flags)
for (criterion in names(flags)) {
  log_info("QC %-18s pass=%d fail=%d", criterion,
           sum(flags[[criterion]]), sum(!flags[[criterion]]))
}
log_info("Common QC kept %d and dropped %d record(s)", sum(pass), sum(!pass))

qc_audit <- qc_audit_qc_rows(library_id, flags, pass, support, db$locus, max_n_fraction)

# Every surviving record has a locus that is one of the seven, because the
# locus_match criterion above only passes a record whose V and J genes name the
# locus it was assigned, and only those seven can be named. So the ordering below
# can index straight into that list.
db <- db[pass, , drop = FALSE]
# The retained table is written in a fixed locus-then-sequence_id order so that
# two runs of one library produce the same rows in the same order.
# set_pipeline_collation() fixes the C collation this sort uses.
db <- db[order(match(db$locus, QC_AUDIT_LOCI), db$sequence_id), , drop = FALSE]

write_qc_audit(qc_audit, qc_audit_file)

write_airr_table(db, output_file)
log_info("Wrote %d common-QC all-locus records to %s", nrow(db), output_file)

finish_stage_log(stage_log)
