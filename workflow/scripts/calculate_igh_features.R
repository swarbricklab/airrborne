#!/usr/bin/env Rscript
# Measure mutation counts, opportunities and consensus-vote support for each
# common-QC IGH record, whole V and by IMGT sub-region.
# build_igh_v_region_feature_table() in igh_features.R makes the measurements;
# params.shm.min_position_votes defines a supported position. An absent alignment
# pair leaves measurements missing. No IGH records produces a header-only table.
# Qualified constant calls supply isotype context through a keyed join.
# Donor summaries apply regional coverage and choose each member's evidence source.

source(file.path(dirname(sub("^--file=", "",
  grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1])),
  "..", "lib", "runtime.R"))
load_pipeline_modules(c(
  "table_io", "gene_calls", "igh_features", "isotype_calls",
  "qc_audit", "shm_regions", "trust4_consensus_votes"))
stage_log <- start_stage_log("calculate_igh_features")

suppressPackageStartupMessages({
  library(optparse)
  library(shazam)
})

script_options <- parse_script_arguments(list(
  make_option("--input_file", type = "character",
              help = "Common-QC all-locus AIRR table for one library. Required."),
  make_option("--votes_file", type = "character",
              help = paste("IGH consensus-vote table for the same library,",
                           "carrying the per-IMGT-position read votes. Required.")),
  make_option("--isotype_calls_file", type = "character",
              help = "Qualified IGH isotype-call sidecar for the same library. Required."),
  make_option("--min_position_votes", type = "integer",
              help = "Fewest agreeing votes for a supported IMGT V position. Required."),
  make_option("--qc_audit_file", type = "character", help = "SHM measurement counts. Required."),
  make_option("--output_file", type = "character",
              help = "Path to write the keyed IGH feature table to. Required."),
  make_option("--library_id", type = "character",
              help = "Library identifier every IGH row must carry. Required.")
))
input_file <- script_options$input_file
votes_file <- script_options$votes_file
isotype_calls_file <- script_options$isotype_calls_file
output_file <- script_options$output_file
library_id <- script_options$library_id

db <- read_pipeline_table(input_file)
require_columns(db, c("sequence_id", "library_id", "locus", "v_call", "sequence_alignment",
                      "germline_alignment_d_mask"), "common AIRR")
db <- db[as.character(db$locus) == "IGH", , drop = FALSE]
if (nrow(db) && anyDuplicated(db$sequence_id)) stop("IGH sequence_id is not unique")
if (nrow(db) && any(as.character(db$library_id) != library_id)) {
  stop("common AIRR contains an IGH library_id other than ", library_id)
}
log_info("Read %d QC-pass IGH record(s) from %s", nrow(db), input_file)

votes <- read_pipeline_table(votes_file)
log_info("Read %d IGH consensus-vote row(s) from %s", nrow(votes), votes_file)

calls <- read_pipeline_table(isotype_calls_file,
                              required_columns = c("library_id", "sequence_id", "qualified_c_call"))
if (nrow(calls) && any(as.character(calls$library_id) != library_id)) {
  stop("Qualified isotype calls contain a library_id other than ", library_id)
}
db <- join_qualified_isotype_calls(db, calls)
log_info("Joined %d qualified IGH calls from %s", nrow(calls), isotype_calls_file)

features <- build_igh_v_region_feature_table(db, votes,
                                            script_options$min_position_votes)

# The measurement's own results supply record and position counts. Read support
# is undefined for positions because one record spans several positions.
audit <- list()
for (region in V_REGION_LEVELS) {
  informative <- features[[paste0("informative_", region)]]
  supported <- features[[paste0("informative_supported_", region)]]
  measured <- !is.na(informative) & informative > 0
  audit[[length(audit) + 1L]] <- qc_audit_rows(
    library_id, "shm_measured_records", nrow(db), sum(measured),
    sum(as_read_support(db$consensus_count)), sum(as_read_support(db$consensus_count)[measured]),
    analysis = "shm", stage_index = 50, v_region = region,
    kind = "eligibility", support_conserved = FALSE)
  audit[[length(audit) + 1L]] <- qc_audit_rows(
    library_id, "shm_supported_positions", sum(informative, na.rm = TRUE),
    sum(supported, na.rm = TRUE), analysis = "shm", stage_index = 51, v_region = region,
    input_unit = "record_position", kind = "eligibility", support_conserved = NA,
    decision_status = if (script_options$min_position_votes == 0) "disabled" else "applied")
}
write_qc_audit(do.call(rbind, audit), script_options$qc_audit_file)

# A record with no alignment pair to compare has nothing to measure, which is why
# an absent informative_v is how that case is counted here. Nothing else is
# withheld: the coverage a given analysis requires is applied where that analysis
# forms its released number.
log_info(
  "IGH V-region mutation rate calculated for %d record(s), no alignment pair=%d",
  sum(!is.na(features$shm_frequency)), sum(is.na(features$informative_v))
)
log_info(
  "V-call germline ambiguity: %s",
  paste(sprintf("%s=%d", names(table(features$v_call_status)),
                as.integer(table(features$v_call_status))), collapse = " ")
)
callability <- summarise_igh_feature_callability(features)
for (index in seq_len(nrow(callability))) {
  log_info(
    paste0("%s: median informative=%s (%s read-supported) of %s ",
           "germline-occupied positions, median majority votes=%s, rate called ",
           "for %d record(s)"),
    callability$region[index],
    callability$median_informative[index],
    callability$median_informative_supported[index],
    callability$median_germline_occupied[index],
    callability$median_trust4_consensus_base_votes[index],
    callability$n_rate_called[index]
  )
}

features <- features[order(features$sequence_id), igh_feature_columns(),
                     drop = FALSE]
# Not write_airr_table(): these are derived rates and counts, not AIRR fields.
write_pipeline_table(features, output_file)
log_info("Wrote %d keyed IGH feature rows to %s", nrow(features), output_file)

finish_stage_log(stage_log)
