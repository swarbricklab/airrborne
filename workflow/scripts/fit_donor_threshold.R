#!/usr/bin/env Rscript
# Fit one donor's IGH junction distances within one scope using SHazaM.
# fit_donor_distance_threshold() in donor_distance_threshold.R owns the distance
# and fit choices. This script reads the prepared representatives and writes the
# one-row fit result and diagnostic histogram. Unusable distances or a failed
# fit leave the threshold missing while retaining the available distance evidence,
# and fit_status records which reason applied.
# The rule runs even with numeric params.lineage.clonal_threshold: fitting records
# evidence; resolve_scope_clone_threshold.R separately chooses the cutoff used.

suppressPackageStartupMessages({
  library(optparse)
  library(shazam)
})

source(file.path(dirname(sub("^--file=", "",
  grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1])),
  "..", "lib", "runtime.R"))
load_pipeline_modules(c(
  "table_io", "donor_distance_threshold", "gene_calls",
  "assembly_backbone", "igh_clones", "igh_grouping", "rearrangement_keys"))
stage_log <- start_stage_log("fit_donor_threshold")

script_options <- parse_script_arguments(list(
  make_option("--input_file", type = "character",
              help = "Donor-pooled QC-pass IGH table for one scope. Required."),
  make_option("--output_file", type = "character",
              help = "Path to write the one-row donor fit to. Required."),
  make_option("--diagnostic_file", type = "character",
              help = paste("Path to write the 100-bin distance histogram and",
                           "fitted mixture to. Required.")),
  make_option("--scope_id", type = "character",
              help = "Analysis scope this donor is being fitted within. Required."),
  make_option("--donor_id", type = "character",
              help = "Donor whose libraries were pooled. Required."),
  make_option("--minimum_n", type = "integer",
              help = paste("Fewest finite nearest-neighbour distances needed",
                           "before a threshold is fitted; at least two.",
                           "Required.")),
  make_option("--fit_edge", type = "double",
              help = "GMM initialisation scan fraction, strictly between 0 and 1. Required."),
  make_option("--fit_specificity", type = "double",
              help = "Specificity for SHazaM's user cutoff, between 0 and 1. Required."),
  make_option("--fit_subsample", type = "character",
              help = "Maximum sampled distances, or null to use all distances. Required."),
  make_option("--grouping_file", type = "character",
              help = "Prepared donor grouping beside the source pool. Required.")
))
input_file <- script_options$input_file
output_file <- script_options$output_file
diagnostic_file <- script_options$diagnostic_file
scope_id <- script_options$scope_id
donor_id <- script_options$donor_id
minimum_n <- script_options$minimum_n
fit_subsample <- if (script_options$fit_subsample == "null") NULL else
  as.numeric(script_options$fit_subsample)

db <- read_pipeline_table(
  input_file,
  required_columns = c("sequence_id", "junction", "v_call", "j_call",
                       "v_gene_status", "j_gene_status", "locus",
                       "consensus_count"),
  context = "Pooled IGH"
)

# Use the same prepared sequences as assignment. Read-support totals are carried
# for consistency but do not weight SHazaM's nearest-neighbour distances.
mapping <- read_pipeline_table(script_options$grouping_file, GROUPING_COLUMNS)
fit_db <- prepared_igh_representatives(db, mapping, inference_support = TRUE)
# The representatives entering the fit, and the libraries the grouping links
# them to.
n_members <- nrow(fit_db)
n_libraries <- length(unique(mapping$library_id[mapping$sequence_id %in% fit_db$sequence_id]))

if (!n_members) {
  write_pipeline_table(build_donor_threshold_row(scope_id, donor_id), output_file)
  write_pipeline_table(build_distance_histogram_diagnostics(scope_id, donor_id),
                       diagnostic_file)
  log_info("No eligible IGH junctions; donor threshold is not fitted")
  finish_stage_log(stage_log)
  quit(save = "no", status = 0)
}

log_info(
  paste0("Fitting %s/%s on %d records: distToNearest first=FALSE; ",
         "findThreshold gmm gamma-norm edge=%s spc=%s subsample=%s"),
  scope_id, donor_id, nrow(fit_db), script_options$fit_edge,
  script_options$fit_specificity, script_options$fit_subsample)

fit <- fit_donor_distance_threshold(fit_db, minimum_n, script_options$fit_edge,
                                    script_options$fit_specificity, fit_subsample)
n_distances <- length(fit$distances)
if (!is.null(fit$error)) {
  log_warn("Threshold fit failed: %s", fit$error)
}
if (fit$fit_status %in% c("too_few_distances", "constant_distances")) {
  log_info("Donor threshold not fitted (%s): %d finite distances (minimum %d)",
           fit$fit_status, n_distances, minimum_n)
}

write_pipeline_table(
  build_donor_threshold_row(scope_id, donor_id, n_libraries, n_members,
                            n_distances, fit$threshold, fit$fit_status,
                            fit$n_distances_fitted, fit$fit),
  output_file)
write_pipeline_table(
  build_distance_histogram_diagnostics(scope_id, donor_id, fit$distances,
                                       fit$fit),
  diagnostic_file)
log_info("Donor threshold: %s (%s; %d finite distances from %d members in %d libraries)",
         if (identical(fit$fit_status, "fitted")) sprintf("%.4f", fit$threshold) else "none",
         fit$fit_status, n_distances, n_members, n_libraries)

finish_stage_log(stage_log)
