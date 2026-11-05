# Choose the IGH distance cutoff shared by separately clustered donors in a scope.
# choose_scope_clone_threshold() in scope_clone_threshold.R chooses a configured
# number, the median usable donor fit or the fallback. This script checks donor
# coverage, saves all calibration evidence and warns on an automatic fallback.
# A configured number changes the operative cutoff without discarding donor fits.
# It does not assign clones.

source(file.path(dirname(sub("^--file=", "",
  grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1])),
  "..", "lib", "runtime.R"))
load_pipeline_modules(c(
  "donor_distance_threshold", "scope_clone_threshold", "table_io"))
stage_log <- start_stage_log("resolve_scope_clone_threshold")

suppressPackageStartupMessages(library(optparse))

script_options <- parse_script_arguments(list(
  make_option("--output_file", type = "character",
              help = paste("Path to write the one-row operative scope",
                           "threshold to. Required.")),
  make_option("--donor_output_file", type = "character",
              help = paste("Path to write the retained per-donor calibration",
                           "table to. Required.")),
  make_option("--diagnostics_output_file", type = "character",
              help = paste("Path to write the pooled per-donor distance",
                           "histograms to. Required.")),
  make_option("--scope_id", type = "character",
              help = "Analysis scope every donor fit must belong to. Required."),
  make_option("--configured_threshold", type = "character",
              help = paste("'auto' to take the median of the usable donor",
                           "fits, or a number from 0 to 1 to fix the cutoff.",
                           "Required.")),
  make_option("--donor_files", type = "character",
              help = paste("Comma-separated one-row donor fit tables, one per",
                           "donor in the scope. Required.")),
  make_option("--diagnostic_files", type = "character",
              help = paste("Comma-separated donor distance-diagnostic tables,",
                           "covering the same donors. Required."))
))
output_file <- script_options$output_file
donor_output_file <- script_options$donor_output_file
diagnostics_output_file <- script_options$diagnostics_output_file
scope_id <- script_options$scope_id
donor_files <- split_path_list(script_options$donor_files)
diagnostic_files <- split_path_list(script_options$diagnostic_files)
configured_threshold <- parse_configured_threshold(
  script_options$configured_threshold)

donors <- if (length(donor_files) > 0L) {
  do.call(rbind, lapply(donor_files, read_donor_threshold_row))
} else {
  make_empty_frame(DONOR_THRESHOLD_COLUMNS)
}
diagnostics <- if (length(diagnostic_files) > 0L) {
  do.call(rbind, lapply(diagnostic_files, read_donor_distance_diagnostics))
} else {
  data.frame(
    scope_id = character(0), donor_id = character(0),
    distance_bin_start = numeric(0), distance_bin_end = numeric(0),
    distance_bin_midpoint = numeric(0), n_distances_in_bin = integer(0),
    observed_density = numeric(0), near_component_density = numeric(0),
    far_component_density = numeric(0), stringsAsFactors = FALSE
  )
}
validate_donor_fit_consistency(donors, diagnostics, scope_id)

if (!all(donors$fit_status %in% DONOR_FIT_STATUSES)) {
  stop("donor threshold table contains a fit_status outside ",
       paste(DONOR_FIT_STATUSES, collapse = ", "))
}
donors$threshold <- suppressWarnings(as.numeric(donors$threshold))
for (column in c("n_libraries", "n_members", "n_distances", "n_distances_fitted")) {
  donors[[column]] <- suppressWarnings(as.integer(donors[[column]]))
}
for (column in names(DONOR_FIT_PARAMETERS)) {
  donors[[column]] <- suppressWarnings(as.numeric(donors[[column]]))
}

# The released table holds every donor column, naming the fit's threshold.
donor_calibration <- donors[, DONOR_THRESHOLD_COLUMNS, drop = FALSE]
names(donor_calibration)[names(donor_calibration) == "threshold"] <- "fitted_threshold"
if (nrow(donor_calibration)) {
  donor_calibration <- donor_calibration[
    order(donor_calibration$donor_id), , drop = FALSE
  ]
  diagnostics <- diagnostics[
    order(diagnostics$donor_id, as.numeric(diagnostics$distance_bin_start)),
    , drop = FALSE
  ]
}
write_pipeline_table(donor_calibration, donor_output_file)
write_pipeline_table(diagnostics, diagnostics_output_file)

# donor_fit_is_usable() also re-checks that a fitted value is a distance.
fitted <- donors$threshold[donor_fit_is_usable(donors$threshold,
                                               donors$fit_status == "fitted")]
choice <- choose_scope_clone_threshold(fitted, configured_threshold)
if (identical(choice$source, "automatic_fallback")) {
  warning(sprintf(
    "scope %s used the automatic fallback clone threshold %.4f; no donor fit was usable",
    scope_id, choice$threshold
  ))
}

# n_donors counts every donor row; n_donors_fitted counts the usable fits, as
# choose_scope_clone_threshold() reports them.
out <- data.frame(
  scope_id = scope_id,
  clone_threshold = choice$threshold,
  threshold_source = choice$source,
  n_donors = nrow(donors),
  n_donors_fitted = choice$n_donors_fitted,
  stringsAsFactors = FALSE
)
write_pipeline_table(out, output_file)

log_info(
  "Scope %s threshold %.4f (%s; %d of %d donors fitted)",
  scope_id, choice$threshold, choice$source, choice$n_donors_fitted, nrow(donors)
)
log_info("Retained %d donor calibration row(s) and %d diagnostic bin row(s)",
         nrow(donor_calibration), nrow(diagnostics))

finish_stage_log(stage_log)
