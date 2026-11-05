# Choose one IGH clone distance cutoff for every donor analysed in a scope.
# Donors are clustered separately but share this cutoff. Donor fits run even with
# a numeric params.lineage.clonal_threshold, retaining evidence for review.
# This module chooses the value; assign_scope_donor_clones.R applies it with SCOPer.

# Automatic mode uses this value only when no donor fit is usable.
AUTOMATIC_CLONE_THRESHOLD_FALLBACK <- 0.1

#' Check for finite junction distance cutoffs from zero to one (vectorised).
#'
#' x is numeric. Missing or non-finite values return FALSE. The configured cutoff,
#' donor fits and saved scope cutoff all use this check: the distance is the
#' number of differing nucleotides divided by junction length.
is_clone_distance <- function(x) {
  is.finite(x) & x >= 0 & x <= 1
}

#' Read the configured threshold argument.
#'
#' `value` is the `clonal_threshold` setting: the literal "auto", or a number from
#' 0 to 1. Returns NA_real_ for "auto", so that "no number was configured" is
#' carried as a value rather than as a second flag, and the number otherwise.
#' Stops on anything else, because a mistyped cutoff would otherwise silently
#' become a fallback.
parse_configured_threshold <- function(value) {
  if (identical(value, "auto")) return(NA_real_)
  numeric_threshold <- suppressWarnings(as.numeric(value))
  if (!is_clone_distance(numeric_threshold)) {
    stop("clonal_threshold must be 'auto' or a number between 0 and 1")
  }
  numeric_threshold
}

#' Which donor fits may inform the scope cutoff (vectorised).
#'
#' `threshold` is the donor fit table's cutoff and `fitted` says whether its
#' `fit_status` is fitted. Returns one logical per donor. A donor counts only when
#' its own stage said the fit succeeded and the value it reported is a distance:
#' an unreadable flag is not treated as a pass.
#'
#' The range is re-checked here rather than trusted from the status because the
#' two come back from a file the donor stage wrote, and a truncated or hand-edited
#' calibration table would otherwise move the scope median unnoticed.
donor_fit_is_usable <- function(threshold, fitted) {
  !is.na(fitted) & fitted & is_clone_distance(threshold)
}

#' Check that the donor fits and their diagnostics describe this scope.
#'
#' `donors` is the donor fit table, `diagnostics` the pooled per-donor histograms,
#' and `scope_id` the scope being frozen. Returns nothing; stops when the inputs
#' cannot be the inputs of this scope.
#'
#' The two tables are assembled from separate lists of files, so this is where a
#' wiring mistake surfaces: a file from another scope, one donor's fit counted
#' twice, or a donor present in one table and absent from the other. Any of those
#' would move the median without being visible in the released row.
validate_donor_fit_consistency <- function(donors, diagnostics, scope_id) {
  if (nrow(donors) && any(as.character(donors$scope_id) != scope_id)) {
    stop("donor threshold table contains a different scope_id")
  }
  if (nrow(diagnostics) && any(as.character(diagnostics$scope_id) != scope_id)) {
    stop("distance diagnostics contain a different scope_id")
  }
  if (anyDuplicated(as.character(donors$donor_id))) {
    stop("donor threshold table contains a duplicate donor_id")
  }
  if (!setequal(as.character(donors$donor_id),
                unique(as.character(diagnostics$donor_id)))) {
    stop("donor threshold and distance diagnostic tables cover different donors")
  }
  invisible(NULL)
}

#' Choose the cutoff used for assignment, separately from fitting donor evidence.
#'
#' donor_thresholds contains only fits accepted by donor_fit_is_usable(). configured
#' is the parsed params.lineage.clonal_threshold: a number wins; NA means auto.
#' Auto uses the median usable donor cutoff, giving each donor one vote. With no
#' usable donor fit it uses AUTOMATIC_CLONE_THRESHOLD_FALLBACK, unless the caller
#' supplies fallback explicitly. Change that constant to change the pipeline's
#' fallback, and update the configuration documentation and threshold tests.
#'
#' Returns threshold, source (configured, automatic_fit or automatic_fallback) and
#' n_donors_fitted. The count describes available fits even with a numeric cutoff.
#' The caller saves this choice with the donor evidence and warns on fallback.
choose_scope_clone_threshold <- function(donor_thresholds, configured,
                                           fallback = AUTOMATIC_CLONE_THRESHOLD_FALLBACK) {
  if (!is.na(configured)) {
    return(list(threshold = configured, source = "configured",
                n_donors_fitted = length(donor_thresholds)))
  }
  if (length(donor_thresholds) > 0L) {
    return(list(threshold = as.numeric(median(donor_thresholds)),
                source = "automatic_fit",
                n_donors_fitted = length(donor_thresholds)))
  }
  list(threshold = fallback, source = "automatic_fallback",
       n_donors_fitted = length(donor_thresholds))
}

#' Read a scope's frozen lineage distance back.
#'
#' `path` is the one-row scope threshold table and `scope_id` the scope the caller
#' believes it is working on. Returns the cutoff as one number.
#'
#' This is the consuming end of the table choose_scope_clone_threshold() feeds.
#' The scope is checked rather than assumed because clone assignment runs once
#' per donor from a path built out of wildcards: handed another scope's table it
#' would cluster happily at the wrong distance.
read_scope_clone_threshold <- function(path, scope_id) {
  threshold_row <- read_pipeline_table(
    path, required_columns = c("scope_id", "clone_threshold"),
    context = "Scope threshold table"
  )
  if (nrow(threshold_row) != 1L) {
    stop("Scope threshold table must contain exactly one row")
  }
  if (as.character(threshold_row$scope_id[1]) != scope_id) {
    stop("Scope threshold table scope_id does not match expected scope_id ", scope_id)
  }
  threshold <- suppressWarnings(as.numeric(threshold_row$clone_threshold[1]))
  if (!is_clone_distance(threshold)) {
    stop("Scope threshold table clone_threshold must be between 0 and 1")
  }
  threshold
}
