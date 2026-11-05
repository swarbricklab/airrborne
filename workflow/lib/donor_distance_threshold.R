# Estimate an IGH clone distance cutoff from one donor's libraries in one scope.
# A scope selects the libraries analysed together; donors are never compared here.
# Preparation in igh_grouping.R first chooses eligible sequence representatives.
# This module measures how different each representative's junction is from its
# closest comparable neighbour, then fits a cutoff between near and far distances.
# It keeps diagnostic distances even when there is too little evidence to fit.
#
# SHazaM 1.3.2 distToNearest measures distances (Gupta et al. 2015,
# doi:10.1093/bioinformatics/btv359). findThreshold fits a gamma-normal mixture
# (Nouri and Kleinstein 2018, doi:10.3389/fimmu.2018.01687). Change the distance
# definition here and in assign_scope_donor_clones.R together. Fit settings come
# from params.lineage; scope_clone_threshold.R chooses the cutoff assignment uses.

DISTANCE_HISTOGRAM_BREAKS <- seq(0, 1, by = 0.01)

# The released name of each fit parameter and the SHazaM 1.3.2 GmmThreshold slot
# holding it. A gamma-norm fit stores the gamma near component's weight, shape and
# scale in a1, b1 and c1, and the normal far component's weight, mean and standard
# deviation in a2, b2 and c2.
DONOR_FIT_PARAMETERS <- c(
  near_weight = "a1", near_shape = "b1", near_scale = "c1",
  far_weight = "a2", far_mean = "b2", far_sd = "c2",
  absolute_log_likelihood = "loglk", fit_sensitivity = "sensitivity"
)

# The one-row donor fit as it is written and read back.
DONOR_THRESHOLD_COLUMNS <- c(
  "scope_id", "donor_id", "n_libraries", "n_members", "n_distances",
  "n_distances_fitted", "fit_status", "threshold", names(DONOR_FIT_PARAMETERS)
)

# Every donor fit ends in one of these states, and only a fitted donor reports a
# threshold. too_few_distances and constant_distances skip the fit; fit_failed
# covers a findThreshold error and a fit without a finite threshold; out_of_range
# is a finite threshold outside [0, 1].
DONOR_FIT_STATUSES <- c("fitted", "too_few_distances", "constant_distances",
                        "fit_failed", "out_of_range")

# The per-bin diagnostic table as it is written and read back.
DONOR_DIAGNOSTIC_COLUMNS <- c(
  "scope_id", "donor_id", "distance_bin_start", "distance_bin_end",
  "distance_bin_midpoint", "n_distances_in_bin", "observed_density",
  "near_component_density", "far_component_density"
)

# ------------------------------------------------------------ determinism ---

#' Replace SHazaM's random seed resets with PIPELINE_SEED for repeatable fits.
#'
#' gmmFit() and rocSpace() each call set.seed(NULL), drawing a new random seed
#' after findThreshold() has sampled the distances. Setting the pipeline seed
#' before findThreshold() alone therefore cannot reproduce mixture initialisation.
#' This function changes only those resets in the process-local shazam namespace;
#' the package's fitting calculations remain unchanged. Returns nothing.
#'
#' The function-body traversal requires exactly one reset in each function and
#' stops otherwise. A SHazaM upgrade needs review of these internals and the
#' repeatability tests in test/components/igh_clones/unit/test-donor_threshold.R.
#' seed remains an argument so those checks can exercise the replacement itself.
override_shazam_gmm_seed <- function(seed) {
  for (function_name in c("gmmFit", "rocSpace")) {
    fit_function <- get(
      function_name, envir = asNamespace("shazam"), inherits = FALSE
    )
    replacements <- 0L
    replace_seed_reset <- function(expr) {
      if (!is.call(expr)) return(expr)
      if (identical(expr[[1L]], as.name("set.seed")) && length(expr) == 2L &&
          identical(expr[[2L]], quote(NULL))) {
        replacements <<- replacements + 1L
        return(call("set.seed", as.integer(seed)))
      }
      for (index in seq_along(expr)) {
        # Single-bracket replacement preserves literal NULL arguments and
        # assignments; double brackets would delete those AST elements.
        expr[index] <- list(replace_seed_reset(expr[[index]]))
      }
      expr
    }
    body(fit_function) <- replace_seed_reset(body(fit_function))
    if (replacements != 1L) {
      stop("Expected one internal SHazaM seed reset in ", function_name,
           "; found ", replacements, " in shazam ",
           as.character(utils::packageVersion("shazam")))
    }
    assignInNamespace(function_name, fit_function, ns = "shazam")
  }
}

# ---------------------------------------------------------------- decisions ---

#' Why the distance evidence cannot support a mixture fit, if it cannot.
#'
#' `distances` are the finite nearest-neighbour distances and `minimum_n` the
#' fewest a fit is attempted on. Returns "too_few_distances" below that count,
#' "constant_distances" when the distances take one value, and NA when a fit can
#' be attempted. The count is checked first. Passing both checks does not
#' guarantee a usable fitted threshold.
distance_evidence_status <- function(distances, minimum_n) {
  if (length(distances) < minimum_n) return("too_few_distances")
  if (length(unique(distances)) < 2L) return("constant_distances")
  NA_character_
}

#' Is a fitted cutoff usable as a clone distance?
#'
#' `threshold` is the finite value taken from the fit. Returns one logical. A
#' normalised Hamming distance outside [0, 1] is not a distance, so a fit that
#' produced one is reported as out_of_range rather than passed on to clustering.
is_threshold_in_range <- function(threshold) {
  is.finite(threshold) && threshold >= 0 && threshold <= 1
}

#' Classify the cutoff an attempted fit returned.
#'
#' `threshold` is the fitted value, or NA when findThreshold raised an error or
#' returned no finite threshold. Returns "fit_failed" for NA, "out_of_range" for
#' a value outside [0, 1] and "fitted" otherwise.
fit_threshold_status <- function(threshold) {
  if (!is.finite(threshold)) return("fit_failed")
  if (!is_threshold_in_range(threshold)) return("out_of_range")
  "fitted"
}

# ------------------------------------------------------------------- the fit ---

#' Measure nearest junction distances and fit a cutoff for one donor in one scope.
#'
#' db contains one representative per distinct retained full sequence within each
#' library, selected by prepared_igh_representatives(). Preparation transfers read
#' support; this fit uses sequences, without weighting distances by that support.
#' minimum_n is params.lineage.min_fit_distances. fit_edge, fit_specificity and
#' fit_subsample pass the corresponding params.lineage settings to findThreshold;
#' a NULL subsample uses every finite distance. PIPELINE_SEED controls sampling
#' and mixture initialisation, including the package's internal seed resets.
#'
#' Hamming distance counts differing nucleotide positions; normalize = "len"
#' divides by junction length. Comparisons stay within equal junction lengths and
#' exact resolved V/J reference-sequence groups. igh_inference_calls() encodes
#' those groups as temporary gene calls for SHazaM and SCOPer, preserving tied
#' gene names only when they resolve to the same complete reference sequence.
#' Both tools receive first = FALSE; the temporary single calls keep these groups
#' separate. symmetry = "avg" averages directional distances. A record without a
#' comparable neighbour has infinite distance and contributes no finite evidence
#' to the fit.
#'
#' Returns finite distances, a usable fit or NULL, threshold, fit_status (one of
#' DONOR_FIT_STATUSES), a caught findThreshold error or NULL, and
#' n_distances_fitted. Too few distances or a constant distribution skips fitting
#' and leaves n_distances_fitted missing; failed or out-of-range fits have a
#' missing threshold. The scope later combines usable donor fits or applies its
#' configured cutoff or fallback.
#'
#' Bulk RNA-seq identifies no individual cell. cellIdColumn stays NULL:
#' assembly_backbone_id names an assembled contig and cannot substitute for it.
fit_donor_distance_threshold <- function(db, minimum_n, fit_edge,
                                         fit_specificity, fit_subsample) {
  db <- igh_inference_calls(db)
  distance_db <- shazam::distToNearest(
    db,
    sequenceColumn = "junction",
    vCallColumn = "v_inference_call",
    jCallColumn = "j_inference_call",
    model = "ham",
    normalize = "len",
    symmetry = "avg",
    first = FALSE,
    VJthenLen = TRUE,
    nproc = 1,
    cellIdColumn = NULL,
    locusColumn = "locus",
    locusValues = "IGH"
  )
  # findThreshold() subsamples by position, and row order follows TRUST4 record
  # numbering, which changes between runs. Sorting makes the fit depend only on
  # the distances.
  distances <- sort(distance_db$dist_nearest[is.finite(distance_db$dist_nearest)])

  skipped <- distance_evidence_status(distances, minimum_n)
  if (!is.na(skipped)) {
    return(list(distances = distances, fit = NULL, threshold = NA_real_,
                fit_status = skipped, error = NULL, n_distances_fitted = NA_integer_))
  }

  set.seed(PIPELINE_SEED)
  override_shazam_gmm_seed(PIPELINE_SEED)
  failure <- NULL
  fit <- tryCatch(
    shazam::findThreshold(
      distances,
      method = "gmm",
      model = "gamma-norm",
      edge = fit_edge,
      cutoff = "user",
      spc = fit_specificity,
      subsample = fit_subsample,
      progress = FALSE
    ),
    error = function(error) {
      failure <<- conditionMessage(error)
      NULL
    }
  )
  threshold <- if (!is.null(fit) && isS4(fit) && "threshold" %in% slotNames(fit) &&
                   length(fit@threshold) == 1L && is.finite(fit@threshold)) {
    as.numeric(fit@threshold)
  } else {
    NA_real_
  }
  # The status is recorded before an unusable value is cleared, so an
  # out-of-range cutoff stays distinguishable from a failed fit.
  fit_status <- fit_threshold_status(threshold)
  fitted <- identical(fit_status, "fitted")
  if (!fitted) threshold <- NA_real_

  # findThreshold() samples min(length(distances), subsample) distances, and a
  # NULL subsample uses them all. Reading the count from these values keeps it
  # available for an attempted fit that returned nothing.
  n_distances_fitted <- if (is.null(fit_subsample)) {
    length(distances)
  } else {
    min(length(distances), fit_subsample)
  }

  list(distances = distances, fit = if (fitted) fit else NULL,
       threshold = threshold, fit_status = fit_status, error = failure,
       n_distances_fitted = as.integer(n_distances_fitted))
}

# ---------------------------------------------------------------- the tables ---

#' Build the one-row donor fit.
#'
#' `scope_id` and `donor_id` name the pool. `n_libraries` and `n_members` count
#' the donor's libraries and the sequence representatives that entered the fit.
#' The remaining arguments default to the values a donor with no records reports.
#' `fit` is the usable GmmThreshold or NULL, and supplies DONOR_FIT_PARAMETERS,
#' which stay missing without one. Returns a one-row frame carrying
#' DONOR_THRESHOLD_COLUMNS; stops on a fit_status outside DONOR_FIT_STATUSES.
build_donor_threshold_row <- function(scope_id, donor_id, n_libraries = 0L,
                                      n_members = 0L, n_distances = 0L,
                                      threshold = NA_real_,
                                      fit_status = "too_few_distances",
                                      n_distances_fitted = NA_integer_,
                                      fit = NULL) {
  if (length(fit_status) != 1L || !fit_status %in% DONOR_FIT_STATUSES) {
    stop("fit_status must be one of ", paste(DONOR_FIT_STATUSES, collapse = ", "))
  }
  row <- data.frame(
    scope_id = scope_id,
    donor_id = donor_id,
    n_libraries = as.integer(n_libraries),
    n_members = as.integer(n_members),
    n_distances = as.integer(n_distances),
    n_distances_fitted = as.integer(n_distances_fitted),
    fit_status = fit_status,
    threshold = as.numeric(threshold),
    stringsAsFactors = FALSE
  )
  for (column in names(DONOR_FIT_PARAMETERS)) {
    row[[column]] <- if (is.null(fit)) {
      NA_real_
    } else {
      as.numeric(slot(fit, DONOR_FIT_PARAMETERS[[column]]))
    }
  }
  row
}

#' Build the donor's distance histogram and, where one was fitted, its mixture.
#'
#' `scope_id` and `donor_id` name the pool, `distances` are the finite
#' nearest-neighbour distances, and `fit` is the GmmThreshold or NULL. Returns
#' one row per bin of DISTANCE_HISTOGRAM_BREAKS, carrying
#' DONOR_DIAGNOSTIC_COLUMNS.
#'
#' The observed histogram is written whether or not a model was fitted, so a
#' donor too sparse to fit still leaves a distance distribution someone can look
#' at. A donor with no finite distances gets every bin at count zero and observed
#' density zero; those zeros record an empty histogram and not a measured density.
#' The two component densities are the fitted near (same-clone, gamma) and far
#' (different-clone, normal) parts evaluated at each bin midpoint, so the fit
#' can be read against the data it was fitted to; they are NA when there is no
#' fit, and when a component evaluates to a non-finite density.
build_distance_histogram_diagnostics <- function(scope_id, donor_id,
                                                 distances = numeric(0),
                                                 fit = NULL) {
  breaks <- DISTANCE_HISTOGRAM_BREAKS
  counts <- if (length(distances)) {
    hist(distances, breaks = breaks, plot = FALSE,
         include.lowest = TRUE, right = FALSE)$counts
  } else {
    integer(length(breaks) - 1L)
  }
  width <- diff(breaks)
  observed_density <- if (sum(counts) > 0L) counts / (sum(counts) * width) else rep(0, length(counts))
  midpoint <- head(breaks, -1L) + width / 2
  near <- far <- rep(NA_real_, length(midpoint))

  if (!is.null(fit) && isS4(fit) && is(fit, "GmmThreshold") &&
      identical(as.character(fit@model), "gamma-norm")) {
    near <- fit@a1 * dgamma(midpoint, shape = fit@b1, scale = fit@c1)
    far <- fit@a2 * dnorm(midpoint, mean = fit@b2, sd = fit@c2)
    near[!is.finite(near)] <- NA_real_
    far[!is.finite(far)] <- NA_real_
  }

  data.frame(
    scope_id = scope_id,
    donor_id = donor_id,
    distance_bin_start = head(breaks, -1L),
    distance_bin_end = tail(breaks, -1L),
    distance_bin_midpoint = midpoint,
    n_distances_in_bin = as.integer(counts),
    observed_density = observed_density,
    near_component_density = near,
    far_component_density = far,
    stringsAsFactors = FALSE
  )
}

#' Read one donor's fit back.
#'
#' `path` is the table build_donor_threshold_row() wrote. Returns its single row.
read_donor_threshold_row <- function(path) {
  row <- read_pipeline_table(
    path,
    required_columns = DONOR_THRESHOLD_COLUMNS,
    context = paste("Donor threshold table", path)
  )
  if (nrow(row) != 1L) stop("Invalid donor threshold table: ", path)
  row
}

#' Read one donor's diagnostic table back.
#'
#' `path` is the table build_distance_histogram_diagnostics() wrote. Returns its
#' rows, projected to DONOR_DIAGNOSTIC_COLUMNS. The bin count is asserted against
#' the same breaks the writer used, because a short histogram means the writer and
#' the reader no longer agree on what one bin is, and the pooled scope table would
#' otherwise silently mix two grids.
read_donor_distance_diagnostics <- function(path) {
  rows <- read_pipeline_table(
    path, required_columns = DONOR_DIAGNOSTIC_COLUMNS,
    context = paste("Donor distance diagnostics table", path)
  )
  if (nrow(rows) != length(DISTANCE_HISTOGRAM_BREAKS) - 1L) {
    stop("Invalid donor distance diagnostics table: ", path)
  }
  rows[, DONOR_DIAGNOSTIC_COLUMNS, drop = FALSE]
}
