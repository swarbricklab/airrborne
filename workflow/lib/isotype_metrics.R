# Calculate how recovered IGH read support is distributed across isotypes.
# build_isotype_metrics() defines the totals and fractions for one library.
# Each fraction uses only support assigned at that level, class or subclass.
# A call can count towards a class while its subclass remains unknown.
# isotype_calls.R checks the call's evidence; isotype_resolution.R names it.
# The IgG:IgA log ratio needs positive support for both classes. No extra count
# is added when either is zero. These helpers use base R.

ISOTYPE_COMPOSITION_LEVELS <- list(
  class = c("IgM", "IgD", "IgG", "IgA", "IgE"),
  subclass = CANONICAL_IGH_ISOTYPES
)

# The released isotype_usage.tsv columns after the library identifier.
ISOTYPE_METRIC_COLUMNS <- c(
  "isotype_level", "isotype", "consensus_count", "callable_consensus_count",
  "support_fraction"
)

ISOTYPE_SUMMARY_COLUMNS <- c(
  "common_qc_consensus_count", "class_callable_consensus_count",
  "class_uncallable_consensus_count", "subclass_callable_consensus_count",
  "class_only_consensus_count", "igm_consensus_count", "igg_consensus_count",
  "iga_consensus_count", "log2_igg_iga", "igg_iga_status"
)

#' Summarise one library's class/subclass composition and finite IgG:IgA ratio.
#'
#' records contains qualified_c_call and consensus_count for all common-QC IGH rows,
#' with the call already masked where its alignment evidence failed. Return a list
#' of three elements: 14 composition rows, one summary row, and the per-record
#' resolutions at each level, which the caller uses for the callable audit rows.
#' A composition row's `callable_consensus_count` is the support resolving at
#' its `isotype_level`, the denominator of its `support_fraction`. Empty input
#' has zero support counts and missing fractions; unresolved calls retain their
#' support in the summary's common-QC total, outside every callable count.
build_isotype_metrics <- function(records) {
  require_columns(records, c("qualified_c_call", "consensus_count"), "Qualified IGH records")
  counts <- as_read_support(records$consensus_count, "consensus_count")
  total <- sum(counts)
  resolved <- lapply(names(ISOTYPE_COMPOSITION_LEVELS), function(resolution) {
    resolve_isotype(records$qualified_c_call, resolution)
  })
  names(resolved) <- names(ISOTYPE_COMPOSITION_LEVELS)

  tables <- lapply(names(resolved), function(resolution) {
    calls <- resolved[[resolution]]
    callable <- sum(counts[!is.na(calls)])
    levels <- ISOTYPE_COMPOSITION_LEVELS[[resolution]]
    support <- vapply(levels, function(isotype) {
      sum(counts[!is.na(calls) & calls == isotype])
    }, numeric(1), USE.NAMES = FALSE)
    data.frame(
      isotype_level = resolution, isotype = levels, consensus_count = support,
      callable_consensus_count = callable,
      support_fraction = if (callable > 0) support / callable else NA_real_,
      stringsAsFactors = FALSE
    )
  })
  composition <- do.call(rbind, tables)
  class_rows <- composition[composition$isotype_level == "class", , drop = FALSE]
  class_support <- stats::setNames(class_rows$consensus_count, class_rows$isotype)
  class_callable <- sum(class_support)
  subclass_callable <- sum(counts[!is.na(resolved$subclass)])
  igg <- unname(class_support["IgG"])
  iga <- unname(class_support["IgA"])
  status <- if (class_callable == 0) {
    "no_class_callable_evidence"
  } else if (igg == 0 && iga == 0) {
    "both_zero"
  } else if (igg == 0) {
    "zero_igg"
  } else if (iga == 0) {
    "zero_iga"
  } else {
    "available"
  }

  summary <- data.frame(
    common_qc_consensus_count = total,
    class_callable_consensus_count = class_callable,
    class_uncallable_consensus_count = total - class_callable,
    subclass_callable_consensus_count = subclass_callable,
    class_only_consensus_count = sum(counts[!is.na(resolved$class) &
                                            is.na(resolved$subclass)]),
    igm_consensus_count = unname(class_support["IgM"]),
    igg_consensus_count = igg,
    iga_consensus_count = iga,
    # Subtract logarithms so even very unequal finite inputs avoid division
    # overflow or underflow before the logarithm is taken.
    log2_igg_iga = if (status == "available") log2(igg) - log2(iga) else NA_real_,
    igg_iga_status = status,
    stringsAsFactors = FALSE
  )
  list(composition = composition[, ISOTYPE_METRIC_COLUMNS, drop = FALSE],
       summary = summary[, ISOTYPE_SUMMARY_COLUMNS, drop = FALSE], resolved = resolved)
}
