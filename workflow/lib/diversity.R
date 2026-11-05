# Calculate observed diversity, dominance and analytical rarefaction from released
# read-support vectors. IGH units are donor/scope clones; other loci use exact
# clonotypes. Interpolation samples support without replacement and estimates no
# unseen groups; the q = 2 asymptote in q2_asymptote() is the one model-based
# estimate. These helpers use base R; method details accompany each estimator.

DIVERSITY_COLS <- c(
  "library_id", "locus", "group_unit", "consensus_count", "richness",
  "shannon_entropy", "inverse_simpson",
  "clonality", "top_1_support_fraction",
  "top_5_support_fraction", "top_10_support_fraction"
)

# The per-library clone-size vector the metrics above are computed from, released
# so a reader can re-derive them or compute an index this pipeline does not.
CLONE_ABUNDANCE_COLS <- c(
  "library_id", "locus", "group_unit", "group_key", "n_members",
  "consensus_count", "support_fraction", "group_rank"
)

# Concentration uses the complete within-library support proportions. Pairwise
# change and observed inverse Simpson share this definition.
simpson_concentration <- function(proportions) {
  sum(proportions^2)
}

# A minimum controls estimate eligibility, while an empty repertoire has no
# fraction denominator. Pair comparisons use the same endpoint decision.
repertoire_support_status <- function(total, minimum_consensus_count = 0) {
  if (length(minimum_consensus_count) != 1L ||
      !is.finite(minimum_consensus_count) || minimum_consensus_count < 0 ||
      minimum_consensus_count != floor(minimum_consensus_count)) {
    stop("minimum_consensus_count must be a non-negative integer")
  }
  if (total == 0) "empty_repertoire" else if (
    total < minimum_consensus_count) "below_minimum_support" else "available"
}

#' Sum one library's read support within each donor-pooled clone.
#'
#' `membership` is a scope's whole clone-membership table, `library_id` names the
#' library to take out of it, and `clone_calls` holds the scope-wide V and J
#' genes `summarise_clone_vj_call_sets()` voted, with the status of each
#' vote. Clones
#' are pooled across a donor's samples, so this is the step that turns a shared
#' clone back into one library's observation of it: member read support is
#' summed within the clone, and the clone then counts once, whatever its
#' member sequences did.
#'
#' Rows never assigned to a clone are outside the observed repertoire and are
#' dropped here rather than inside a metric. as_read_support() checks support
#' before aggregation. compute_diversity_row() also rejects zero-support groups,
#' which cannot represent an observed clone.
#'
#' Returns one row per clone the library contributed to (`clone_id`, the summed
#' `consensus_count`, how many members `n_members` the library contributed,
#' the code-point smallest junction `min_junction` among those members, and the
#' voted V and J gene with each vote's status) and no rows when it contributed
#' none. `compute_diversity_row()`, `build_group_abundance_rows()` and
#' `count_gene_usage()` all read this same aggregation, which is also what
#' puts every one of the library's clones into exactly one usage or exclusion
#' count.
aggregate_clone_support_per_library <- function(membership, library_id, clone_calls) {
  library_members <- membership[as.character(membership$library_id) == library_id, ,
                                drop = FALSE]
  if (nrow(library_members)) {
    library_members$consensus_count <- as_read_support(
      library_members$consensus_count, "consensus_count")
    library_members <- library_members[
      !is.na(library_members$clone_id) & library_members$clone_id != "", ,
      drop = FALSE]
  }
  if (!nrow(library_members)) {
    return(data.frame(
      clone_id = character(0), consensus_count = numeric(0),
      n_members = integer(0), min_junction = character(0),
      v_gene = character(0), j_gene = character(0),
      v_vote_status = character(0), j_vote_status = character(0),
      stringsAsFactors = FALSE
    ))
  }

  # Clones come out in code-point order, which set_pipeline_collation() fixes
  # for the whole pipeline, so the order depends on the identifiers and not on
  # which image's locale ran the stage.
  clone_id <- as.character(library_members$clone_id)
  members <- split(
    seq_along(clone_id),
    factor(clone_id, levels = sort(unique(clone_id)))
  )
  lineages <- data.frame(
    clone_id = as.character(names(members)),
    consensus_count = vapply(
      members,
      function(rows) sum(library_members$consensus_count[rows]),
      numeric(1), USE.NAMES = FALSE),
    # Distinct retained full sequences of the clone in this library, the unit
    # n_members has for clonotypes too.
    n_members = vapply(members, length, integer(1), USE.NAMES = FALSE),
    # The content tie-break build_group_abundance_rows() orders equally
    # supported clones by, since SCOPer's clone labels are not content.
    min_junction = vapply(
      members,
      function(rows) min(as.character(library_members$junction[rows])),
      character(1), USE.NAMES = FALSE),
    stringsAsFactors = FALSE
  )
  matched <- match(lineages$clone_id, clone_calls$clone_id)
  lineages$v_gene <- clone_calls$v_gene[matched]
  lineages$j_gene <- clone_calls$j_gene[matched]
  lineages$v_vote_status <- clone_calls$v_vote_status[matched]
  lineages$j_vote_status <- clone_calls$j_vote_status[matched]
  lineages
}

#' Return the fixed diversity row for a library with no eligible clones.
make_empty_diversity_row <- function(sample_name, locus, level) {
  data.frame(
    library_id = sample_name,
    locus = locus,
    group_unit = level,
    consensus_count = 0,
    richness = 0L,
    shannon_entropy = NA_real_,
    inverse_simpson = NA_real_,
    clonality = NA_real_,
    top_1_support_fraction = NA_real_,
    top_5_support_fraction = NA_real_,
    top_10_support_fraction = NA_real_,
    stringsAsFactors = FALSE
  )
}

#' Summarise one already-aggregated repertoire.
#'
#' `df` is one row per observed clone with its summed `consensus_count`, and every row
#' must carry positive support. Which clones count as observed is decided by the caller
#' that aggregates them (aggregate_clone_support_per_library() for IGH clones,
#' build_rearrangement_group_table() for the non-IGH groups), so a reader can see the
#' inclusion rule beside the records instead of having it applied invisibly inside the
#' metric.
#'
#' Let p be each clone's share of the library's total read support at this locus.
#' The published metrics are then:
#'   richness         the number of clones observed.
#'   shannon_entropy  -sum(p * log(p)), in nats, because the logarithm is
#'                    natural. Shannon 1948, doi:10.1002/j.1538-7305.1948.tb01338.x
#'   inverse_simpson  1 / sum(p^2), the Hill number of order 2, i.e. the number
#'                    of equally supported clones that would give this Simpson
#'                    concentration. Hill 1973, doi:10.2307/1934352
#'   clonality
#'                    1 - Shannon/log(richness), i.e. one minus Pielou's evenness
#'                    J. Pielou 1966, doi:10.1016/0022-5193(66)90013-0.
#'                    Evenness needs two clones to mean anything, so this is NA
#'                    for a repertoire holding one.
#'   top_N_support_fraction
#'                    the summed share of up to N largest observed groups.
#'                    It is one when a nonempty repertoire has fewer than N.
#'
#' Returns one row in DIVERSITY_COLS order.
compute_diversity_row <- function(df, sample_name, locus, level) {
  if (nrow(df) == 0L) return(make_empty_diversity_row(sample_name, locus, level))

  # Zero-support groups inflate observed richness. Reject them instead of
  # silently changing the population supplied by the caller.
  counts <- as_read_support(df$consensus_count, "consensus_count")
  if (any(counts == 0)) {
    stop(sprintf("[%s/%s] every clone must carry finite positive read support",
                 locus, level))
  }

  consensus_count <- sum(counts)
  richness <- length(counts)
  p <- counts / consensus_count
  shannon <- -sum(p * log(p))
  ranked_p <- sort(p, decreasing = TRUE)

  top_support_fraction <- function(n) {
    sum(ranked_p[seq_len(min(n, richness))])
  }

  # Clonality is defined on [0, 1], but Shannon is summed independently of
  # log(richness) here, so a perfectly even repertoire can land a fraction of an
  # ulp either side of an endpoint and be written out in full by the wide scipen
  # setting. Clamping to the two theoretical bounds keeps the released column
  # inside the range the data dictionary declares, which is what a reader passing
  # it to logit(), sqrt() or a beta model relies on.
  clonality <- function() {
    if (richness < 2L) return(NA_real_)
    min(1, max(0, 1 - shannon / log(richness)))
  }

  data.frame(
    library_id = sample_name,
    locus = locus,
    group_unit = level,
    consensus_count = consensus_count,
    richness = richness,
    shannon_entropy = shannon,
    inverse_simpson = 1 / simpson_concentration(p),
    clonality = clonality(),
    top_1_support_fraction = top_support_fraction(1L),
    top_5_support_fraction = top_support_fraction(5L),
    top_10_support_fraction = top_support_fraction(10L),
    stringsAsFactors = FALSE
  )[, DIVERSITY_COLS, drop = FALSE]
}

#' The clone-size vector one diversity row was computed from.
#'
#' `df` is the same already-aggregated repertoire compute_diversity_row() is
#' handed, `group_unit` names what one row of it is (`clone` for IGH or
#' `clonotype` elsewhere) and `key_column` says which of its columns identifies
#' the clone. Returns one row per clone in CLONE_ABUNDANCE_COLS order.
#'
#' `group_rank` is rank(-consensus_count, ties.method = "min"): 1 for the
#' best-supported clone, and equally supported clones share the smaller rank,
#' so ranks run 1, 2, 2, 4. A top-k selection by rank can therefore hold more
#' than k clones. The top_N_support_fraction values do not depend on the tie,
#' because tied clones carry equal support. Rows are ordered by group_rank,
#' then by `order_columns` as text, which gives tied clones an order that
#' follows their content. The non-IGH caller orders by clonotype_key; the IGH
#' caller orders by each clone's smallest member junction, then clone_id.
#'
#' The support vector is released so these measurements can be reconstructed.
#' `df$n_members` counts the distinct retained member sequences of each clone or
#' clonotype. It is descriptive context and is not a cell count or an
#' alternative sampling unit for these support-weighted measurements.
build_group_abundance_rows <- function(df, sample_name, locus, group_unit,
                                       key_column, order_columns = key_column) {
  if (nrow(df) == 0L) return(make_empty_group_abundance_table())

  counts <- as_read_support(df$consensus_count, "consensus_count")
  group_rank <- rank(-counts, ties.method = "min")
  rows <- do.call(order, c(list(group_rank),
                           lapply(order_columns, function(column) {
                             as.character(df[[column]])
                           })))
  out <- data.frame(
    library_id = sample_name,
    locus = locus,
    group_unit = group_unit,
    group_key = as.character(df[[key_column]]),
    group_rank = group_rank,
    n_members = as.integer(df$n_members),
    consensus_count = counts,
    support_fraction = counts / sum(counts),
    stringsAsFactors = FALSE
  )[rows, CLONE_ABUNDANCE_COLS, drop = FALSE]
  rownames(out) <- NULL
  out
}

#' The clone-abundance table of a library that observed no clone at a locus.
#'
#' Typed rather than make_empty_frame()'s all-character, because the counts and
#' the support fraction are numbers everywhere else this table appears.
make_empty_group_abundance_table <- function() {
  data.frame(
    library_id = character(),
    locus = character(),
    group_unit = character(),
    group_key = character(),
    n_members = integer(),
    consensus_count = numeric(),
    support_fraction = numeric(),
    group_rank = integer(),
    stringsAsFactors = FALSE
  )
}

#
# Hill diversity interpolated within a library's retained group support.
#
# At depth m, q = 0 is expected richness, q = 1 is exp(expected Shannon
# entropy), and q = 2 is reciprocal expected concentration. These are the
# size-based interpolation conventions of Chao et al. (2014), not averages of
# transformed diversity values over draws. Calculations use the observed
# support vector without replacement; they do not predict cell diversity or
# recovery from a shallower FASTQ experiment.
#
# Method: Chao et al. 2014 (doi:10.1890/13-0133.1).
# Reference implementation: iNEXT::RTD in src/subfuns.cpp, iNEXT 3.0.2.
# Unlike iNEXT, no curve point or target lies above observed support: that
# would be extrapolation, which estimates unseen groups, and this module
# estimates none. The depth grid holds every integer depth when total support
# is at most RAREFACTION_GRID_POINTS units, and RAREFACTION_GRID_POINTS
# log-spaced depths otherwise.
# stats::phyper() and stats::dhyper() evaluate hypergeometric probabilities
# directly, avoiding ratios of binomial coefficients that can overflow.

# Forty log-spaced depths show the early change while bounding each library's
# curve. Exact targets and the full-depth endpoint are added independently.
RAREFACTION_GRID_POINTS <- 40L

# Each probability calculation holds at most this many possible group counts.
# Chunking bounds memory without dropping tails or changing the calculation.
RAREFACTION_CHUNK_SIZE <- 65536L

#' Check the positive integer support vector of already defined groups.
#'
#' An empty vector represents an empty repertoire. R stores these counts as
#' doubles; totals above its consecutive-integer range cannot support exact
#' without-replacement count arithmetic.
validate_rarefaction_counts <- function(consensus_count) {
  counts <- as_read_support(consensus_count, "consensus_count")
  if (any(counts <= 0 | counts != floor(counts))) {
    stop("consensus_count must contain positive integer group counts")
  }
  if (sum(counts) > 2^53 - 1) {
    stop("total consensus_count exceeds R's consecutive-integer range")
  }
  counts
}

#' Expected Shannon contribution from one group in a draw without replacement.
#'
#' The count X has a hypergeometric distribution. Sum -(X/m) log(X/m) over
#' its full positive support; X = 0 contributes zero. Equal original group
#' counts share this calculation in interpolate_hill_diversity().
expected_group_shannon <- function(consensus_count, total, depth) {
  first <- max(1, depth - (total - consensus_count))
  last <- min(consensus_count, depth)
  entropy <- 0
  while (first <= last) {
    end <- min(last, first + RAREFACTION_CHUNK_SIZE - 1)
    drawn <- seq(from = first, to = end)
    proportion <- drawn / depth
    probability <- stats::dhyper(drawn, consensus_count,
                                 total - consensus_count, depth)
    entropy <- entropy - sum(probability * proportion * log(proportion))
    first <- end + 1
  }
  entropy
}

#' Interpolate Hill diversity at specified integer support depths.
#'
#' consensus_count is one positive integer per retained group. depths may
#' range from zero to total support. Depth zero has zero richness and missing
#' q = 1 and q = 2. At full depth the values equal observed Hill diversity.
#' No minimum is applied here because curves remain available for diagnostics.
interpolate_hill_diversity <- function(consensus_count, depths) {
  counts <- validate_rarefaction_counts(consensus_count)
  total <- sum(counts)
  if (any(!is.finite(depths) | depths != floor(depths) |
          depths < 0 | depths > total)) {
    stop("depths must be integers between zero and total consensus_count")
  }
  observed <- c(richness = length(counts), hill_shannon = NA_real_,
                inverse_simpson = NA_real_)
  concentration <- NA_real_
  if (total > 0) {
    p <- counts / total
    observed["hill_shannon"] <- exp(-sum(p * log(p)))
    observed["inverse_simpson"] <- 1 / sum(p^2)
  }
  if (total > 1) {
    concentration <- sum((counts / total) * ((counts - 1) / (total - 1)))
  }
  sizes <- sort(unique(counts))
  multiplicity <- tabulate(match(counts, sizes), nbins = length(sizes))
  values <- vapply(depths, function(depth) {
    if (depth == 0) return(c(0, NA_real_, NA_real_))
    if (depth == total) return(unname(observed))
    if (depth == 1 || length(counts) == 1) return(rep(1, 3))
    if (all(counts == 1)) return(rep(depth, 3))
    richness <- sum(multiplicity * stats::phyper(
      0, sizes, total - sizes, depth, lower.tail = FALSE))
    entropy <- sum(multiplicity * vapply(sizes, expected_group_shannon,
                                         numeric(1), total = total,
                                         depth = depth))
    c(richness, exp(entropy), depth / (1 + (depth - 1) * concentration))
  }, numeric(3))
  data.frame(
    support_depth = as.numeric(depths), richness = values[1, ],
    hill_shannon = values[2, ], inverse_simpson = values[3, ],
    row.names = NULL
  )
}

#' Choose a bounded interpolation grid including each available target.
#'
#' An empty repertoire has a single zero-depth endpoint. Targets above total
#' support are omitted; adding a target never extends a curve beyond its data.
rarefaction_depth_grid <- function(consensus_count, targets = numeric(),
                                   grid_points = RAREFACTION_GRID_POINTS) {
  total <- sum(validate_rarefaction_counts(consensus_count))
  if (length(grid_points) != 1L || !is.finite(grid_points) ||
      grid_points < 2 || grid_points != floor(grid_points)) {
    stop("grid_points must be an integer of at least two")
  }
  if (any(!is.finite(targets) | targets < 1 | targets != floor(targets))) {
    stop("rarefaction targets must be positive integers")
  }
  if (total == 0) return(0)
  grid <- if (total <= grid_points) seq_len(total) else
    pmin(total, round(exp(seq(0, log(total), length.out = grid_points))))
  sort(unique(c(1, grid, total, targets[targets <= total])))
}

#' Build a library's diagnostic curve with exact target and endpoint markers.
build_rarefaction_curve <- function(consensus_count, targets = numeric(),
                                    grid_points = RAREFACTION_GRID_POINTS) {
  depths <- rarefaction_depth_grid(consensus_count, targets, grid_points)
  curve <- interpolate_hill_diversity(consensus_count, depths)
  curve$is_target <- curve$support_depth %in% targets
  curve$is_endpoint <- curve$support_depth == sum(consensus_count)
  curve
}

#' Return the model-based q = 2 asymptote and its availability.
#'
#' The asymptote A is the reciprocal of the corrected concentration
#' C = sum(n (n - 1)) / (N (N - 1)). It is an estimated plateau, not an unbiased
#' diversity estimate. interpolate_hill_diversity() gives inverse Simpson at
#' depth m as m / (1 + (m - 1) C) = m / (1 + (m - 1) / A), which reaches 95% of
#' A at the continuous depth 19 (A - 1). That depth can exceed observed support
#' and is not a required FASTQ depth. One observed group carrying at least two
#' support units gives A = 1. Below two support units the asymptote is
#' insufficient_support, and all-singleton support (C = 0) has no finite
#' asymptote and is no_repeated_support.
q2_asymptote <- function(consensus_count) {
  counts <- validate_rarefaction_counts(consensus_count)
  total <- sum(counts)
  concentration <- if (total >= 2) {
    sum((counts / total) * ((counts - 1) / (total - 1)))
  } else NA_real_
  status <- if (total < 2) "insufficient_support" else if (concentration == 0)
    "no_repeated_support" else "available"
  data.frame(
    inverse_simpson_asymptote = if (status == "available") 1 / concentration else NA_real_,
    asymptote_status = status
  )
}

#' Summarise target-depth eligibility and the three rarefied Hill values.
#'
#' The minimum gates Shannon entropy, Hill q = 1 and q = 2, dominance and target
#' values at this locus. It leaves observed richness, support and curves intact.
#' An absent target has `rarefaction_status=not_configured`. With a target,
#' empty and below-minimum statuses precede the below-target check. The output
#' retains total support and the configured target.
summarise_rarefaction <- function(consensus_count, target = NULL,
                                  minimum_consensus_count = 0) {
  counts <- validate_rarefaction_counts(consensus_count)
  total <- sum(counts)
  diversity_status <- repertoire_support_status(total, minimum_consensus_count)
  if (!is.null(target) && (length(target) != 1L || !is.finite(target) ||
                          target < 1 || target != floor(target))) {
    stop("target must be a positive integer or NULL")
  }
  rarefaction_status <- if (is.null(target)) "not_configured" else if (
    diversity_status != "available") diversity_status else if (
      total < target) "below_target_support" else "available"
  values <- if (rarefaction_status == "available") {
    interpolate_hill_diversity(counts, target)
  } else {
    data.frame(richness = NA_real_, hill_shannon = NA_real_,
               inverse_simpson = NA_real_)
  }
  data.frame(
    consensus_count = total, richness = length(counts),
    min_consensus_count_for_diversity = minimum_consensus_count,
    diversity_status = diversity_status,
    rarefaction_target = if (is.null(target)) NA_real_ else target,
    rarefaction_status = rarefaction_status,
    rarefied_richness = values$richness,
    rarefied_hill_shannon = values$hill_shannon,
    rarefied_inverse_simpson = values$inverse_simpson
  )
}

#
# Observed and interpolated diversity from one authoritative support vector.
# The minimum gates Shannon entropy, Hill q = 1 and q = 2 and dominance.
# Observed richness, support and curves remain visible below that minimum.

# The released diversity row: support and its minimum and status first, then
# the observed values, the target-depth values and the q = 2 asymptote. The
# row's total support is grouping-eligible support, so it is released as
# grouped_consensus_count; compute_diversity_row() keeps
# consensus_count, the name of the support vector it sums.
DIVERSITY_DETAIL_COLUMNS <- c(
  "library_id", "locus", "group_unit", "grouped_consensus_count",
  "min_consensus_count_for_diversity", "diversity_status", "richness",
  "shannon_entropy", "hill_shannon", "inverse_simpson", "clonality",
  "top_1_support_fraction", "top_5_support_fraction", "top_10_support_fraction",
  "rarefaction_target", "rarefaction_status", "rarefied_richness",
  "rarefied_hill_shannon", "rarefied_inverse_simpson",
  "inverse_simpson_asymptote", "asymptote_status"
)

RAREFACTION_CURVE_COLUMNS <- c(
  "library_id", "locus", "group_unit", "support_depth", "is_target", "is_endpoint",
  "richness", "hill_shannon", "inverse_simpson"
)

#' Build one detailed row and its optional support-depth curve.
#'
#' groups has one row per retained clone or single-chain clonotype. An absent
#' target disables interpolation for this locus. No grouping is performed here.
#' The detail row follows DIVERSITY_DETAIL_COLUMNS, and hill_shannon there is
#' exp(shannon_entropy).
build_repertoire_diversity <- function(groups, library_id, locus, group_unit,
                                       target = NULL, minimum_consensus_count = 0) {
  counts <- validate_rarefaction_counts(groups$consensus_count)
  observed <- compute_diversity_row(groups, library_id, locus, group_unit)
  observed$hill_shannon <- exp(observed$shannon_entropy)
  summary <- summarise_rarefaction(counts, target, minimum_consensus_count)
  if (summary$diversity_status != "available") {
    for (column in c("shannon_entropy", "hill_shannon", "inverse_simpson",
                     "clonality", "top_1_support_fraction",
                     "top_5_support_fraction", "top_10_support_fraction")) {
      observed[[column]] <- NA_real_
    }
  }
  names(observed)[names(observed) == "consensus_count"] <-
    "grouped_consensus_count"
  summary <- summary[, setdiff(names(summary), c("consensus_count", "richness")),
                     drop = FALSE]
  detail <- cbind(observed, summary, q2_asymptote(counts))
  curve <- if (is.null(target)) {
    make_empty_frame(RAREFACTION_CURVE_COLUMNS)
  } else {
    data.frame(library_id = library_id, locus = locus, group_unit = group_unit,
               build_rarefaction_curve(counts, target), stringsAsFactors = FALSE)
  }
  list(detail = detail[, DIVERSITY_DETAIL_COLUMNS, drop = FALSE],
       curve = curve[, RAREFACTION_CURVE_COLUMNS, drop = FALSE])
}

#' Read an optional per-locus setting from Snakemake's named list.
locus_setting <- function(settings, locus, default = NULL) {
  if (locus %in% names(settings)) as.numeric(settings[[locus]]) else default
}
