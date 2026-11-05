# Apply configured regional coverage requirements to clone-assigned IGH
# measurements, then summarise by clone, library and qualified isotype group.
# Retain counts and explicit unavailable states. Whole-V and regional results
# help assess unequal V coverage. These helpers use base R.

WHOLE_V_REGION <- "whole_v"

# The donor member shard keeps mutation measurements for whole V, CDR2 and FWR3.
# Its eligibility flags and evidence IDs cover all five subregions and whole V,
# as do the member SHM rows, the library SHM table and clone SHM detail.
CLONE_SHM_V_REGIONS <- c(WHOLE_V_REGION, "cdr2", "fwr3")

# The released SHM table, in the order the data dictionary declares it. One row per
# library, isotype group and V region.
SHM_COLUMNS <- c("scope_id", "library_id", "donor_id", "tissue_id", "isotype_group",
                 "v_region", "n_members_measured", "n_members_eligible",
                 "n_members_zero_shm", "n_clones_eligible",
                 "n_clones_zero_median_shm", "fraction_clones_zero_median_shm",
                 "mean_clone_median_shm", "median_clone_median_shm",
                 "q25_clone_median_shm", "q75_clone_median_shm", "iqr_clone_median_shm",
                 "median_trust4_consensus_base_votes")

# The SHM table reports class rows before subclass rows. Unresolved subclass
# calls can contribute to a class, so subclass rows need not sum to that class.
SHM_ISOTYPE_GROUPS <- c(ANALYSIS_IGH_ISOTYPE_GROUPS, ANALYSIS_IGH_SUBCLASS_GROUPS)

#' Name the V regions reported, in reporting order.
#'
#' `region_levels` are the IMGT V sub-regions in 5' to 3' order. Returns them
#' preceded by the whole-V region.
shm_v_regions <- function(region_levels) {
  c(WHOLE_V_REGION, region_levels)
}

#' Which columns of the joined membership describe one V region's mutation load.
#'
#' `v_region` is one of shm_v_regions(). Returns the rate column, the mutation
#' count, the count of positions that contributed to the rate, the count of those
#' that are also vote-supported, the count of positions the germline occupies at
#' all and the consensus-base vote medians. For the whole V region the per-region
#' names are returned instead, because the whole V is summarised as its five
#' sub-regions rather than as a window of its own.
shm_region_columns <- function(v_region, region_levels) {
  regions <- if (identical(v_region, WHOLE_V_REGION)) region_levels else v_region
  list(rate = if (identical(v_region, WHOLE_V_REGION)) "shm_frequency" else
         paste0("shm_frequency_", v_region),
       mutations = paste0("mutations_", regions),
       informative = paste0("informative_", regions),
       supported = paste0("informative_supported_", regions),
       germline_occupied = paste0("germline_occupied_", regions),
       votes = paste0("median_majority_votes_", regions))
}

#' How many mutations each member carries in one V region.
#'
#' `members` is clone membership with the IGH feature columns joined on and
#' `v_region` is one of shm_v_regions(). The five sub-regions partition the whole V
#' window, so the whole-V count is their sum; a sub-region names one column, and the
#' same code covers both.
#'
#' Whether a sequence is unmutated is decided on this count rather than on the rate
#' beside it, so the answer is an integer zero rather than a rate that happens to
#' compare equal to zero. Returns a numeric vector as long as `members` has rows, NA
#' where the member has no count for a sub-region the region is made of.
region_mutation_counts <- function(members, v_region, region_levels) {
  columns <- shm_region_columns(v_region, region_levels)$mutations
  rowSums(matrix(vapply(columns, function(column) as.numeric(members[[column]]),
                        numeric(nrow(members))),
                 nrow = nrow(members), ncol = length(columns)))
}

#' Test one V region's coverage on each source record independently.
#'
#' `members` carries the IGH feature columns, and `v_region` names the region.
#' A source needs a finite mutation rate and a positive germline-occupied
#' denominator. The fraction of vote-supported positions must reach
#' `min_supported_subregion_fraction` in every subregion being tested.
#' Whole V also needs at least `min_germline_occupied_v_positions` occupied
#' positions. Returns one logical value per source record.
shm_region_eligible <- function(members, v_region,
                                min_supported_subregion_fraction,
                                min_germline_occupied_v_positions,
                                region_levels) {
  columns <- shm_region_columns(v_region, region_levels)
  # vapply() returns a bare vector when there is one member row and a matrix
  # otherwise, so the shape is imposed here rather than left to the row count.
  numeric_matrix <- function(names) {
    matrix(vapply(names, function(name) as.numeric(members[[name]]),
                  numeric(nrow(members))),
           nrow = nrow(members), ncol = length(names))
  }
  supported <- numeric_matrix(columns$supported)
  occupied <- numeric_matrix(columns$germline_occupied)
  covered <- is.finite(supported) & is.finite(occupied) & occupied > 0 &
    supported / occupied >= min_supported_subregion_fraction
  eligible <- is.finite(as.numeric(members[[columns$rate]])) &
    rowSums(!covered) == 0L
  if (identical(v_region, WHOLE_V_REGION)) {
    eligible <- eligible &
      rowSums(occupied) >= min_germline_occupied_v_positions
  }
  eligible
}

#' Select one evidence source for each inference member and V region.
#'
#' `members` are the unique-sequence inference members with the IGH feature
#' columns joined on. `candidates` are the retained backbone representatives
#' carrying each member's exact full DNA, with their own features and original
#' `consensus_count`, and `targets` names the member each candidate belongs to.
#' Returns one row per member and region: the selected `evidence_sequence_id`,
#' `eligible`, that source's `median_trust4_consensus_base_votes` (for whole V, the mean of
#' its five sub-region medians), and `n_vote_supported_positions`.
#'
#' For each region every candidate is tested alone against shm_region_eligible()
#' first. Among the candidates that pass, the one with the most supported
#' positions in the region wins, then the one with the larger original record
#' support, then the first `sequence_id` in sort order. Whole V requires one
#' candidate to pass all five sub-regions and the occupancy floor. Vote matrices
#' and partial coverage are never combined across candidates. When no candidate
#' passes, the member is not eligible for that region and its evidence ID and
#' votes are NA.
#'
#' `n_vote_supported_positions` is the supported position count of the candidate
#' this ranking places first, whether or not it passes, because that count
#' explains the coverage decision. It is NA where the count is not finite.
choose_shm_evidence <- function(members, candidates, targets,
                                min_supported_subregion_fraction,
                                min_germline_occupied_v_positions, region_levels) {
  if (length(targets) != nrow(candidates) || anyDuplicated(members$sequence_id) ||
      anyDuplicated(candidates$sequence_id) || anyNA(targets) ||
      !setequal(targets, members$sequence_id)) stop("SHM evidence must cover the merged members")
  if (!nrow(members)) return(data.frame(sequence_id = character(), v_region = character(),
    evidence_sequence_id = character(), eligible = logical(), median_trust4_consensus_base_votes = numeric(),
    n_vote_supported_positions = numeric()))
  support <- as_read_support(candidates$consensus_count)
  target_index <- match(targets, members$sequence_id)
  if (any(candidates$library_id != members$library_id[target_index])) {
    stop("SHM evidence crosses libraries")
  }
  # Equal DNA must have equal sequence-derived measurements and isotype context.
  # Position support can differ because it comes from the original assembly.
  agreement <- c("qualified_c_call", "shm_frequency", "informative_v",
    paste0("shm_frequency_", region_levels), paste0("mutations_", region_levels),
    paste0("informative_", region_levels), paste0("germline_occupied_", region_levels))
  for (column in agreement) {
    actual <- as.character(candidates[[column]])
    expected <- as.character(members[[column]][target_index])
    same <- (is.na(actual) & is.na(expected)) |
      (!is.na(actual) & !is.na(expected) & actual == expected)
    if (any(!same)) stop("Identical-sequence SHM evidence disagrees in ", column)
  }
  frames <- lapply(shm_v_regions(region_levels), function(region) {
    columns <- shm_region_columns(region, region_levels)
    eligible <- shm_region_eligible(candidates, region,
      min_supported_subregion_fraction, min_germline_occupied_v_positions, region_levels)
    matrix_columns <- function(fields) {
      matrix(vapply(fields, function(field) as.numeric(candidates[[field]]),
                    numeric(nrow(candidates))), nrow = nrow(candidates), ncol = length(fields))
    }
    coverage <- rowSums(matrix_columns(columns$supported))
    votes <- rowMeans(matrix_columns(columns$votes))
    coverage[!is.finite(coverage)] <- -Inf
    ranked <- order(targets, -as.integer(eligible), -coverage, -support, candidates$sequence_id)
    winners <- ranked[!duplicated(targets[ranked])]
    chosen <- winners[match(members$sequence_id, targets[winners])]
    passes <- eligible[chosen]
    data.frame(sequence_id = as.character(members$sequence_id), v_region = rep(region, nrow(members)),
      evidence_sequence_id = ifelse(passes, as.character(candidates$sequence_id[chosen]), NA_character_),
      eligible = passes, median_trust4_consensus_base_votes = ifelse(passes, votes[chosen], NA_real_),
      n_vote_supported_positions = ifelse(is.finite(coverage[chosen]), coverage[chosen], NA_real_),
      stringsAsFactors = FALSE)
  })
  do.call(rbind, frames)
}

# With no `evidence`, each member is its own only candidate source, as in direct
# calls from the unit suites. The donor stage always supplies the shared
# multi-source decisions, so all outputs use the same selected evidence.
member_region_evidence <- function(members, region, min_subregion_fraction,
                                    min_germline_occupied_v, region_levels,
                                    evidence = NULL) {
  if (is.null(evidence)) {
    evidence <- choose_shm_evidence(members, members, members$sequence_id,
      min_subregion_fraction, min_germline_occupied_v, region_levels)
  }
  selected <- evidence[evidence$v_region == region, , drop = FALSE]
  if (anyDuplicated(selected$sequence_id) || !setequal(selected$sequence_id, members$sequence_id)) {
    stop("SHM evidence differs from the member population")
  }
  selected[match(members$sequence_id, selected$sequence_id), , drop = FALSE]
}

#' Return the donor member shard's mutation columns in column order.
#'
#' The whole-V rate and count come first, then each subregion of
#' CLONE_SHM_V_REGIONS. Eligibility and evidence IDs follow for all five
#' subregions and whole V. The release reports member SHM from
#' build_clone_member_shm_table(), so no column here enters clone_members.tsv.
#'
#' Subregion rates keep their names from the IGH feature table. The whole-V
#' rate and count take a `_whole_v` suffix, as every other whole-V member column
#' does: `shm_frequency_whole_v` is the feature table's `shm_frequency`, and
#' `mutations_whole_v` sums the five subregion counts, matching the whole-V
#' denominator `informative_v`.
member_shm_columns <- function() {
  subregions <- setdiff(CLONE_SHM_V_REGIONS, WHOLE_V_REGION)
  regions <- c(WHOLE_V_REGION, "fwr1", "cdr1", "fwr2", "cdr2", "fwr3")
  c("shm_frequency_whole_v", "mutations_whole_v",
    paste0("shm_frequency_", subregions),
    paste0("mutations_", subregions),
    paste0("shm_eligible_", regions), paste0("shm_evidence_sequence_id_", regions))
}

#' Attach the donor member shard's mutation columns.
#'
#' `members` already carries the IGH features joined by
#' attach_igh_features_to_membership(). This function adds mutation counts,
#' eligibility flags and evidence IDs from the same source decisions used by
#' `shm.tsv`. Source coverage counts stay in the per-library feature tables.
#' Mutation rates remain available before coverage eligibility in the donor
#' shard, so a finite rate can accompany a false eligibility flag there.
#'
#' Returns `members` with every column named by member_shm_columns().
add_member_shm_columns <- function(members, min_supported_subregion_fraction,
                                   min_germline_occupied_v_positions,
                                   region_levels, evidence = NULL) {
  for (v_region in CLONE_SHM_V_REGIONS) {
    members[[paste0("mutations_", v_region)]] <- region_mutation_counts(
      members, v_region, region_levels)
  }
  members$shm_frequency_whole_v <- members$shm_frequency
  for (region in shm_v_regions(region_levels)) {
    selected <- member_region_evidence(members, region, min_supported_subregion_fraction,
      min_germline_occupied_v_positions, region_levels, evidence)
    members[[paste0("shm_eligible_", region)]] <- selected$eligible
    members[[paste0("shm_evidence_sequence_id_", region)]] <- selected$evidence_sequence_id
  }
  members
}

# The released member SHM table, in the order the data dictionary declares it.
# One row per IGH member and V region with a measured rate.
CLONE_MEMBER_SHM_COLUMNS <- c(
  "scope_id", "library_id", "tissue_id", "donor_id", "clone_id", "sequence_id",
  "v_region", "shm_eligible", "n_informative_positions",
  "n_germline_occupied_positions", "n_vote_supported_positions", "n_mutations",
  "shm_frequency", "median_trust4_consensus_base_votes"
)

#' Build the donor shard of clone_member_shm.tsv from enriched member records.
#'
#' `records` are the inference members with the IGH feature columns joined on and
#' their scope, library, tissue and donor. `evidence` supplies the shared per-region
#' source decisions, so `shm_eligible` is the decision shm.tsv uses. A row exists
#' where the member has a finite rate in the region, so the rows of one library,
#' isotype group and region number `n_members_measured` in shm.tsv.
#'
#' The position counts stay filled for ineligible members, because they explain the
#' coverage decision. `n_informative_positions` is the rate's denominator. Whole V
#' takes `informative_v` and sums `germline_occupied_` over the five subregions.
#' `n_mutations` and `shm_frequency` keep their measured values in this shard, and
#' scope collection empties them where `shm_eligible` is not TRUE.
#'
#' Returns the rows in CLONE_MEMBER_SHM_COLUMNS order, with V regions in the order
#' of `v_regions`, or the released header alone when there are no members.
build_clone_member_shm_table <- function(records, min_subregion_fraction,
                                         min_germline_occupied_v, region_levels,
                                         v_regions = shm_v_regions(region_levels),
                                         evidence = NULL) {
  if (!nrow(records)) return(make_empty_frame(CLONE_MEMBER_SHM_COLUMNS))
  rows <- lapply(v_regions, function(region) {
    columns <- shm_region_columns(region, region_levels)
    selected <- member_region_evidence(records, region, min_subregion_fraction,
                                       min_germline_occupied_v, region_levels, evidence)
    rate <- as.numeric(records[[columns$rate]])
    informative <- if (identical(region, WHOLE_V_REGION)) "informative_v" else columns$informative
    occupied <- matrix(vapply(columns$germline_occupied, function(column) as.numeric(records[[column]]),
                              numeric(nrow(records))),
                       nrow = nrow(records), ncol = length(columns$germline_occupied))
    data.frame(
      scope_id = as.character(records$scope_id),
      library_id = as.character(records$library_id),
      tissue_id = as.character(records$tissue_id),
      donor_id = as.character(records$donor_id),
      clone_id = as.character(records$clone_id),
      sequence_id = as.character(records$sequence_id),
      v_region = rep(region, nrow(records)),
      shm_eligible = selected$eligible,
      n_informative_positions = as.numeric(records[[informative]]),
      n_germline_occupied_positions = rowSums(occupied),
      n_vote_supported_positions = selected$n_vote_supported_positions,
      n_mutations = region_mutation_counts(records, region, region_levels),
      shm_frequency = rate,
      median_trust4_consensus_base_votes = selected$median_trust4_consensus_base_votes,
      stringsAsFactors = FALSE
    )[is.finite(rate), , drop = FALSE]
  })
  table <- do.call(rbind, rows)
  rownames(table) <- NULL
  table[, CLONE_MEMBER_SHM_COLUMNS, drop = FALSE]
}

#' Select the clone members that can contribute to one region's SHM summary.
#'
#' `membership` carries the joined IGH features. `evidence` supplies the shared
#' per-region source decisions from choose_shm_evidence(). When it is NULL,
#' member_region_evidence() tests each member as its only source.
#'
#' Returns eligible members with library, sequence and clone IDs, the qualified
#' C call, regional mutation rate and count, selected-source vote median and
#' resolved isotype class and subclass. Class and subclass resolution can differ.
#' The `selected_rows` attribute identifies the retained input rows. The
#' `measured` attribute keeps every input member's library, isotype resolutions
#' and finite-rate flag before coverage eligibility.
select_shm_eligible_members <- function(membership, v_region,
                                        min_supported_subregion_fraction,
                                        min_germline_occupied_v_positions,
                                        region_levels, evidence = NULL) {
  columns <- shm_region_columns(v_region, region_levels)
  selected_evidence <- member_region_evidence(membership, v_region,
                              min_supported_subregion_fraction,
                              min_germline_occupied_v_positions, region_levels, evidence)
  keep <- selected_evidence$eligible
  # The whole V region's vote medians average over the five sub-regions that
  # partition it, and a sub-region names one column, so the same code covers both.
  # Its mutation count is region_mutation_counts(), which the sequence-level tables
  # report from the same definition.
  eligible <- data.frame(
    library_id = as.character(membership$library_id[keep]),
    sequence_id = as.character(membership$sequence_id[keep]),
    clone_id = as.character(membership$clone_id[keep]),
    qualified_c_call = as.character(membership$qualified_c_call[keep]),
    shm_frequency = as.numeric(membership[[columns$rate]][keep]),
    mutations = region_mutation_counts(membership, v_region, region_levels)[keep],
    median_trust4_consensus_base_votes = selected_evidence$median_trust4_consensus_base_votes[keep],
    stringsAsFactors = FALSE
  )
  eligible$isotype_class <- resolve_analysis_igh_isotype(eligible$qualified_c_call, "class")
  eligible$isotype_subclass <- resolve_analysis_igh_isotype(eligible$qualified_c_call,
                                                            "subclass")
  # The census carries the same two resolutions the eligible rows do, rather than
  # the raw call, so shm_group_mask() decides a group's membership the one way for
  # both and no consumer has to resolve a call a second time.
  attr(eligible, "selected_rows") <- which(keep)
  attr(eligible, "measured") <- data.frame(
    library_id = as.character(membership$library_id),
    isotype_class = resolve_analysis_igh_isotype(membership$qualified_c_call, "class"),
    isotype_subclass = resolve_analysis_igh_isotype(membership$qualified_c_call, "subclass"),
    measured = is.finite(as.numeric(membership[[columns$rate]])),
    stringsAsFactors = FALSE
  )
  eligible
}

#' Which members belong to one released isotype group.
#'
#' `members` are rows carrying the two isotype resolutions
#' select_shm_eligible_members() attached (either the eligible members or the
#' pre-eligibility census), and `isotype_group` is one of SHM_ISOTYPE_GROUPS.
#' Returns a logical mask.
#'
#' Every member belongs to `all_IGH`, regardless of its C call. Other groups
#' require resolution at their class or subclass level. A call that resolves to
#' IgG can remain unresolved at subclass level. Summary rows and their measured
#' denominators use this same mask.
shm_group_mask <- function(members, isotype_group) {
  if (identical(isotype_group, "all_IGH")) return(rep(TRUE, nrow(members)))
  column <- if (identical(analysis_igh_isotype_level(isotype_group), "class")) {
    "isotype_class"
  } else {
    "isotype_subclass"
  }
  resolved <- members[[column]]
  !is.na(resolved) & resolved == isotype_group
}

#' The released SHM statistic for one library, isotype group and V region.
#'
#' Each clone contributes the median mutation frequency of its eligible members.
#' The mean, median and quartiles give each clone equal weight, and
#' `iqr_clone_median_shm` is `q75_clone_median_shm` minus `q25_clone_median_shm`.
#' A clone can contribute to several isotype groups when its eligible members
#' represent those groups.
#'
#' `n_clones_zero_median_shm` counts clones whose median regional mutation count
#' is zero. `n_members_zero_shm` counts eligible members with zero mutations.
#' Their denominators are `n_clones_eligible` and `n_members_eligible`. A zero
#' median does not require every member to have zero mutations. Read both counts
#' with `n_members_measured`, which includes members before coverage eligibility.
#' These counts describe recovered receptors, not cells.
#'
#' `members` are one library's eligible members of one isotype group for one
#' region, already selected by build_scope_shm_table() from
#' select_shm_eligible_members() and shm_group_mask(). `library_id` and
#' `isotype_group` name what was selected and `v_region` the region, all three
#' of which stamp the row rather than filter it; `n_measured` is that cell's
#' pre-eligibility count from count_measured_members().
#'
#' Returns one row. Every statistic is undefined below one clone rather than
#' left to the arithmetic, because mean() and quantile() of nothing are not the
#' released null.
summarise_lineage_median_shm <- function(members, library_id, isotype_group, v_region,
                                         n_measured = NA_integer_) {
  # One factor over the members serves both statistics, and both read a vector
  # rather than the frame, because this runs once per library, group and region.
  lineage <- factor(as.character(members$clone_id))
  mutations <- as.numeric(members$mutations)
  lineage_medians <- vapply(split(as.numeric(members$shm_frequency), lineage),
                            stats::median, numeric(1), USE.NAMES = FALSE)
  # A clone enters n_clones_zero_median_shm when its members' median mutation count
  # is zero. shm_frequency is zero exactly when the mutation count is, so this
  # agrees with testing each clone's entry of lineage_medians against zero; it is
  # taken on the integer count so that no floating-point rate is compared against
  # zero. The sequence-level count below is the same test without the median, so
  # the two denominators in the row (n_clones_eligible and n_members_eligible) each have their own
  # zero-mutation numerator.
  clone_zero_median <- vapply(split(mutations, lineage),
                              function(x) stats::median(x) == 0, logical(1),
                              USE.NAMES = FALSE)
  quantile_of <- function(p) {
    if (!length(lineage_medians)) return(NA_real_)
    unname(stats::quantile(lineage_medians, probs = p, names = FALSE))
  }
  q25 <- quantile_of(0.25)
  q75 <- quantile_of(0.75)
  n_sequences <- length(lineage)
  data.frame(
    library_id = library_id,
    isotype_group = isotype_group,
    v_region = v_region,
    n_members_measured = n_measured,
    n_members_eligible = n_sequences,
    n_members_zero_shm = sum(mutations == 0),
    n_clones_eligible = length(lineage_medians),
    n_clones_zero_median_shm = sum(clone_zero_median),
    fraction_clones_zero_median_shm = if (length(lineage_medians)) {
      mean(clone_zero_median)
    } else NA_real_,
    mean_clone_median_shm = if (length(lineage_medians) > 0L) {
      mean(lineage_medians)
    } else {
      NA_real_
    },
    median_clone_median_shm = if (length(lineage_medians) > 0L) {
      stats::median(lineage_medians)
    } else {
      NA_real_
    },
    q25_clone_median_shm = q25,
    q75_clone_median_shm = q75,
    iqr_clone_median_shm = q75 - q25,
    # Votes are not read counts. TRUST4 counts a merged pair twice and writes a
    # synthetic 1 where a tally came out zero, so this is a technical covariate to
    # check the rate against, not a depth. With min_supported_subregion_fraction
    # at 1, every position of an eligible region clears the support floor, so this
    # describes how far above it the surviving evidence sits.
    median_trust4_consensus_base_votes = if (n_sequences > 0L) {
      as.numeric(stats::median(as.numeric(members$median_trust4_consensus_base_votes),
                               na.rm = TRUE))
    } else {
      NA_real_
    },
    stringsAsFactors = FALSE
  )
}

#' How many members of each library had a rate before eligibility, in one group.
#'
#' `measured` is the census attached by select_shm_eligible_members().
#' `isotype_group` names the group, and `library_ids` lists every scope library.
#' Returns one named count per library in that order, or NA when the census is
#' absent. Compare these counts with `n_members_eligible` to assess coverage retention.
count_measured_members <- function(measured, isotype_group, library_ids) {
  if (is.null(measured)) {
    return(stats::setNames(rep(NA_integer_, length(library_ids)), library_ids))
  }
  counted <- measured$measured & shm_group_mask(measured, isotype_group)
  stats::setNames(
    tabulate(factor(measured$library_id[counted], levels = library_ids),
             nbins = length(library_ids)),
    library_ids)
}

#' Build one scope's released SHM table.
#'
#' `membership` is the scope's clone membership with the IGH feature columns joined
#' on, `manifest` names every library the scope reports with its donor and tissue, and
#' `scope_id` stamps each row. `v_regions` fixes both which regions are reported and
#' their order, `isotype_groups` the same for isotype. `evidence` supplies the shared
#' per-region source decisions. The per-library feature tables retain rates and
#' coverage counts before eligibility for auditing or a different coverage rule.
#'
#' Every library in the manifest gets a row per group per region, including one that
#' contributed no eligible sequence. That result does not establish biological
#' absence. Returns the table in SHM_COLUMNS
#' order, sorted by library, then group, then region.
build_scope_shm_table <- function(membership, manifest, scope_id, v_regions,
                                 min_supported_subregion_fraction,
                                 min_germline_occupied_v_positions,
                                 region_levels,
                                 isotype_groups = SHM_ISOTYPE_GROUPS, audit_counts = FALSE,
                                 evidence = NULL) {
  library_ids <- as.character(manifest$library_id)
  # The libraries index both the members and the census, and the manifest is one
  # row per library, so a repeat is a broken manifest rather than something to
  # summarise twice.
  if (anyDuplicated(library_ids)) {
    stop("scope manifest names a library more than once: ",
         paste(unique(library_ids[duplicated(library_ids)]), collapse = ", "))
  }
  # Eligibility is a property of the region and group membership a property of the
  # C call, so each is decided once for every member and then read per library
  # rather than re-decided inside each row's summary. Deciding them per row would
  # re-resolve every member's isotype for each of the ten groups and six regions in
  # turn.
  audit <- list()
  region_rows <- lapply(v_regions, function(v_region) {
    eligible <- select_shm_eligible_members(membership, v_region,
                                            min_supported_subregion_fraction,
                                            min_germline_occupied_v_positions,
                                            region_levels, evidence)
    measured <- attr(eligible, "measured")
    if (audit_counts) {
      support <- as_read_support(membership$consensus_count)
      selected <- attr(eligible, "selected_rows")
      for (group in isotype_groups) {
        group_mask <- shm_group_mask(measured, group)
        for (i in seq_along(library_ids)) {
          before <- which(membership$library_id == library_ids[i] & group_mask)
          after <- intersect(before, selected)
          context <- list(scope_id = scope_id, donor_id = manifest$donor_id[i],
                          tissue_id = manifest$tissue_id[i], isotype_group = group, v_region = v_region)
          audit[[length(audit) + 1L]] <<- do.call(qc_audit_rows, c(list(
            library_id = library_ids[i], stage = "shm_eligible_members", n_in = length(before),
            n_out = length(after), support_in = sum(support[before]), support_out = sum(support[after]),
            analysis = "shm", stage_index = 60, kind = "eligibility",
            input_unit = "clone_member", support_conserved = FALSE), context))
          audit[[length(audit) + 1L]] <<- do.call(qc_audit_rows, c(list(
            library_id = library_ids[i], stage = "shm_eligible_clones",
            n_in = length(unique(membership$clone_id[after])),
            support_in = sum(support[after]), analysis = "shm", stage_index = 61,
            input_unit = "library_clone_observation"), context))
        }
      }
    }
    do.call(rbind, lapply(isotype_groups, function(group) {
      in_group <- eligible[shm_group_mask(eligible, group), , drop = FALSE]
      by_library <- split(seq_len(nrow(in_group)),
                          factor(in_group$library_id, levels = library_ids))
      n_measured <- count_measured_members(measured, group, library_ids)
      do.call(rbind, lapply(library_ids, function(library_id) {
        summarise_lineage_median_shm(
          in_group[by_library[[library_id]], , drop = FALSE],
          library_id, group, v_region, n_measured[[library_id]])
      }))
    }))
  })
  shm <- do.call(rbind, region_rows)
  shm <- merge(manifest, shm, by = "library_id", all.y = TRUE, sort = FALSE)
  shm$scope_id <- scope_id
  shm <- shm[, SHM_COLUMNS]
  shm <- shm[order(shm$library_id, match(shm$isotype_group, isotype_groups),
                   match(shm$v_region, v_regions)), , drop = FALSE]
  if (audit_counts) attr(shm, "qc_audit") <- do.call(rbind, audit)
  shm
}

#
# Library-local clone support and the mutation states of its observed members.
#
# Each row describes one clone within a library, V region and isotype subset.
# n_members and consensus_count include all members of that subset. Mutation
# counts, median frequency and eligible_consensus_count use only members whose
# selected SHM evidence source passes shm_region_eligible(), the same rule the
# library SHM summary uses. A high-support clone with no evaluable SHM remains in the donor shard
# with its support and member count; scope collection leaves that row out of the
# released clone_shm.tsv, and group_abundance.tsv keeps the clone's support.
#
# Mutation states describe zero or positive scored mutations in the specified
# region. They do not establish exact germline identity or the history of a
# biological clone. A mixture can have median SHM zero. Isotype subsets overlap;
# mixed-isotype clones receive the applicable subset rows, not one chosen label.
# This analysis has no diversity support minimum.

CLONE_SHM_COLUMNS <- c(
  "scope_id", "donor_id", "tissue_id", "library_id", "clone_id", "isotype_group",
  "v_region", "n_members", "n_members_eligible", "n_members_zero_shm",
  "n_members_positive_shm", "median_shm_frequency", "mutation_state",
  "consensus_count", "eligible_consensus_count"
)

#' Build clone SHM detail from enriched donor records before member projection.
#'
#' records uses the existing enriched membership schema, including library_id,
#' clone_id, consensus_count and raw per-region feature counts. The output uses
#' the released CLONE_SHM_COLUMNS names. The default isotype groups match the
#' library SHM summaries; unresolved calls and IgD/IgE remain in all_IGH. A
#' subset with no observed member produces no row, while observed but
#' unevaluable subsets remain here; scope collection leaves those rows out of
#' the released table. Rows are sorted by scope, donor and library,
#' then by isotype group and V region in the order of `isotype_groups` and
#' `v_regions`, then by clone.
build_clone_shm_table <- function(records, min_subregion_fraction,
                                   min_germline_occupied_v, region_levels,
                                   v_regions = shm_v_regions(region_levels),
                                   isotype_groups = SHM_ISOTYPE_GROUPS, evidence = NULL) {
  require_columns(records, c("scope_id", "donor_id", "tissue_id", "library_id",
                             "clone_id", "qualified_c_call", "consensus_count"),
                  "Enriched IGH records")
  if (!nrow(records)) return(make_empty_frame(CLONE_SHM_COLUMNS))
  support <- as_read_support(records$consensus_count, "consensus_count")
  isotypes <- data.frame(
    isotype_class = resolve_analysis_igh_isotype(records$qualified_c_call, "class"),
    isotype_subclass = resolve_analysis_igh_isotype(records$qualified_c_call, "subclass")
  )
  # The key is never parsed. The separator cannot occur in manifest identifiers.
  group_key <- paste(records$scope_id, records$donor_id, records$library_id,
                     records$clone_id, sep = "\x1f")
  region_data <- lapply(v_regions, function(region) {
    list(
      eligible = member_region_evidence(records, region, min_subregion_fraction,
                                       min_germline_occupied_v, region_levels, evidence)$eligible,
      mutations = region_mutation_counts(records, region, region_levels),
      frequency = as.numeric(records[[shm_region_columns(region, region_levels)$rate]])
    )
  })
  tables <- list()
  for (isotype_group in isotype_groups) {
    rows <- which(shm_group_mask(isotypes, isotype_group))
    units <- split(rows, group_key[rows])
    if (!length(units)) next
    first <- vapply(units, function(unit) unit[1], integer(1), USE.NAMES = FALSE)
    context <- data.frame(
      scope_id = as.character(records$scope_id[first]),
      donor_id = as.character(records$donor_id[first]),
      tissue_id = as.character(records$tissue_id[first]),
      library_id = as.character(records$library_id[first]),
      clone_id = as.character(records$clone_id[first]),
      isotype_group = isotype_group,
      consensus_count = vapply(units, function(unit) sum(support[unit]),
                                numeric(1), USE.NAMES = FALSE),
      n_members = lengths(units, use.names = FALSE),
      stringsAsFactors = FALSE
    )
    for (index in seq_along(v_regions)) {
      region <- region_data[[index]]
      summaries <- t(vapply(units, function(unit) {
        evaluable <- unit[region$eligible[unit]]
        if (!length(evaluable)) return(c(0, 0, 0, NA_real_, 0))
        mutations <- region$mutations[evaluable]
        c(length(evaluable), sum(mutations == 0), sum(mutations > 0),
          stats::median(region$frequency[evaluable]), sum(support[evaluable]))
      }, numeric(5), USE.NAMES = FALSE))
      state <- ifelse(summaries[, 1] == 0, "no_evaluable_members",
                       ifelse(summaries[, 2] == summaries[, 1], "all_zero",
                         ifelse(summaries[, 3] == summaries[, 1], "all_positive", "mixed")))
      tables[[length(tables) + 1L]] <- data.frame(
        context,

        v_region = sub("^whole_v$", "whole_v", v_regions[index]),
        n_members_eligible = as.integer(summaries[, 1]),
        n_members_zero_shm = as.integer(summaries[, 2]),
        n_members_positive_shm = as.integer(summaries[, 3]),
        median_shm_frequency = summaries[, 4], mutation_state = state,
        eligible_consensus_count = summaries[, 5],
        stringsAsFactors = FALSE
      )
    }
  }
  if (!length(tables)) return(make_empty_frame(CLONE_SHM_COLUMNS))
  detail <- do.call(rbind, tables)
  detail <- detail[order(detail$scope_id, detail$donor_id, detail$library_id,
                         match(detail$isotype_group, isotype_groups),
                         match(detail$v_region, v_regions), detail$clone_id),
                    CLONE_SHM_COLUMNS, drop = FALSE]
  rownames(detail) <- NULL
  detail
}

# Keep every member-region combination, including regions with no measured rate.
# The evidence table is the decision returned by choose_shm_evidence().
MEMBER_SHM_COVERAGE_COLUMNS <- c(
  "scope_id", "library_id", "tissue_id", "donor_id", "clone_id", "sequence_id",
  "v_region", "measurement_status", "shm_eligible", "evidence_sequence_id",
  "germline_reconstruction_status")

build_member_shm_coverage <- function(records, evidence, region_levels) {
  if (!nrow(records)) return(make_empty_frame(MEMBER_SHM_COVERAGE_COLUMNS))
  rows <- lapply(shm_v_regions(region_levels), function(region) {
    chosen <- evidence[evidence$v_region == region, , drop = FALSE]
    chosen <- chosen[match(records$sequence_id, chosen$sequence_id), , drop = FALSE]
    if (anyNA(chosen$sequence_id)) stop("SHM evidence does not cover every member")
    measured <- is.finite(as.numeric(records[[shm_region_columns(region, region_levels)$rate]]))
    if (any(chosen$eligible & !measured)) stop("SHM eligibility requires a measured rate")
    table <- records[, c("scope_id", "library_id", "tissue_id", "donor_id", "clone_id", "sequence_id"), drop = FALSE]
    table$v_region <- region
    table$measurement_status <- ifelse(!measured, "unmeasured",
      ifelse(chosen$eligible, "eligible", "measured_ineligible"))
    table$shm_eligible <- chosen$eligible
    table$evidence_sequence_id <- chosen$evidence_sequence_id
    germline <- as.character(records$germline_alignment_d_mask)
    table$germline_reconstruction_status <- ifelse(!is.na(germline) & nzchar(germline),
      "available", "unavailable_other")
    table
  })
  do.call(rbind, rows)[, MEMBER_SHM_COVERAGE_COLUMNS, drop = FALSE]
}
