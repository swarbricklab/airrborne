# Measure V-region mutations, opportunities and position support for every
# common-QC IGH record with SHazaM 1.3.2 in the pinned image, through the
# shm_regions.R definitions. Copy qualified_c_call from the supplied sidecar.
# Regional coverage gates belong to summary consumers;
# these measurements preserve records and release their denominators.

#' The columns this stage releases, in the order it releases them.
#'
#' Column order is 5' to 3' within each family, so the table reads in sequence
#' order. The stage writes exactly these columns in exactly this order, which is
#' the contract the release validator checks against the data dictionary, and
#' which test/components/shm/unit/test-features.R asserts against this same
#' definition rather than against a copy of it.
#'
#' There is no whole-V mutation count, because it is the exact integer sum of the
#' five sub-region counts and a derivable column is one more thing that can
#' disagree with itself.
#'
#' `region_levels` supplies the sub-region labels when this function is called.
#' Callers using the default must load shm_regions.R before calling it.
igh_feature_columns <- function(region_levels = V_REGION_LEVELS) {
  c(
    "library_id", "sequence_id", "qualified_c_call",
    "v_call_status", "v_germline_allele",
    "shm_frequency", "informative_v",
    paste0("shm_frequency_", region_levels),
    paste0("mutations_", region_levels),
    paste0("informative_", region_levels),
    paste0("germline_occupied_", region_levels),
    paste0("median_majority_votes_", region_levels),
    paste0("informative_supported_", region_levels)
  )
}

#' Measure each record's V-region mutation load and vote support, whole and by
#' sub-region.
#'
#' `records` are one library's common-QC IGH rows, `votes` the library's IGH
#' consensus-vote table keyed on the same `sequence_id`, and `region_of_position`
#' the IMGT sub-region partition, taken once so the loop does not rebuild it per
#' record. Returns one row per record in igh_feature_columns() order, with `qualified_c_call`
#' copied from the qualified isotype-call sidecar joined by the caller.
#'
#' A rate requires at least one informative position. Numerators, occupancy
#' and support let summary consumers apply their stated coverage requirements.
#'
#' CreateGermlines uses one V allele. v_germline_allele records its identity
#' and v_call_status retains call ambiguity. Alternative tied alleles can
#' place mismatches differently within V, affecting regional measurements.
#'
#' SHazaM method documentation:
#' https://shazam.readthedocs.io/en/stable/topics/observedMutations/
build_igh_v_region_feature_table <- function(records, votes, min_position_votes,
                                             region_of_position =
                                               imgt_v_region_of_position()) {
  n_records <- nrow(records)
  n_positions <- length(region_of_position)
  vote_index <- match_vote_rows_to_records(records, votes)

  features <- data.frame(
    library_id = as.character(records$library_id),
    sequence_id = as.character(records$sequence_id),
    qualified_c_call = as.character(records$qualified_c_call),
    v_call_status = vj_call_ambiguity_status(records$v_call),
    v_germline_allele = first_vj_allele(records$v_call),
    shm_frequency = rep(NA_real_, n_records),
    informative_v = rep(NA_integer_, n_records),
    stringsAsFactors = FALSE
  )

  # One SHazaM call per record, over the sub-region partition. The whole-V rate is
  # then the summed counts divided by the summed opportunity, which reproduces
  # shazam::observedMutations(regionDefinition = IMGT_V) exactly, because the two
  # region definitions describe the same 312-position window.
  shm_rows <- which(!is_airr_null(records$sequence_alignment) &
                      !is_airr_null(records$germline_alignment_d_mask))
  # A record with no alignment pair is left absent in every one of these, rather
  # than zero: nothing was measured, which is not the same as measuring nothing.
  # A record whose votes are all absent does reach the loop, and comes out with a
  # supported count of zero, which is the honest answer for it.
  empty_region_matrix <- function(value = NA_integer_) {
    matrix(value, nrow = n_records, ncol = length(V_REGION_LEVELS),
           dimnames = list(NULL, V_REGION_LEVELS))
  }
  mutations <- empty_region_matrix()
  informative <- empty_region_matrix()
  occupied <- empty_region_matrix()
  median_votes <- empty_region_matrix(NA_real_)
  supported <- empty_region_matrix()

  for (row in shm_rows) {
    observed <- summarise_observed_mutations_by_region(
      records$sequence_alignment[row], records$germline_alignment_d_mask[row]
    )
    mutations[row, ] <- observed$mutations[V_REGION_LEVELS]
    informative[row, ] <- observed$informative[V_REGION_LEVELS]
    occupied[row, ] <- count_germline_occupied_by_region(
      records$germline_alignment_d_mask[row], region_of_position
    )[V_REGION_LEVELS]

    total_votes <- parse_vote_vector(
      votes$total_votes_by_imgt_position[vote_index[row]], n_positions)
    dissenting_votes <- parse_vote_vector(
      votes$non_consensus_votes_by_imgt_position[vote_index[row]], n_positions)
    support <- summarise_vote_support_by_region(
      position_majority_votes(total_votes, dissenting_votes),
      imgt_v_informative_positions(records$sequence_alignment[row],
                                   records$germline_alignment_d_mask[row],
                                   n_positions),
      min_position_votes, region_of_position
    )
    median_votes[row, ] <- support$median[V_REGION_LEVELS]
    supported[row, ] <- support$n_supported[V_REGION_LEVELS]
  }

  features$informative_v <- as.integer(rowSums(informative))
  features$shm_frequency <- ifelse(features$informative_v > 0L,
                             rowSums(mutations) / features$informative_v, NA_real_)

  for (region in V_REGION_LEVELS) {
    features[[paste0("shm_frequency_", region)]] <-
      as.double(ifelse(informative[, region] > 0L,
                       mutations[, region] / informative[, region], NA_real_))
    features[[paste0("mutations_", region)]] <- as.integer(mutations[, region])
    features[[paste0("informative_", region)]] <- as.integer(informative[, region])
    features[[paste0("germline_occupied_", region)]] <- as.integer(occupied[, region])
    features[[paste0("median_majority_votes_", region)]] <-
      as.double(median_votes[, region])
    features[[paste0("informative_supported_", region)]] <-
      as.integer(supported[, region])
  }
  # The loop above adds five columns per sub-region, so the frame comes out region
  # by region. Project it onto the released order here rather than leaving a caller
  # to discover that the table and its own column contract disagree.
  features[, igh_feature_columns(), drop = FALSE]
}

#' Find each IGH record's row in the consensus-vote table.
#'
#' `records` are the library's common-QC IGH rows and `votes` its released
#' consensus-vote table. Returns one row index of `votes` per row of `records`.
#'
#' Both tables are built from the same filtered AIRR file, so the vote table holds
#' exactly one row per IGH record. A missing one means the two are out of step: a
#' stale reannotation paired with current votes runs clean and produces
#' plausible-looking numbers, so the mismatch has to be caught here rather than
#' noticed later.
match_vote_rows_to_records <- function(records, votes) {
  require_columns(votes, c("sequence_id", "total_votes_by_imgt_position",
                           "non_consensus_votes_by_imgt_position"),
                  "IGH consensus votes")
  if (anyDuplicated(as.character(votes$sequence_id))) {
    stop("sequence_id is not unique in the IGH consensus-vote table")
  }
  index <- match(as.character(records$sequence_id), as.character(votes$sequence_id))
  if (anyNA(index)) {
    stop(sum(is.na(index)), " IGH record(s) have no consensus-vote row; the votes ",
         "and the reannotation are out of step")
  }
  index
}

#' Report how much of each IMGT V sub-region this library actually measured.
#'
#' `features` is the table build_igh_v_region_feature_table() returned. Returns one
#' row per sub-region: the median number of positions that contributed a mutation
#' opportunity, the median number of those that are also vote-supported, the median
#' number the germline occupies at all, the median vote support behind the
#' contributing positions, and how many records were given a rate there.
#'
#' The three medians read together are the diagnostic. A sub-region whose median
#' coverage is well below its median occupancy was measured on a fraction of the
#' region. A lower median supported count also reflects positions with insufficient
#' or missing votes. The SHM gate applies coverage per source record, so these
#' separate medians do not determine how many records qualify. Counts are reported
#' rather than ratios because the denominator differs per record.
summarise_igh_feature_callability <- function(features) {
  # median() of an empty or all-NA integer column returns an integer NA, so each
  # value is coerced before vapply() is asked for a double.
  region_median <- function(prefix) {
    vapply(V_REGION_LEVELS, function(region) {
      as.numeric(stats::median(features[[paste0(prefix, region)]], na.rm = TRUE))
    }, numeric(1), USE.NAMES = FALSE)
  }
  data.frame(
    region = V_REGION_LEVELS,
    median_informative = region_median("informative_"),
    median_informative_supported = region_median("informative_supported_"),
    median_germline_occupied = region_median("germline_occupied_"),
    median_trust4_consensus_base_votes = region_median("median_majority_votes_"),
    n_rate_called = vapply(V_REGION_LEVELS, function(region) {
      sum(!is.na(features[[paste0("shm_frequency_", region)]]))
    }, integer(1), USE.NAMES = FALSE),
    stringsAsFactors = FALSE
  )
}
