# Measure mutation counts, opportunities and support in the IMGT V partition
# with SHazaM 1.3.2 calcObservedMutations over IMGT_V_BY_REGIONS (Gupta et al.
# 2015, doi:10.1093/bioinformatics/btv359) in the pinned image. Germline
# occupancy defines regional coverage; shared numbering gaps provide no
# opportunity. The five region counts sum to the whole-V numerator and
# denominator. calculate_igh_features.R and the donor IGH summaries load these
# same definitions.

V_REGION_LEVELS <- c("fwr1", "cdr1", "fwr2", "cdr2", "fwr3")

# The characters SHazaM counts as an observed base when it builds `nonN`: the four
# unambiguous bases plus the ten IUPAC ambiguity codes. `N`, `.` and `-` are all
# absent, which is why an unobserved position is not a mutation opportunity.
SHAZAM_OBSERVED_BASE_CODES <- c("A", "C", "G", "T", "M", "R", "W", "S", "Y", "K",
                                "V", "H", "D", "B")



# --------------------------------------------------------------- partition ---

#' Sub-region label for every nucleotide position of the IMGT V region.
#'
#' Returns a character vector as long as the IMGT V window, holding one of
#' `V_REGION_LEVELS` per position. Raises if SHazaM's partition no longer covers
#' the window or no longer uses these five labels, so a package upgrade that
#' changed the numbering cannot pass silently.
imgt_v_region_of_position <- function() {
  boundaries <- as.character(shazam::IMGT_V_BY_REGIONS@boundaries)
  expected_length <- as.integer(shazam::IMGT_V@seqLength)
  if (length(boundaries) != expected_length) {
    stop("shazam::IMGT_V_BY_REGIONS covers ", length(boundaries),
         " position(s) but shazam::IMGT_V covers ", expected_length)
  }
  if (!setequal(unique(boundaries), V_REGION_LEVELS)) {
    stop("shazam::IMGT_V_BY_REGIONS labels its regions ",
         paste(sort(unique(boundaries)), collapse = ", "), "; expected ",
         paste(V_REGION_LEVELS, collapse = ", "))
  }
  boundaries
}

#' One IMGT-aligned sequence as an upper-case character vector of fixed length.
#'
#' Longer input is truncated to the window and shorter input is padded with `N`,
#' which is how SHazaM treats a sequence that stops inside the V region. A missing
#' or empty sequence gives all `N`, so it contributes no opportunity anywhere.
imgt_v_alignment_characters <- function(alignment, n_positions) {
  characters <- rep("N", n_positions)
  value <- toupper(as.character(alignment))
  if (length(value) != 1L) {
    stop("imgt_v_alignment_characters expects one sequence; got ", length(value))
  }
  if (is.na(value) || !nzchar(value)) return(characters)
  observed <- strsplit(substr(value, 1L, n_positions), "", fixed = TRUE)[[1]]
  characters[seq_along(observed)] <- observed
  characters
}

# ------------------------------------------------------------- denominators ---

#' Positions of each sub-region that the germline actually occupies.
#'
#' Returns one count per IMGT V sub-region. Germline-occupied positions define
#' the coverage denominator; numbering spacers are excluded.
count_germline_occupied_by_region <- function(germline_alignment,
                                              region_of_position =
                                                imgt_v_region_of_position()) {
  characters <- imgt_v_alignment_characters(germline_alignment,
                                            length(region_of_position))
  occupied <- characters %in% SHAZAM_OBSERVED_BASE_CODES
  counts <- vapply(V_REGION_LEVELS, function(region) {
    sum(occupied & region_of_position == region)
  }, integer(1))
  counts
}

#' Which IMGT V positions contributed a mutation opportunity, position by position.
#'
#' `sequence_alignment` and `germline_alignment` are one record's IMGT-gapped
#' strings, the germline being the D-masked one. Returns a logical vector as long
#' as the IMGT V window, `TRUE` where both strings carry an observed base.
#'
#' This is SHazaM's `nonN` test written out per position rather than counted, so
#' `sum()` of it equals the record's `informative_v`, and summing it within a
#' sub-region equals that sub-region's `informative_<region>`. It exists because
#' the per-sub-region vote summary below needs the positions themselves and not
#' just their count: it must average over exactly the positions the rate was
#' measured on.
imgt_v_informative_positions <- function(sequence_alignment, germline_alignment,
                                         n_positions =
                                           as.integer(shazam::IMGT_V@seqLength)) {
  observed <- imgt_v_alignment_characters(sequence_alignment, n_positions)
  germline <- imgt_v_alignment_characters(germline_alignment, n_positions)
  observed %in% SHAZAM_OBSERVED_BASE_CODES &
    germline %in% SHAZAM_OBSERVED_BASE_CODES
}

#' Read support behind each IMGT V sub-region of one record.
#'
#' `majority_votes` are the per-position votes agreeing with the assembled base
#' (`position_majority_votes()`), `informative` the mask
#' `imgt_v_informative_positions()` returned, and `region_of_position` the
#' partition. All three must describe the same window. Returns a list of three
#' named vectors over `V_REGION_LEVELS`: `median`, the median vote count; `n`, how
#' many positions that median was taken over; and `n_supported`, how many of those
#' positions carry at least the configured `min_position_votes`.
#'
#' A position counts only when it is both informative and carries a vote. The
#' vote statistic therefore uses a subset of the mutation-rate positions when
#' votes are missing. A missing vote is not a zero: it means no
#' contig coordinate mapped to that position, because of an IMGT numbering spacer,
#' a truncated end of the V region, a base TRUST4's substituted CDR3 rewrote, or a
#' record assembled in the opposite orientation to its contig. Treating those as
#' zero would manufacture a support deficit at the five prime end, which is where
#' bulk RNA-seq coverage is genuinely thinnest. `median` is `NA` where `n` is zero.
#'
#' `n_supported` is what the released SHM eligibility rule is written against. It is
#' a count rather than a fraction so that the denominator stays visible: compared
#' against `count_germline_occupied_by_region()` it says whether the whole
#' sub-region was both seen and read-supported, and compared against `n` it says how
#' much of what was seen rests on a single vote. `median` cannot answer either
#' question, because a median above the floor is compatible with individual
#' positions below it.
#'
#' The counts are votes, not read depth. TRUST4 writes a literal 1 on the
#' consensus base wherever a position's tally came out zero, so a median of 1 pools
#' one-read evidence with none. See `trust4_consensus_votes.R` for the assembler's
#' behaviour and the config's min_position_votes default for why the floor is two.
summarise_vote_support_by_region <- function(majority_votes, informative,
                                            min_position_votes,
                                            region_of_position =
                                              imgt_v_region_of_position()) {
  if (length(majority_votes) != length(region_of_position) ||
      length(informative) != length(region_of_position)) {
    stop("vote, informative and region vectors must all cover the ",
         length(region_of_position), "-position IMGT V window")
  }
  counted <- informative & !is.na(majority_votes)
  supported <- counted & majority_votes >= min_position_votes
  count_by_region <- function(mask) {
    vapply(V_REGION_LEVELS, function(region) {
      sum(mask & region_of_position == region)
    }, integer(1))
  }
  median_votes <- vapply(V_REGION_LEVELS, function(region) {
    values <- majority_votes[counted & region_of_position == region]
    if (!length(values)) return(NA_real_)
    as.numeric(stats::median(values))
  }, numeric(1))
  list(median = median_votes, n = count_by_region(counted),
       n_supported = count_by_region(supported))
}

# ------------------------------------------------------------ mutation load ---

#' Mutations and their opportunity within each IMGT V sub-region.
#'
#' Runs one `shazam::calcObservedMutations()` over `shazam::IMGT_V_BY_REGIONS` and
#' returns a list of two named integer vectors over `V_REGION_LEVELS`:
#' `mutations`, the replacement plus silent count SHazaM scored in that
#' sub-region, and `informative`, SHazaM's `nonN` count of positions carrying an
#' observed base in both the record and its germline. Dividing one by the other
#' gives the sub-region mutation frequency; summing both across sub-regions gives
#' the whole-V frequency, because the two region definitions describe the same
#' window.
#'
#' `sequence_alignment` and `germline_alignment` are single IMGT-gapped strings;
#' the germline is the D-masked one, matching the released whole-V rate.
#'
#' SHazaM scores a candidate mutation by substituting that one position into the
#' germline codon, and keeps it only if it is replacement or silent. A
#' substitution whose isolated codon would be a stop codon, or whose codon cannot
#' be translated, is therefore counted by neither. That is intrinsic to the
#' function and is deliberately kept: it is what makes these counts continuous
#' with the whole-V rate.
summarise_observed_mutations_by_region <- function(sequence_alignment,
                                                   germline_alignment) {
  raw <- shazam::calcObservedMutations(
    as.character(sequence_alignment),
    as.character(germline_alignment),
    regionDefinition = shazam::IMGT_V_BY_REGIONS,
    mutationDefinition = NULL,
    ambiguousMode = "eitherOr",
    returnRaw = TRUE,
    frequency = FALSE
  )

  # `nonN` is named per sub-region. Index by name rather than by position: SHazaM
  # returns it in alphabetical order, not in the 5' to 3' order used here.
  informative <- stats::setNames(integer(length(V_REGION_LEVELS)), V_REGION_LEVELS)
  named <- intersect(names(raw$nonN), V_REGION_LEVELS)
  informative[named] <- as.integer(round(as.numeric(raw$nonN[named])))

  # `pos` holds one row per kept mutation, with its replacement and silent
  # scores. Summing those two is how SHazaM itself forms a per-region count, so
  # counting them rather than counting rows keeps an ambiguous base weighted the
  # way SHazaM weights it. When nothing was kept, `pos` is a bare NA rather than
  # an empty frame.
  mutations <- stats::setNames(integer(length(V_REGION_LEVELS)), V_REGION_LEVELS)
  if (is.data.frame(raw$pos) && nrow(raw$pos) > 0L) {
    scored <- as.numeric(raw$pos$r) + as.numeric(raw$pos$s)
    region <- factor(as.character(raw$pos$region), levels = V_REGION_LEVELS)
    if (anyNA(region)) {
      stop("shazam returned a mutation outside the five IMGT V sub-regions")
    }
    totals <- tapply(scored, region, sum)
    totals[is.na(totals)] <- 0
    mutations[V_REGION_LEVELS] <- as.integer(round(totals[V_REGION_LEVELS]))
  }

  list(mutations = mutations, informative = informative)
}
