# Estimate BASELINe selection with SHazaM 1.3.2 in the pinned Immcantation
# image: collapseClones, then calcBaseline (focused test, HH_S5F targeting
# model), groupBaseline and summarizeBaseline.
# Each library/isotype/clone contributes one collapsed sequence. The module also
# masks uncovered positions and incomplete codons before estimation.
# Methods: Yaari et al. 2012 (doi:10.1093/nar/gks457), with the S5F model from
# Yaari et al. 2013 (doi:10.3389/fimmu.2013.00358).

#' Mask observed alignment gaps that cover a germline nucleotide.
#'
#' IMGT-aligned AIRR strings use gaps both as shared numbering spacers and as
#' padding where the observed sequence does not cover a germline position. A
#' shared gap must remain aligned, but an observed `.` or `-` opposite germline
#' A/C/G/T is not a callable point-mutation opportunity. Convert only those
#' positions to N. Inputs may be equal-length vectors or one may have length one;
#' each observed/germline string pair must itself have equal aligned length.
mask_unobserved_alignment_gaps <- function(sequence, germline) {
  n <- max(length(sequence), length(germline))
  if (!n) return(character(0))
  if (!length(sequence) %in% c(1L, n) || !length(germline) %in% c(1L, n)) {
    stop("sequence and germline must have equal vector lengths or length one")
  }
  sequence <- rep_len(as.character(sequence), n)
  germline <- rep_len(as.character(germline), n)

  vapply(seq_len(n), function(i) {
    if (is_airr_null(sequence[i]) || is_airr_null(germline[i])) {
      stop("sequence and germline alignments must be non-missing")
    }
    sequence_chars <- strsplit(toupper(sequence[i]), "", fixed = TRUE)[[1]]
    germline_chars <- strsplit(toupper(germline[i]), "", fixed = TRUE)[[1]]
    if (length(sequence_chars) != length(germline_chars)) {
      stop("sequence and germline alignments must have equal string lengths")
    }
    unobserved <- sequence_chars %in% c(".", "-") &
      germline_chars %in% c("A", "C", "G", "T")
    sequence_chars[unobserved] <- "N"
    paste(sequence_chars, collapse = "")
  }, character(1), USE.NAMES = FALSE)
}

#' Mask incomplete IMGT codons before selection counts and opportunities.
#'
#' A codon contributes only when both sequences contain three A/C/G/T bases.
#' Shared numbering gaps retain their positions. Codon boundaries follow the
#' aligned frame from position one, including gaps, and a terminal partial codon
#' is masked. Germline strings stay unchanged for the targeting model.
mask_incomplete_selection_codons <- function(sequence, germline) {
  sequence <- mask_unobserved_alignment_gaps(sequence, germline)
  germline <- rep_len(as.character(germline), length(sequence))
  vapply(seq_along(sequence), function(i) {
    observed <- strsplit(sequence[i], "", fixed = TRUE)[[1]]
    reference <- strsplit(toupper(germline[i]), "", fixed = TRUE)[[1]]
    for (start in seq.int(1L, length(observed), by = 3L)) {
      positions <- seq.int(start, min(start + 2L, length(observed)))
      complete <- length(positions) == 3L &&
        all(observed[positions] %in% c("A", "C", "G", "T")) &&
        all(reference[positions] %in% c("A", "C", "G", "T"))
      if (!complete) {
        shared_gap <- observed[positions] %in% c(".", "-") &
          reference[positions] %in% c(".", "-")
        observed[positions[!shared_gap]] <- "N"
      }
    }
    paste(observed, collapse = "")
  }, character(1), USE.NAMES = FALSE)
}

#
# Antigen-selection pressure on a library's B-cell repertoire, by BASELINe.
#
# BASELINe asks whether a V region carries more or fewer replacement mutations
# than the underlying mutation process would produce on its own, separately in the
# complementarity determining regions and in the frameworks around them. The
# answer is a Sigma per region: positive is selection for replacement, negative is
# selection against it. Each sequence contributes a probability density rather
# than a point estimate, and the densities of a group are convolved into one.
#
# Two decisions shape what a "sequence" is here. Clones are inferred from all of
# a donor's in-scope libraries, so the unit entering the convolution is one
# effective clone sequence per (library, isotype group, clone): a clonal expansion
# counts once, not once per member, and a clone seen in two libraries or in two
# isotype groups speaks once in each. Members are not weighted by read support,
# because each member casts one vote in the consensus. Membership resolves BCR
# backbone variants and then merges identical retained DNA within each library.
# Each distinct retained full sequence therefore casts one vote.
#
# The several CDR3s TRUST4 can report for one contig come from its AIRR row
# construction, in trust-airr.pl at v1.1.10 lines 398 to 455:
# https://github.com/liulab-dfci/TRUST4/blob/v1.1.10/trust-airr.pl#L398-L455

# The two halves of the IMGT V region a focused BASELINe test compares, in the
# order they are reported. CDR3 and FWR4 lie outside IMGT_V and are not reported.
SELECTION_REGIONS <- c("cdr", "fwr")

#' The released columns describing one region's estimate and evidence.
#'
#' `region` is one of SELECTION_REGIONS. Returns the column names in released
#' order: the effective clone sequences with a callable position in the region,
#' the usable-density count, availability status, then Sigma and its 95%
#' interval.
selection_region_columns <- function(region) {
  c(paste0("n_clone_consensus_sequences_with_positions_", region),
    paste0("n_clone_consensus_sequences_used_", region),
    paste0("sigma_", region, "_status"),
    paste0("sigma_", region),
    paste0("sigma_", region, "_ci_low"),
    paste0("sigma_", region, "_ci_high"))
}

#' Count the effective clone sequences with a callable position in one region.
#'
#' A position is callable when the effective sequence and its reconstructed
#' germline both hold A, C, G or T there after selection masking, so an
#' uncovered or incomplete-codon position is not. Each effective sequence counts
#' once, whether or not it supplies a usable density. Returns a one-row frame
#' holding `n_clone_consensus_sequences_with_positions_<region>`.
summarise_selection_context <- function(representatives, region,
                                        region_definition = selection_region_definition()) {
  n_with_positions <- 0L
  if (!is.null(representatives) && nrow(representatives)) {
    positions <- which(as.character(region_definition@boundaries) == region)
    for (i in seq_len(nrow(representatives))) {
      observed <- strsplit(representatives$clonal_sequence[i], "", fixed = TRUE)[[1]]
      germline <- strsplit(representatives$clonal_germline[i], "", fixed = TRUE)[[1]]
      callable <- germline[positions] %in% c("A", "C", "G", "T") &
        observed[positions] %in% c("A", "C", "G", "T")
      if (any(callable)) n_with_positions <- n_with_positions + 1L
    }
  }
  row <- data.frame(n_with_positions)
  names(row) <- paste0("n_clone_consensus_sequences_with_positions_", region)
  row
}

# The released library x isotype schema: the identifiers, how many effective
# clone sequences the library contributed to the group, then each region's
# evidence counts, status and estimate. Built from selection_region_columns() so
# the schema and the frame the estimator fills cannot spell a column differently.
SELECTION_COLUMNS <- c(
  "scope_id", "library_id", "donor_id", "tissue_id", "isotype_group", "n_clone_consensus_sequences",
  unlist(lapply(SELECTION_REGIONS, selection_region_columns), use.names = FALSE)
)

#' The region partition both the collapse and the estimate are made over.
#'
#' Returns `shazam::IMGT_V`, the V region split into its CDR and FWR halves.
#' A focused test needs exactly those two regions, and calcBaseline() otherwise
#' falls back silently to a local test, so the guard is asserted here rather than
#' left to the reader of either call site.
selection_region_definition <- function() {
  region_definition <- shazam::IMGT_V
  stopifnot(setequal(region_definition@regions, c("cdr", "fwr")))
  region_definition
}

# --------------------------------------------------------------- eligibility ---

#' Which records selection can say anything about (vectorised).
#'
#' `clone_id`, `sequence_alignment` and `germline_alignment_d_mask` are the three
#' columns a density needs. Returns one logical per record.
#'
#' Records reaching this stage have already passed the pipeline's common
#' post-IgBLAST gate: alignment N fraction, locus-concordant V and J calls, a
#' locus-concordant D call where one is present, productivity, a junction length
#' that is a multiple of three and positive read support. Selection adds
#' only what is intrinsic to the analysis (a scoped clone to be collapsed
#' within, and a D-masked germline to compare against) and no second length or
#' codon cutoff, because SHazaM already records zero usable densities where an
#' estimate cannot be made.
selection_record_eligible <- function(clone_id, sequence_alignment,
                                      germline_alignment_d_mask) {
  !is_airr_null(clone_id) &
    !is_airr_null(sequence_alignment) &
    !is_airr_null(germline_alignment_d_mask)
}

#' Stratify the eligible records into the groups Sigma is reported for.
#'
#' `records` are the eligible records, carrying the `isotype_group` their C call
#' resolved to, or NA where it did not resolve. Returns the same columns with
#' `isotype_group` rewritten, and with the records that resolved appearing twice.
#'
#' Every record enters `all_IGH`, which is the C-independent library endpoint, and
#' a record whose call resolved to IgM, IgG or IgA also enters that group. This is
#' the study's stratification: the isotype rows are not a partition of the all-IGH
#' row but separate readings of the same repertoire, so a clone can yield one
#' effective sequence in all-IGH and one in each isotype it contains, and none of
#' them is weighted by how many members it had. A record whose C call is absent or
#' ambiguous still counts in all-IGH, where the C call is not used.
#'
#' The stop() below is a precondition, not a policy on empty repertoires: a donor
#' with no eligible record is reported through build_empty_selection_rows()
#' instead, and never reaches the stratification. Reaching it means a caller
#' skipped that branch, which is worth saying here rather than letting an empty
#' frame travel on into the collapse and the convolution.
build_isotype_grouped_selection_input <- function(records) {
  if (!nrow(records)) {
    stop("no eligible IGH record to stratify: this donor has no clone-assigned ",
         "record with a D-masked germline in this scope")
  }
  all_eligible <- records
  all_eligible$isotype_group <- "all_IGH"
  resolved_eligible <- records[!is.na(records$isotype_group), , drop = FALSE]
  rbind(all_eligible, resolved_eligible)
}

# ------------------------------------------------------------- the estimate ---

#' Require a finite consensus frequency in (0, 1].
#'
#' At zero, SHazaM's thresholdedFreq method admits unobserved characters and
#' its deterministic ordering can replace every position with A. Validate before
#' either collapse or the donor entry point's empty-input exit.
validate_selection_collapse_frequency <- function(collapse_min_freq) {
  if (!is.numeric(collapse_min_freq) || is.complex(collapse_min_freq) ||
      length(collapse_min_freq) != 1L || !is.finite(collapse_min_freq) ||
      collapse_min_freq <= 0 || collapse_min_freq > 1) {
    stop("collapse_min_freq must be one finite number in (0, 1]")
  }
  invisible(collapse_min_freq)
}

#' Reduce each clone to one effective sequence.
#'
#' `contigs` are the grouped records, `collapse_min_freq` the frequency a base must
#' reach to be called, in (0, 1], `threads` the cores handed to shazam, and
#' `region_definition` the partition. Returns `effective_sequences`, one row per
#' (library, isotype group, clone) with its `clonal_sequence` and
#' `clonal_germline`, and `n_masked`, how many changed after masking uncovered
#' positions or incomplete codons.
#'
#' `fields` is the package-supported way to keep a donor-pooled clone separate by
#' library and isotype group, so one clone yields one effective sequence in each
#' rather than one for the whole donor. Ties are not broken stochastically and
#' ambiguous positions are not included, so the same input gives the same
#' consensus.
#'
#' The masking afterwards is not cosmetic. SHazaM treats a dot in the observed
#' sequence as a permitted IMGT gap when calculating expected mutation opportunity,
#' but IgBLAST also uses gaps to pad sequence that was not observed, which a
#' truncated bulk RNA-seq assembly produces routinely. Left alone, an unobserved
#' CDR would be counted as germline that failed to mutate, and the region would
#' report negative selection that is really an absence of data. Only the second
#' case is masked: an observed gap opposite an A/C/G/T germline base. Gaps shared
#' by sequence and germline are structural IMGT spacers and must stay aligned.
#' Incomplete codons are then masked in the aligned frame, so observed mutations
#' and expected opportunities use the same complete-codon evidence.
collapse_clones_to_effective_sequences <- function(
    contigs, collapse_min_freq, threads,
    region_definition = selection_region_definition()) {
  validate_selection_collapse_frequency(collapse_min_freq)
  keys <- c("library_id", "isotype_group", "clone_id")
  collapse_input <- contigs |>
    dplyr::select(
      dplyr::all_of(keys), sequence_alignment, germline_alignment_d_mask
    ) |>
    dplyr::arrange(
      library_id, isotype_group, clone_id,
      sequence_alignment, germline_alignment_d_mask
    )

  collapsed <- shazam::collapseClones(
    as.data.frame(collapse_input),
    cloneColumn = "clone_id",
    sequenceColumn = "sequence_alignment",
    germlineColumn = "germline_alignment_d_mask",
    regionDefinition = region_definition,
    method = "thresholdedFreq",
    minimumFrequency = collapse_min_freq,
    includeAmbiguous = FALSE,
    breakTiesStochastic = FALSE,
    fields = c("library_id", "isotype_group"),
    nproc = threads
  )
  collapsed <- collapsed |>
    dplyr::select(
      library_id, isotype_group, clone_id, clonal_sequence, clonal_germline
    )

  masked <- mask_incomplete_selection_codons(
    collapsed$clonal_sequence, collapsed$clonal_germline
  )
  n_masked <- sum(masked != collapsed$clonal_sequence)
  collapsed$clonal_sequence <- masked

  list(effective_sequences = collapsed, n_masked = n_masked)
}

#' Convolve one isotype group's effective sequences into a Sigma per library.
#'
#' `representatives` are the effective sequences of every group, `isotype_group`
#' the one being estimated, `threads` the cores, and `region_definition` the
#' partition. Returns NULL when the group holds nothing, otherwise the
#' `summary` rows summarizeBaseline() produced and the `pdf_counts` behind them.
#'
#' calcBaseline() derives the observed and expected mutation counts itself, since
#' this minimal frame supplies none. The rows are ordered by their sequence
#' content rather than by the clone label because the convolution is order
#' sensitive at machine precision, and a clone label carries its scope namespace.
#' Ordering on content keeps the input order independent of that namespace.
#' groupBaseline()@numbOfSeqs is the number of densities actually retained per
#' region, which is not the number of clones. The focused test counts, for one
#' region, its replacement and silent mutations plus the silent mutations of the
#' other region. A sequence contributes no density to a region when that count is
#' zero, or when the region has no expected replacement opportunity, as in a
#' region it did not observe. A sequence whose only mutation is silent and in FWR
#' therefore still contributes a CDR density.
estimate_baseline_selection <- function(representatives, isotype_group, threads,
                                        region_definition = selection_region_definition()) {
  if (is.null(representatives) || !nrow(representatives)) return(NULL)
  representatives <- representatives[
    representatives$isotype_group == isotype_group, , drop = FALSE
  ]
  if (!nrow(representatives)) return(NULL)
  baseline_input <- representatives[, c(
    "library_id", "clone_id",
    "clonal_sequence", "clonal_germline"
  ), drop = FALSE]
  baseline_input <- baseline_input |>
    dplyr::arrange(
      library_id, clonal_sequence, clonal_germline, clone_id
    )

  baseline <- shazam::calcBaseline(
    baseline_input,
    sequenceColumn = "clonal_sequence",
    germlineColumn = "clonal_germline",
    testStatistic = "focused",
    regionDefinition = region_definition,
    targetingModel = shazam::HH_S5F,
    mutationDefinition = NULL,
    nproc = threads
  )

  grouped <- shazam::groupBaseline(
    baseline, groupBy = "library_id", nproc = threads
  )
  summary <- as.data.frame(
    shazam::summarizeBaseline(grouped, returnType = "df", nproc = threads)
  )
  require_columns(summary, c("library_id", "region", "baseline_sigma",
                             "baseline_ci_lower", "baseline_ci_upper"),
                  "summarizeBaseline output")

  n_pdfs <- grouped@numbOfSeqs
  pdf_counts <- data.frame(
    library_id = as.character(grouped@db$library_id),
    n_pdf_cdr = as.integer(n_pdfs[, "cdr"]),
    n_pdf_fwr = as.integer(n_pdfs[, "fwr"]),
    stringsAsFactors = FALSE
  )

  list(summary = summary, pdf_counts = pdf_counts)
}

# ----------------------------------------------------------- the output rows ---

#' Read one library's estimate for one region off the group's results.
#'
#' `baseline_result` is what estimate_baseline_selection() returned for the isotype
#' group, or NULL when the group held nothing; `library_id` names the library;
#' `region` is one of SELECTION_REGIONS; `n_clone_consensus_sequences` is how
#' many effective sequences the library contributed to the group. Returns the
#' usable-density count, availability status, Sigma and interval limits for that
#' region, named as selection_region_columns() names them.
#'
#' A non-finite Sigma counts as no estimate rather than as a value, and a Sigma
#' with no retained density behind it is not an estimate either.
#'
#' The three statuses are three different findings and a reader has to tell them
#' apart. "estimated" means only that a finite Sigma with usable densities
#' exists; it is not a quality judgement. "no_eligible_clones" means the library contributed no effective clone
#' sequence to this group: nothing of that isotype was recovered, or what was
#' recovered carried no clone assignment or no D-masked germline. Non-recovery
#' from bulk RNA-seq does not establish that cells of that isotype were absent.
#' "no_usable_mutation_pdfs" means the library did contribute, but the region has
#' no finite estimate: no effective clone sequence yielded a density for the
#' region (see estimate_baseline_selection()), or SHazaM returned a non-finite
#' Sigma.
#'
#' The used count is NA rather than zero for a library that contributed clones
#' yet is absent from the estimate: that is a case nothing was reported about, and
#' zero would read as a measured answer. A library with no clones really did use
#' zero densities, which is a known count.
summarise_region_for_library <- function(baseline_result, library_id, region,
                                         n_clone_consensus_sequences) {
  estimate <- NULL
  n_used <- if (n_clone_consensus_sequences) NA_integer_ else 0L
  if (!is.null(baseline_result)) {
    estimate <- baseline_result$summary[
      as.character(baseline_result$summary$library_id) == library_id &
        as.character(baseline_result$summary$region) == region, , drop = FALSE]
    counts <- baseline_result$pdf_counts[
      as.character(baseline_result$pdf_counts$library_id) == library_id, , drop = FALSE]
    if (nrow(counts)) n_used <- as.integer(counts[[paste0("n_pdf_", region)]][1])
  }

  has_density <- is.finite(n_used) && n_used > 0L
  has_sigma <- !is.null(estimate) && nrow(estimate) &&
    is.finite(estimate$baseline_sigma[1])
  status <- if (has_sigma && has_density) {
    "estimated"
  } else if (!n_clone_consensus_sequences) {
    "no_eligible_clones"
  } else {
    "no_usable_mutation_pdfs"
  }

  row <- data.frame(
    n_used,
    status,
    if (status == "estimated") estimate$baseline_sigma[1] else NA_real_,
    if (status == "estimated") estimate$baseline_ci_lower[1] else NA_real_,
    if (status == "estimated") estimate$baseline_ci_upper[1] else NA_real_,
    stringsAsFactors = FALSE
  )
  names(row) <- c(paste0("n_clone_consensus_sequences_used_", region), paste0("sigma_", region, "_status"),
                  paste0("sigma_", region), paste0("sigma_", region, "_ci_low"),
                  paste0("sigma_", region, "_ci_high"))
  row
}

#' Build the donor's released rows: one per library, per isotype group.
#'
#' `baseline_results` is the estimate for each group, named by group;
#' `representatives` are the effective sequences, or NULL; `samples` are the
#' donor's libraries from the scope manifest, with their donor and tissue;
#' `scope_id` stamps every row; `isotype_groups` fixes the groups and their order.
#' Returns the frame carrying SELECTION_COLUMNS.
#'
#' Every library gets a row in every group whether or not anything was estimated
#' for it. An absent row and a row saying nothing was estimated read differently
#' downstream, and only the row states what this stage actually found.
build_selection_result_rows <- function(baseline_results, representatives, samples,
                                        scope_id,
                                        isotype_groups = ANALYSIS_IGH_ISOTYPE_GROUPS) {
  rows <- list()
  for (i in seq_len(nrow(samples))) {
    library_id <- as.character(samples$library_id[i])
    for (group in isotype_groups) {
      members <- if (is.null(representatives)) NULL else representatives[
        representatives$library_id == library_id &
          representatives$isotype_group == group, , drop = FALSE]
      n_clone_consensus_sequences <- if (is.null(members)) 0L else nrow(members)
      row <- data.frame(
        scope_id = scope_id,
        library_id = library_id,
        donor_id = as.character(samples$donor_id[i]),
        tissue_id = as.character(samples$tissue_id[i]),
        isotype_group = group,
        n_clone_consensus_sequences = as.integer(n_clone_consensus_sequences),
        stringsAsFactors = FALSE
      )
      for (region in SELECTION_REGIONS) {
        row <- cbind(row, summarise_region_for_library(
          baseline_results[[group]], library_id, region, n_clone_consensus_sequences),
          summarise_selection_context(members, region))
      }
      rows[[length(rows) + 1L]] <- row
    }
  }
  do.call(rbind, rows)[, SELECTION_COLUMNS, drop = FALSE]
}

#' Build the donor's released rows when there is nothing to estimate.
#'
#' `samples` are the donor's libraries from the scope manifest; `scope_id` stamps
#' every row; `isotype_groups` fixes the groups and their order. Returns the frame
#' carrying SELECTION_COLUMNS, one row per library and group, each reporting no
#' clone and no Sigma.
#'
#' A donor whose pooled records hold no clone-assigned sequence with a D-masked
#' germline has an empty repertoire as far as selection is concerned, which is a
#' valid reading rather than an error. It is reported through the same row builder
#' the estimate uses, with no group result to read from, so "nothing was
#' estimated" has one definition here regardless of why nothing was.
build_empty_selection_rows <- function(samples, scope_id,
                                       isotype_groups = ANALYSIS_IGH_ISOTYPE_GROUPS) {
  no_results <- stats::setNames(vector("list", length(isotype_groups)),
                                isotype_groups)
  build_selection_result_rows(no_results, NULL, samples, scope_id, isotype_groups)
}
