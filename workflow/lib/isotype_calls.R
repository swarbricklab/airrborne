# Decide whether an IGH constant-region gene call has enough evidence for isotype
# analysis. The constant region identifies the antibody class or subclass.
# c_call_qc_flags() checks the evidence; mask_low_evidence_c_calls() blanks failed
# calls. All records stay in the table, with their original c_call preserved.
# isotype_resolution.R names surviving calls; isotype_metrics.R counts support.
# These helpers use base R. Common record QC is in record_qc.R.

#' Flag constant calls with sufficient aligned evidence, without removing rows.
#'
#' Record fields must have equal lengths or length one. Three empty call, locus
#' and position vectors return zero rows; completeness may then be empty or scalar.
#'
#' MakeDb supplies the C call, which can list several equally good gene matches.
#' reannotated_airr.R reads the first C alignment in IgBLAST's best-first output:
#' c_informative_positions counts its aligned A, C, G and T query bases.
#' A missing call has nothing to mask and passes
#' `locus_match`. A present call is masked when it is cross-locus or
#' wrong-segment, its candidate list is incomplete, or its count is missing or below
#' `min_informative_c_alignment_positions`
#' (`params.airr_qc.min_informative_c_alignment_positions`). Zero removes the
#' length floor, and a finite count is still required.
c_call_qc_flags <- function(c_call, locus, c_informative_positions,
                            min_informative_c_alignment_positions, candidates_complete = TRUE) {
  if (length(min_informative_c_alignment_positions) != 1L ||
      !is.finite(min_informative_c_alignment_positions) ||
      min_informative_c_alignment_positions < 0) {
    stop("min_informative_c_alignment_positions must be one non-negative number")
  }
  n <- max(length(c_call), length(locus), length(c_informative_positions))
  lengths <- c(length(c_call), length(locus), length(c_informative_positions),
               length(candidates_complete))
  if (any(!lengths %in% c(1L, n))) {
    stop("C-call fields must have equal lengths or length one")
  }
  c_call <- rep_len(c_call, n)
  locus <- rep_len(locus, n)
  positions <- suppressWarnings(as.numeric(rep_len(c_informative_positions, n)))
  present <- !is_airr_null(c_call)
  locus_match <- rep(TRUE, n)
  if (any(present)) {
    locus_match[present] <- gene_call_locus_compatible(
      c_call[present], locus[present], "C"
    )
  }
  enough_evidence <- is.finite(positions) &
    positions >= min_informative_c_alignment_positions
  candidates_complete <- rep_len(candidates_complete, n)
  data.frame(
    candidates_complete = candidates_complete,
    present = present,
    locus_match = locus_match,
    enough_evidence = enough_evidence,
    callable = present & locus_match & enough_evidence & candidates_complete,
    stringsAsFactors = FALSE
  )
}

#' Blank constant-region calls that failed the evidence checks.
#'
#' `c_call` is the record's constant-region call and `callable` the decision
#' c_call_qc_flags() already made about it. Returns `c_call` with every
#' non-callable entry set to NA. A surviving call carries enough IGH
#' constant-region evidence, but it still needs a name from resolve_isotype() in
#' isotype_resolution.R. A tie across classes or the IGHGP pseudogene can pass
#' these evidence checks and still have no resolved isotype.
mask_low_evidence_c_calls <- function(c_call, callable) {
  c_call[!callable] <- NA_character_
  c_call
}

#' Join qualified calls to one library's common-QC IGH records by sequence_id.
#'
#' Both tables must have unique keys and describe exactly the same records.
#' The caller checks library_id. Return df in input order with qualified_c_call
#' attached and the original c_call unchanged. Empty tables return an empty column.
join_qualified_isotype_calls <- function(df, calls) {
  require_columns(calls, c("sequence_id", "qualified_c_call"), "Qualified isotype calls")
  airr_ids <- as.character(df$sequence_id)
  call_ids <- as.character(calls$sequence_id)
  if (anyDuplicated(airr_ids)) {
    stop("Common-QC IGH AIRR has duplicate sequence_id keys")
  }
  if (anyDuplicated(call_ids)) {
    stop("Qualified isotype calls have duplicate sequence_id keys")
  }
  if (!setequal(airr_ids, call_ids)) {
    stop(sprintf(
      "Qualified isotype-call keys do not exactly match common-QC IGH AIRR (missing=%d, extra=%d)",
      length(setdiff(airr_ids, call_ids)),
      length(setdiff(call_ids, airr_ids))
    ))
  }
  df$qualified_c_call <- as.character(calls$qualified_c_call[match(airr_ids, call_ids)])
  df
}

# The shared per-library qualification record. Calls that fail keep their
# original spelling beside the missing qualified call and its evidence flags.
# qualify_isotype_calls.R appends isotype_class and isotype_subclass to the
# sidecar it writes.
ISOTYPE_CALL_COLUMNS <- c(
  "library_id", "sequence_id", "c_call", "qualified_c_call", "c_call_present",
  "c_call_locus_match", "c_call_enough_evidence", "c_call_candidates_complete", "c_call_callable"
)

#' Qualify constant-region calls once for every common-QC IGH record.
#'
#' records is the common-QC all-locus table for one library. Retain all its IGH
#' rows, including missing or failed calls. Qualification masks the call only;
#' it does not exclude a record from other analyses.
build_qualified_isotype_calls <- function(records, sample_label,
                                          min_c_positions) {
  require_columns(records, c("library_id", "sequence_id", "locus", "c_call",
                             "c_informative_positions"), "Common-QC AIRR")
  igh <- records[as.character(records$locus) == "IGH", , drop = FALSE]
  if (nrow(igh) && any(as.character(igh$library_id) != sample_label)) {
    stop("Common-QC AIRR contains an IGH library_id other than ", sample_label)
  }
  if (any(is_airr_null(igh$sequence_id)) || anyDuplicated(igh$sequence_id)) {
    stop("Common-QC IGH AIRR must have unique, non-missing sequence_id keys")
  }
  flags <- c_call_qc_flags(igh$c_call, igh$locus,
                           igh$c_informative_positions, min_c_positions,
                           if ("c_call_candidates_complete" %in% names(igh)) {
                             parse_airr_boolean(igh$c_call_candidates_complete) == "true"
                           } else rep(TRUE, nrow(igh)))
  calls <- data.frame(
    library_id = as.character(igh$library_id),
    sequence_id = as.character(igh$sequence_id),
    c_call = as.character(igh$c_call),
    qualified_c_call = mask_low_evidence_c_calls(as.character(igh$c_call), flags$callable),
    c_call_present = flags$present,
    c_call_locus_match = flags$locus_match,
    c_call_enough_evidence = flags$enough_evidence,
    c_call_candidates_complete = flags$candidates_complete,
    c_call_callable = flags$callable,
    stringsAsFactors = FALSE
  )
  calls <- calls[order(calls$sequence_id), ISOTYPE_CALL_COLUMNS, drop = FALSE]
  rownames(calls) <- NULL
  calls
}

