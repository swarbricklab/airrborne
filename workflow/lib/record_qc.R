# Check which annotated records satisfy every common QC requirement.
# filter_airr.R applies these checks to one library's IgBLAST/MakeDb rows and
# writes their results to the QC audit. A failure removes the whole record.
# Grouping, C-call qualification and regional SHM apply further requirements
# to their own measurements. Base R only.

#' Check the shared requirements for keeping a record after IgBLAST annotation.
#'
#' Retain productive, locus-consistent records whose junction length is
#' divisible by three, with positive TRUST4 read support and an aligned N
#' fraction at or below `max_alignment_n_fraction`
#' (`params.airr_qc.max_alignment_n_fraction`). The N fraction counts N and n in
#' `sequence_alignment` over its positions other than the gap characters `.` and
#' `-`, and an alignment with no such position fails at every setting. The AIRR
#' `vj_in_frame` flag is not tested. Every called V, J or D gene must agree with
#' the assigned locus. A missing D call is allowed; a D call on a locus without
#' D genes is not. A V/J call naming several genes can pass when every candidate
#' agrees on locus. Later grouping checks whether the gene assignment resolves
#' through gene_calls.R; common QC does not require one assigned V or J gene.
#'
#' Constant-region qualification is separate because junction and V analyses
#' do not require a callable constant region. Isotype composition, SHM and
#' selection read the separate table of qualified constant calls written by
#' qualify_isotype_calls.R when their measurement needs it.
#'
#' These checks impose no minimum length for the whole alignment. Truncated
#' bulk RNA assemblies can retain a usable junction or measured subregion.
#' SHM and isotype apply coverage requirements to the regions they inspect;
#' junction analyses require their own grouping eligibility.
#'
#' Junction length is checked separately from productivity because the length
#' comes from the aligned junction and productivity comes from IgBLAST's summary.
#' A record must pass every check to remain in the common-QC table.
#'
#' Returns one row per record with six logical columns, TRUE where the record
#' passes: `locus_match`, `d_locus_match`, `low_n`, `productive`,
#' `junction_mod3` and `positive_support`. The QC audit names each criterion row
#' `qc_<column>`. Every criterion is evaluated on the same input, so a record can
#' fail several at once.
airr_qc_flags <- function(df, max_alignment_n_fraction = 0.1) {
  required <- c("locus", "v_call", "j_call", "d_call", "sequence_alignment",
                "productive", "junction_length", "consensus_count")
  missing <- setdiff(required, names(df))
  if (length(missing)) {
    stop("post-IgBLAST QC requires: ", paste(missing, collapse = ", "))
  }
  if (length(max_alignment_n_fraction) != 1L ||
      !is.finite(max_alignment_n_fraction) ||
      max_alignment_n_fraction < 0 || max_alignment_n_fraction > 1) {
    stop("max_alignment_n_fraction must be between 0 and 1")
  }

  alignment <- as.character(df$sequence_alignment)
  aligned <- nchar(gsub("[-.]", "", alignment))
  n_bases <- nchar(gsub("[^Nn]", "", alignment))
  n_fraction <- ifelse(!is.na(aligned) & aligned > 0L, n_bases / aligned, Inf)

  junction_length <- suppressWarnings(as.numeric(df$junction_length))
  support <- as_read_support(df$consensus_count, "consensus_count")
  locus <- as.character(df$locus)

  # Missing D calls pass because IgBLAST can leave D unresolved. Present calls
  # must match the assigned locus, including its ability to carry a D gene.
  # gene_call_locus_compatible() owns the allowed loci and rejects TRA/TRD
  # disagreements rather than assigning the record to either locus here.
  d_absent <- is_airr_null(df$d_call)
  d_locus_match <- rep(TRUE, length(d_absent))
  if (any(!d_absent)) {
    d_locus_match[!d_absent] <-
      gene_call_locus_compatible(df$d_call[!d_absent], locus[!d_absent], "D")
  }

  data.frame(
    locus_match = gene_call_locus_compatible(df$v_call, locus, "V") &
      gene_call_locus_compatible(df$j_call, locus, "J"),
    d_locus_match = d_locus_match,
    low_n = is.finite(n_fraction) & n_fraction <= max_alignment_n_fraction,
    productive = parse_airr_boolean(df$productive) == "true",
    junction_mod3 = is.finite(junction_length) & junction_length >= 0 &
      junction_length %% 3 == 0,
    positive_support = support > 0,
    stringsAsFactors = FALSE
  )
}

#' Count failed common-QC criteria per record.
#'
#' flags is the logical, non-missing frame returned by airr_qc_flags(). Returns
#' one integer per record, with zero denoting a pass. Stops on non-logical or
#' missing flags so an unknown result cannot introduce missing rows when filtering.
airr_qc_fail_counts <- function(flags) {
  not_logical <- names(flags)[!vapply(flags, is.logical, logical(1))]
  if (length(not_logical)) {
    stop("QC flags must be logical; these are not: ",
         paste(not_logical, collapse = ", "))
  }
  has_na <- names(flags)[vapply(flags, anyNA, logical(1))]
  if (length(has_na)) {
    stop("QC flags must not be NA; these are: ", paste(has_na, collapse = ", "))
  }
  as.integer(rowSums(!as.matrix(flags)))
}

#' TRUE where a record satisfies every common QC criterion.
#'
#' `flags` is airr_qc_flags()' output. Returns one logical per record.
airr_qc_pass <- function(flags) airr_qc_fail_counts(flags) == 0L
