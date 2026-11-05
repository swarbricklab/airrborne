# Identify TRUST4 records whose CDR3 contains bases inferred by the assembler.
# CDR3 is the variable junction region later used to distinguish rearrangements.
# parse_trust4_airr.R removes flagged records before IgBLAST annotation. Other
# records retain their IDs and recovered read support, including alternative
# CDR3s reported on one assembly. These helpers use base R.

TRUST4_IMPUTED_CDR3_SCORE <- 0.01

# How close a score has to be to the sentinel to be it. The score is read back
# out of a decimal text field, so it is compared within a tolerance rather than
# for exact floating-point equality. The filter assumes that a CDR3 TRUST4
# actually observed does not score within this tolerance of the sentinel, so
# widening the tolerance would start removing observed CDR3s.
TRUST4_IMPUTED_CDR3_TOLERANCE <- 1e-6

#' Read the contig CDR3 annotation score repeated on each cdr3.out row.
#'
#' `path` is one library's cdr3.out. The file has no header; columns 1-2 are
#' `<consensus_id>` and `<index>`, whose joined spelling matches AIRR identifiers.
#' AIRR suffixes rank reported variants and need not retain the cdr3.out index.
#' In pinned TRUST4 1.1.10, column 10 repeats the same contig annotation score
#' for every variant, so this difference does not change imputation status.
#' Having no header is also why this cannot go through read_pipeline_table(),
#' whose first line is always the column names; the typing decisions are the
#' same.
#'
#' Returns one row per reported CDR3, with the assembly's `consensus_id`, the
#' joined `cdr3_key`, and the numeric `cdr3_score`. A duplicated key is an error,
#' because the score it carries would then be ambiguous for the sequence_id it
#' joins to.
#'
#' TRUST4 writes thirteen columns in cdr3.out, and all thirteen are required
#' here: a shorter row means the format changed and column 10 is no longer the
#' CDR3 score.
#'
#' TRUST4 cdr3.out column format (CDR3 score):
#'   https://github.com/liulab-dfci/TRUST4
read_trust4_cdr3_scores <- function(path) {
  # A library with no reported CDR3 has a zero-byte, headerless file. Missing
  # files and nonempty malformed rows still fail through the reader below.
  if (isTRUE(file.info(path)$size == 0)) {
    return(data.frame(consensus_id = character(), cdr3_key = character(),
                      cdr3_score = numeric(), stringsAsFactors = FALSE))
  }
  cdr3 <- read.delim(
    path, header = FALSE, stringsAsFactors = FALSE, check.names = FALSE,
    colClasses = "character", na.strings = character(0)
  )
  if (ncol(cdr3) < 13L) {
    stop("TRUST4 cdr3.out has fewer than 13 columns: ", path)
  }
  cdr3_key <- paste0(cdr3[[1L]], "_", cdr3[[2L]])
  if (anyDuplicated(cdr3_key)) {
    stop("TRUST4 cdr3.out has duplicate <consensus_id>_<index> keys")
  }
  data.frame(
    consensus_id = cdr3[[1L]],
    cdr3_key     = cdr3_key,
    cdr3_score   = suppressWarnings(as.numeric(cdr3[[10L]])),
    stringsAsFactors = FALSE
  )
}

#' Which records carry a CDR3 TRUST4 imputed rather than observed.
#'
#' `sequence_ids` are the AIRR sequence_ids to classify; `cdr3_keys` and
#' `cdr3_scores` are the two columns read_trust4_cdr3_scores() returned;
#' `sentinel` and `tolerance` are the imputed-score constants above.
#'
#' In the pinned TRUST4 release the imputation routine returns early unless the
#' V gene name starts with "T" (ImputeCDR3 in SeqSet.hpp at v1.1.10:
#' https://github.com/liulab-dfci/TRUST4/blob/v1.1.10/SeqSet.hpp#L5867-L5868),
#' so this identifies TRUST4's inferred TCR CDR3 completions. The test below
#' reads the score alone, so it is applied at
#' every locus and the caller reports the loci it fired on. Because
#' the identity-forming junction is re-derived downstream by IgBLAST from the
#' contig's own `sequence`, an imputed base would otherwise be treated as an
#' observed exact junction (the max_n gate counts only non-ACGT symbols).
#' Removing these rows before reannotation is what keeps imputed TCR junctions
#' out of every rearrangement group and overlap comparison.
#'
#' Returns a logical vector as long as `sequence_ids`, TRUE where the CDR3 was
#' imputed. Deciding only: the caller drops the flagged rows and accounts for
#' them. A sequence_id with no score, or with a non-numeric one, is an error,
#' because the rule cannot be applied to it and keeping it would let an imputed
#' junction through unexamined.
flag_imputed_cdr3_records <- function(sequence_ids, cdr3_keys, cdr3_scores,
                                      sentinel = TRUST4_IMPUTED_CDR3_SCORE,
                                      tolerance = TRUST4_IMPUTED_CDR3_TOLERANCE) {
  score_idx <- match(sequence_ids, cdr3_keys)
  if (any(is.na(score_idx))) {
    stop(sprintf(
      "TRUST4 cdr3.out is missing a CDR3 score for %d sequence_id(s), e.g. %s",
      sum(is.na(score_idx)),
      paste(utils::head(sequence_ids[is.na(score_idx)], 3L), collapse = ", ")))
  }
  matched_score <- cdr3_scores[score_idx]
  if (any(is.na(matched_score))) {
    stop("TRUST4 cdr3.out has a non-numeric CDR3 score for a matched sequence_id")
  }
  abs(matched_score - sentinel) < tolerance
}
