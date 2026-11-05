# Interpret V/J gene assignments without choosing an arbitrary tied gene.
# An AIRR call lists alleles that IgBLAST could not distinguish. A gene entity
# is either one gene or a reference-defined group of genes with identical complete
# sequences. vj_tie_kind() decides which applies; vj_canonical_token() names it.
# reannotated_airr.R writes these results as v_gene_status/j_gene_status and
# v_gene/j_gene. Grouping reads those columns rather than the first listed call.
# These helpers use base R.

.split_calls <- function(cell) {
  toks <- trimws(strsplit(as.character(cell), "[,;|]")[[1]])
  toks[!is_airr_null(toks) & toks != "." & toks != "*"]
}

#' Remove an allele suffix while preserving every genomic gene name.
#'
#' A duplicated-location D marker does not establish sequence identity.
canonical_gene_name <- function(x) {
  sub("\\*.*$", "", as.character(x))
}

# Return the distinct gene names in one call string, keeping their listed order.
.gene_set <- function(cell) {
  g <- canonical_gene_name(.split_calls(cell))
  unique(g[g != ""])
}

# Return the locus or loci compatible with one canonical gene token and segment.
# The TRA/TRD V pool supports either chain, resolved from consistent J, D and C
# evidence. Other V/J tokens must have a complete IMGT-style name,
# not merely the right prefix; IGH constant genes use their isotype letter
# instead of C.
.gene_locus_options <- function(gene, segment) {
  if (segment == "V" &&
      grepl("^TR[AD]V[0-9][A-Z0-9/-]*$", gene)) {
    return(c("TRA", "TRD"))
  }
  ordinary_vj <- grepl(
    paste0("^(IG[HKL]|TR[ABDG])", segment,
           "(?=[A-Z0-9/-]*[0-9])[A-Z0-9/-]+$"),
    gene, perl = TRUE
  )
  special_vj <- (segment == "V" && gene == "TRGVA") ||
    (segment == "J" && gene == "TRGJP")
  if (segment %in% c("V", "J") && (ordinary_vj || special_vj)) {
    return(substr(gene, 1L, 3L))
  }
  if (segment == "D") {
    # A true D segment always carries a number straight after the D, which is
    # what separates IGHD3-10 (a heavy-chain D gene) from IGHD (the delta
    # constant gene) and covers the IGHD1/OR15-1a orphans as well.
    if (grepl("^IGHD[0-9]", gene)) return("IGH")
    if (grepl("^TRBD[0-9]+$", gene)) return("TRB")
    if (grepl("^TRDD[0-9]+$", gene)) return("TRD")
    return(character(0))
  }
  if (segment == "C") {
    # IGHD is both the delta constant gene and the prefix of true IGH D
    # segments, so heavy-chain constants need a whole-token match. The pinned
    # human IMGT C database also contains IGHGP and the IGLJ-C/OR18 orphan gene.
    if (grepl("^IGH(?:M|D|E|A[12]|G[1-4]A?|GP)$", gene, perl = TRUE)) {
      return("IGH")
    }
    if (gene == "IGLJ-C/OR18") return("IGL")
    if (gene == "IGKC" || grepl("^IGLC[0-9]+[A-Z]?$", gene)) return(substr(gene, 1L, 3L))
    if (gene %in% c("TRAC", "TRDC") || grepl("^TR[BG]C[0-9]+$", gene)) {
      return(substr(gene, 1L, 3L))
    }
  }
  character(0)
}

#' Test whether every gene in an AIRR call is compatible with a locus/segment.
#'
#' Missing, malformed, wrong-segment and cross-locus calls return FALSE. The
#' shared TRAV/DV genes are valid for either TRA or TRD when used as V calls.
#' A missing D call is normal rather than incoherent, so callers testing D
#' should treat absence separately (see `airr_qc_flags()`).
gene_call_locus_compatible <- function(call, locus, segment = c("V", "J", "D", "C")) {
  segment <- match.arg(segment)
  n <- max(length(call), length(locus))
  if (!length(call) || !length(locus)) return(logical(0))
  if (!length(call) %in% c(1L, n) || !length(locus) %in% c(1L, n)) {
    stop("call and locus must have equal lengths or length one")
  }
  call <- rep_len(call, n)
  locus <- rep_len(as.character(locus), n)
  vapply(seq_len(n), function(i) {
    genes <- .gene_set(call[i])
    length(genes) > 0L && !is.na(locus[i]) && nzchar(locus[i]) &&
      all(vapply(genes, function(gene) {
        locus[i] %in% .gene_locus_options(gene, segment)
      }, logical(1)))
  }, logical(1), USE.NAMES = FALSE)
}

#' Sort a call's distinct gene names without changing duplicated-location names.
parse_vj_call_set <- function(x) {
  vapply(x, function(cell) paste(sort(.gene_set(cell)), collapse = ","),
         character(1), USE.NAMES = FALSE)
}

#' Classify complete candidate lists before reducing named alleles to genes.
#'
#' Within-gene allele ties resolve at gene level. A multi-gene call qualifies
#' only when every named allele has the same complete, unambiguous reference
#' sequence class. identity_classes maps original allele names to those classes;
#' read_germline_identity_classes() in reannotated_airr.R reads the reference table.
#' Complete candidates means IgBLAST's reporting limit did not leave the set of
#' equally best hits potentially truncated. attach_candidate_completeness() in
#' reannotated_airr.R checks this separately from reference-sequence completeness.
#' An incomplete list cannot establish either result. Return one status per call:
#' missing, incomplete_candidates, resolved, assigned_group or ambiguous.
vj_tie_kind <- function(x, identity_classes, candidates_complete = TRUE) {
  complete <- rep_len(candidates_complete, length(x))
  vapply(seq_along(x), function(i) {
    alleles <- unique(.split_calls(x[i]))
    if (!length(alleles)) return("missing")
    if (!isTRUE(complete[i])) return("incomplete_candidates")
    genes <- unique(canonical_gene_name(alleles))
    if (length(genes) == 1L) return("resolved")
    classes <- unname(identity_classes[alleles])
    if (all(!is.na(classes) & nzchar(classes)) && length(unique(classes)) == 1L) {
      "assigned_group"
    } else "ambiguous"
  }, character(1), USE.NAMES = FALSE)
}

#' Name a single gene or a complete-sequence composite as a distinct entity.
#'
#' A composite does not join any single-gene member. Its sequence-derived token
#' also keeps different composites apart when their gene labels overlap. Return
#' one token per call, or an empty string for missing, incomplete or ambiguous
#' calls. Grouping eligibility also checks the status from vj_tie_kind().
vj_canonical_token <- function(x, identity_classes, candidates_complete = TRUE) {
  kinds <- vj_tie_kind(x, identity_classes, candidates_complete)
  vapply(seq_along(x), function(i) {
    if (kinds[i] == "resolved") return(.gene_set(x[i])[1L])
    if (kinds[i] == "assigned_group") return(unname(identity_classes[.split_calls(x[i])[1L]]))
    ""
  }, character(1), USE.NAMES = FALSE)
}

#' Describe the analytical entity without replacing its stable identifier.
vj_entity_label <- function(x, identity_classes, candidates_complete = TRUE) {
  tokens <- vj_canonical_token(x, identity_classes, candidates_complete)
  labels <- attr(identity_classes, "labels")
  matched <- if (length(labels)) unname(labels[tokens]) else rep(NA_character_, length(tokens))
  ifelse(!is.na(matched), matched, tokens)
}

#' Flag multi-gene ambiguity and potentially truncated candidate lists.
vj_calls_ambiguous <- function(v_gene_status, j_gene_status) {
  as.character(v_gene_status) %in% c("ambiguous", "incomplete_candidates") |
    as.character(j_gene_status) %in% c("ambiguous", "incomplete_candidates")
}

#' Give the pinned inference tools one distinct parseable call per V/J entity.
#'
#' Alakazam strips duplicated-location markers and parses comma-separated calls.
#' These temporary calls keep the pipeline's exact entity boundaries in both
#' distance fitting and clone assignment while preserving every original call.
igh_inference_calls <- function(records) {
  require_columns(records, c("v_gene", "j_gene"), "IGH inference representatives")
  for (segment in c("v", "j")) {
    entity <- as.character(records[[paste0(segment, "_gene")]])
    if (any(is_airr_null(entity))) stop("IGH inference requires resolved V/J entities")
    prefix <- if (segment == "v") "IGHV1-" else "IGHJ"
    records[[paste0(segment, "_inference_call")]] <-
      paste0(prefix, match(entity, sort(unique(entity))), "*01")
  }
  records
}

#' Report whether a constant-region call can accompany a grouped record.
#'
#' A missing C call remains eligible. A present call must match the locus and
#' have complete candidates; failures return named reasons joined by semicolons.
#' This check applies in both grouping branches. It does not qualify isotypes:
#' c_call_qc_flags() in isotype_calls.R applies the separate alignment evidence gate.
grouping_constant_status <- function(records) {
  call <- column_or_default(records, "c_call")
  complete <- if ("c_call_candidates_complete" %in% names(records)) {
    parse_airr_boolean(records$c_call_candidates_complete) == "true"
  } else rep(TRUE, nrow(records))
  present <- !is_airr_null(call)
  locus_failed <- present & !gene_call_locus_compatible(call, records$locus, "C")
  incomplete <- present & !complete
  vapply(seq_len(nrow(records)), function(i) {
    reasons <- c(if (locus_failed[i]) "c_locus_conflict",
                 if (incomplete[i]) "c_call_incomplete_candidates")
    if (length(reasons)) paste(reasons, collapse = ";") else "eligible"
  }, character(1))
}

#' Return the first listed V/J gene as a descriptive label.
#'
#' Equal-score calls have no preferred first gene, so this label does not define
#' grouping. vj_canonical_token() names an accepted gene assignment;
#' parse_vj_call_set() returns the sorted set when listing order should not matter.
#' Missing calls return ""; an empty vector returns character(0).
first_vj_gene <- function(x) {
  vapply(x, function(cell) {
    genes <- .gene_set(cell)
    if (length(genes)) genes[1L] else ""
  }, character(1), USE.NAMES = FALSE)
}

#' Return the first V/J call exactly as IgBLAST wrote it, allele included.
#'
#' Unlike `first_vj_gene()`, the allele suffix is kept, because this answers a
#' different question: which single sequence a germline was reconstructed from.
#' Change-O's `CreateGermlines.py` takes the first allele of `v_call` and builds
#' `germline_alignment_d_mask` from that one reference sequence, so this is the
#' allele every mutation count is measured against. Empty input returns "".
first_vj_allele <- function(x) {
  vapply(x, function(cell) {
    alleles <- .split_calls(cell)
    if (length(alleles)) alleles[1L] else ""
  }, character(1), USE.NAMES = FALSE)
}

#' How resolved a call is for the purpose of reconstructing one germline.
#'
#' Returns one of "resolved" (a single allele), "ambiguous_multi_allele" (one gene
#' named with several alleles), "ambiguous_multi_gene" (more than one gene) or
#' "missing". Vectorised.
#'
#' This is a finer question than `vj_gene_status()` answers, and the two are
#' separate. Grouping resolves a same-gene multi-allele call at gene level, while
#' reconstructing a germline selects one reference allele. Differences between
#' candidate alleles can affect mutation counts. These descriptive statuses do
#' not filter records; grouping eligibility uses `vj_tie_kind()` instead.
#'
#' Read "resolved" as "IgBLAST reported one candidate", not as "one candidate
#' exists": this status describes the window the aligner reported through and not
#' the sequence, for the reason `vj_tie_kind()` sets out.
vj_call_ambiguity_status <- function(x) {
  vapply(x, function(cell) {
    alleles <- .split_calls(cell)
    if (!length(alleles)) return("missing")
    if (length(.gene_set(cell)) > 1L) return("ambiguous_multi_gene")
    if (length(unique(alleles)) > 1L) return("ambiguous_multi_allele")
    "resolved"
  }, character(1), USE.NAMES = FALSE)
}

#' Return the genes every one of several V/J calls admits.
#'
#' This asks whether a group already defined some other way has a single gene
#' consistent with all of it. An empty result means no one gene is named by every
#' member.
#'
#' It runs the opposite way to the two counts below it, which tally everything
#' the members named between them: one member's ambiguous call raises those and
#' cannot raise this. A member naming no gene at all empties the intersection,
#' because a member that says nothing cannot confirm a gene, so a group holding an
#' unannotated member has no shared gene.
#'
#' `x` is a character vector of AIRR V or J calls. Returns the shared genes joined
#' by commas in sorted order, or an empty string where there are none and for an
#' empty input.
intersect_vj_call_sets <- function(x) {
  if (!length(x)) return("")
  genes <- Reduce(intersect, lapply(x, .gene_set))
  genes <- sort(genes)
  if (length(genes)) paste(genes, collapse = ",") else ""
}

#' Count the distinct V/J gene families several calls name between them.
#'
#' A family is the part of an IMGT gene name before its first hyphen, so IGHV3-23
#' and IGHV3-30 are both IGHV3 while IGHV1-2 is not, and IGKV1/OR2-0 is
#' IGKV1/OR2. It takes the family level alakazam::getFamily() reports, restated
#' here because this module loads in an image that does not carry alakazam; the
#' two are not checked against each other. Genes IgBLAST could not separate are
#' paralogs of one another and usually share a family, so a count above one says
#' the calls disagree about more than which duplicated gene was seen.
#'
#' A name with no hyphen has no family level below itself and counts as its own,
#' which is what heavy-chain J genes look like. That is why the released column
#' counts V families rather than J ones.
#'
#' `x` is a character vector of AIRR V or J calls. Returns the number of distinct
#' families, and zero where no call names a gene.
count_vj_gene_families <- function(x) {
  genes <- unique(unlist(lapply(x, .gene_set), use.names = FALSE))
  length(unique(sub("-.*$", "", genes)))
}

#' Count the distinct V/J genes several calls name between them.
#'
#' The gene-level counterpart of count_vj_gene_families(), over the same union of the
#' calls' canonical gene sets. It is what separates two members called to different
#' genes of one family (IGHV3-23 against IGHV3-30, which the family count reports as
#' one) from members that agree.
#'
#' Like the family count it is a union, so it cannot on its own tell a genuine
#' disagreement between members from a single member's ambiguous call naming
#' several genes. Read it beside the shared-gene flag, which is the intersection
#' and answers that question directly.
#'
#' `x` is a character vector of AIRR V or J calls. Returns the number of distinct
#' genes, and zero where no call names one.
count_vj_genes <- function(x) {
  length(unique(unlist(lapply(x, .gene_set), use.names = FALSE)))
}

#' Count whether each call names one gene, several genes or no gene.
#'
#' Return resolved, ambiguous_multi_gene or missing. Alleles of the same gene
#' count as one gene. Representative ranking uses this to break support ties.
#' This does not decide grouping: vj_tie_kind() checks whether tied genes share
#' identical complete reference sequences. vj_call_ambiguity_status() separately
#' describes allele ambiguity for germline reconstruction.
vj_gene_status <- function(x) {
  vapply(x, function(cell) {
    g <- .gene_set(cell)
    if (!length(g)) "missing" else if (length(g) > 1L) "ambiguous_multi_gene" else "resolved"
  }, character(1), USE.NAMES = FALSE)
}
