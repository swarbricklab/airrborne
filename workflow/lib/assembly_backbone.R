# Prepare the sequence rows used for grouping after common QC.
# TRUST4 can report several CDR3 variants from one assembled contig, called an
# assembly backbone here. A representative is the selected row used in later
# analysis; it retains links to the original rows and their recovered read support.
# prepare_assembly_mapping() chooses backbone representatives, then merges those
# with identical full DNA within each library and locus. These counts describe
# recovered observations, not cells. The helpers use base R.

#' Parse each distinct call string once, then return results in input order.
#'
#' The canonical gene-status parser walks one call string at a time, and a
#' library's frame repeats a small set of distinct V/J calls across many rows.
#' Parsing the distinct values and expanding the result is the same result
#' without parsing repeated strings again.
.by_distinct <- function(x, parser) {
  values <- as.character(x)
  distinct <- unique(values)
  parser(distinct)[match(values, distinct)]
}

#' Rank rows so the first row of each backbone is its representative.
#'
#' This is the vectorised form of `select_rearrangement_group_representative()`:
#' the same preference ladder (read support, then unambiguous, then V identity,
#' then usable alignment length, then junction, full sequence and V, J and C
#' calls, then `sequence_id`), evaluated once for the whole frame instead of once
#' per group. Calling that function per backbone repeats work. Both
#' implementations must choose the same row; the backbone tests compare their
#' results.
#'
#' The content keys come before `sequence_id` because TRUST4 numbers the CDR3
#' variants of a contig differently between runs. `sequence_id` decides only
#' between rows with identical content.
#'
#' Missing columns do not affect the ranking, matching that function, with one
#' exception: `sequence_alignment` falls back to the junction, and siblings of one
#' backbone can carry different junctions, so the length criterion then ranks on
#' junction length. Common QC requires `sequence_alignment`, so production never
#' takes that path.
.backbone_representative_order <- function(df, group_key) {
  n <- nrow(df)
  col <- function(name, default) if (name %in% names(df)) df[[name]] else rep_len(default, n)
  junc <- as.character(col("junction", NA_character_))
  unamb <- !is.na(junc) & grepl("^[ACGTacgt]+$", junc) &
    .by_distinct(col("v_call", NA_character_), vj_gene_status) == "resolved" &
    .by_distinct(col("j_call", NA_character_), vj_gene_status) == "resolved"
  support <- if ("consensus_count" %in% names(df)) {
    as_read_support(df$consensus_count, "consensus_count")
  } else {
    rep(0, n)
  }
  identity <- suppressWarnings(as.numeric(col("v_identity", NA)))
  identity[is.na(identity)] <- -Inf
  usable_len <- nchar(as.character(col("sequence_alignment", junc)))
  usable_len[is.na(usable_len)] <- 0L
  content <- lapply(c("sequence", "v_call", "j_call", "c_call"),
                    function(name) as.character(col(name, NA_character_)))
  seqid <- as.character(col("sequence_id", seq_len(n)))
  do.call(order, c(list(group_key, -support, -as.integer(unamb), -identity,
                        -usable_len, junc), content, list(seqid)))
}

#' Find the assembled contig named by each TRUST4 AIRR sequence_id.
#'
#' TRUST4 adds a trailing _<CDR3 index> for each reported CDR3 variant. Remove
#' only that numeric suffix to recover the contig ID used by _final.out and
#' _annot.fa. A pooled <library_id>~assemble<N>_<M> keeps its library prefix.
#' This assumes TRUST4 converter IDs; another producer's <name>_<number> IDs
#' could merge unrelated rows and require a different naming contract.
trust4_assembly_backbone_id <- function(sequence_id) {
  sub("_[0-9]+$", "", as.character(sequence_id))
}

#' Map each AIRR row to the representative row of its assembly backbone.
#'
#' Rows with the same contig ID share a representative even when their junctions
#' or V/J annotations differ. This groups reported assembly variants; it does not
#' establish that the variants came from one biological rearrangement.
#' Choose the row with greatest original read support, then use
#' .backbone_representative_order() to break ties. The caller applies eligibility
#' separately. Return one row index per df row; a representative maps to itself.
#' Empty input returns integer(0).
trust4_assembly_representative_index <- function(df) {
  if (nrow(df) == 0L) return(integer(0))
  # V/J calls break ties between records with equal support.
  required <- c("sequence_id", "v_call", "j_call", "consensus_count")
  missing <- setdiff(required, names(df))
  if (length(missing)) {
    stop("TRUST4 assembly backbone mapping requires: ",
         paste(missing, collapse = ", "))
  }

  backbone <- trust4_assembly_backbone_id(df$sequence_id)
  # Numbering the backbones makes their order independent of string sorting.
  group_key <- match(backbone, unique(backbone))

  ord <- .backbone_representative_order(df, group_key)
  winner_row <- ord[!duplicated(group_key[ord])]
  winner_by_group <- integer(max(group_key))
  winner_by_group[group_key[winner_row]] <- winner_row
  winner_by_group[group_key]
}

GROUPING_COLUMNS <- c(
  "library_id", "sequence_id", "locus", "assembly_backbone_id",
  "grouping_eligible", "grouping_status",
  "backbone_representative_sequence_id", "is_backbone_representative",
  "inference_representative_sequence_id",
  "is_inference_representative", "inference_consensus_count"
)

make_empty_grouping <- function() {
  data.frame(
    library_id = character(), sequence_id = character(), locus = character(),
    assembly_backbone_id = character(), grouping_eligible = logical(),
    grouping_status = character(), backbone_representative_sequence_id = character(),
    is_backbone_representative = logical(), inference_representative_sequence_id = character(),
    is_inference_representative = logical(), inference_consensus_count = numeric(),
    stringsAsFactors = FALSE
  )
}

#' Choose the dominant common-QC BCR variant before grouping eligibility, then
#' merge identical full sequences among retained representatives within each
#' library/locus. Keep both links and every source row; transfer each eligible
#' source's support once.
#'
#' `records` are one library's common-QC rows, `library_id` the library they
#' must all carry, `eligible` one logical per row from the caller's grouping
#' rule and `status` the matching `grouping_status` text. `collapse` is TRUE
#' where eligible variants of one assembly backbone collapse to one
#' representative (IGH, IGK and IGL) and FALSE where every eligible variant is
#' its own backbone representative (TCR loci).
#'
#' Within each locus, the backbone representative is chosen from all common-QC
#' variants by `.backbone_representative_order()` on original record support.
#' An ineligible dominant variant excludes the whole backbone. Eligible minor
#' siblings never replace it. Backbone
#' representatives with identical full `sequence` then merge into one inference
#' representative, chosen by the same ladder on each backbone representative's
#' own original support. Merged rows must agree on each of `junction`,
#' `junction_aa`, `v_call`, `j_call`, `v_gene`, `j_gene`, `v_gene_status`,
#' `j_gene_status`, `c_call`, `sequence_alignment` and
#' `germline_alignment_d_mask` that the input carries, or the stage stops.
#' `inference_consensus_count` sums every eligible source record's original
#' support onto its inference representative once and is zero elsewhere.
#' Ineligible rows keep no links.
#'
#' Returns one row per input row with GROUPING_COLUMNS, validated by
#' validate_grouping_mapping(). Empty input returns make_empty_grouping().
#' Callers own eligibility: build_rearrangement_group_table() in
#' rearrangement_groups.R for non-IGH, prepare_igh_grouping() in igh_grouping.R
#' for IGH. Neither call replaces the earlier airr_qc_flags() record gate.
prepare_assembly_mapping <- function(records, library_id, eligible, status, collapse) {
  require_columns(records, c("library_id", "sequence_id", "sequence", "locus", "consensus_count"),
                  "Common-QC grouping input")
  if (!is.character(library_id) || length(library_id) != 1L ||
      is_airr_null(library_id)) stop("Grouping requires one nonempty string library_id")
  if (any(is_airr_null(records$sequence_id)) || anyDuplicated(records$sequence_id)) {
    stop("Grouping requires unique, nonempty source sequence_id values")
  }
  if (any(is_airr_null(records$sequence))) stop("Grouping requires complete sequence strings")
  if (any(is.na(records$library_id) | records$library_id != library_id)) stop("Input belongs to another library")
  if (anyNA(records$locus) ||
      any(!records$locus %in% c("IGH", "IGK", "IGL", "TRA", "TRB", "TRD", "TRG"))) {
    stop("Grouping requires one of the seven receptor loci")
  }
  support <- as_read_support(records$consensus_count, "consensus_count")
  if (any(support <= 0)) stop("Common-QC input contains non-positive read support")
  if (!is.logical(eligible) || anyNA(eligible) || length(eligible) != nrow(records) ||
      length(status) != nrow(records)) stop("Invalid prepared eligibility")
  if (!nrow(records)) return(make_empty_grouping())
  out <- data.frame(
    library_id = as.character(records$library_id), sequence_id = as.character(records$sequence_id),
    locus = as.character(records$locus),
    assembly_backbone_id = trust4_assembly_backbone_id(records$sequence_id),
    grouping_eligible = eligible, grouping_status = status,
    backbone_representative_sequence_id = NA_character_, is_backbone_representative = FALSE,
    inference_representative_sequence_id = NA_character_,
    is_inference_representative = FALSE, inference_consensus_count = 0,
    stringsAsFactors = FALSE
  )
  for (locus in unique(out$locus)) {
    locus_rows <- which(out$locus == locus)
    dominant <- if (collapse) {
      locus_rows[trust4_assembly_representative_index(records[locus_rows, , drop = FALSE])]
    } else locus_rows
    blocked <- eligible[locus_rows] & !eligible[dominant]
    if (any(blocked)) {
      eligible[locus_rows[blocked]] <- FALSE
      out$grouping_eligible[locus_rows[blocked]] <- FALSE
      out$grouping_status[locus_rows[blocked]] <- "ineligible_backbone_representative"
    }
    at <- locus_rows[eligible[locus_rows]]
    if (!length(at)) next
    backbone_links <- dominant[eligible[locus_rows]]
    out$backbone_representative_sequence_id[at] <- out$sequence_id[backbone_links]
    backbone_rows <- unique(backbone_links)
    out$is_backbone_representative[backbone_rows] <- TRUE
    exact <- match(records$sequence[backbone_rows], unique(records$sequence[backbone_rows]))
    annotation_columns <- intersect(c("junction", "junction_aa", "v_call", "j_call",
      "v_gene", "j_gene", "v_gene_status", "j_gene_status", "c_call",
      "sequence_alignment", "germline_alignment_d_mask"), names(records))
    for (group in split(backbone_rows, exact)) {
      if (length(group) < 2L) next
      for (column in annotation_columns) {
        if (length(unique(records[[column]][group])) != 1L) {
          stop("Identical sequences disagree in ", column, ": ",
               paste(records$sequence_id[group], collapse = ", "))
        }
      }
    }
    ranked <- .backbone_representative_order(records[backbone_rows, , drop = FALSE], exact)
    winners <- ranked[!duplicated(exact[ranked])]
    sequence_links <- backbone_rows[winners[match(exact, exact[winners])]]
    links <- sequence_links[match(backbone_links, backbone_rows)]
    out$inference_representative_sequence_id[at] <- out$sequence_id[links]
    out$is_inference_representative[unique(links)] <- TRUE
    totals <- vapply(split(support[at], links), sum, numeric(1))
    out$inference_consensus_count[as.integer(names(totals))] <- unname(totals)
  }
  validate_grouping_mapping(records, out)
}

# Validate links and support against the source pool without reapplying gates,
# ranking representatives or collapsing observations. Returns `mapping` with its
# three flags parsed to logical and `inference_consensus_count` made numeric.
validate_grouping_mapping <- function(records, mapping) {
  require_columns(records, c("library_id", "sequence_id", "sequence", "locus", "consensus_count"),
                  "Prepared grouping source records")
  require_columns(mapping, GROUPING_COLUMNS, "Prepared grouping")
  if (nrow(mapping) != nrow(records) || anyDuplicated(mapping$sequence_id) ||
      anyDuplicated(records$sequence_id) ||
      !identical(as.character(mapping$sequence_id), as.character(records$sequence_id))) {
    stop("Prepared grouping must have identical ordered source keys")
  }
  for (field in c("library_id", "locus")) {
    if (!identical(as.character(records[[field]]), as.character(mapping[[field]]))) {
      stop("Prepared grouping disagrees on ", field)
    }
  }
  for (field in c("grouping_eligible", "is_backbone_representative", "is_inference_representative")) {
    parsed <- parse_airr_boolean(mapping[[field]])
    if (any(!parsed %in% c("true", "false"))) stop("Incomplete prepared grouping flags")
    mapping[[field]] <- parsed == "true"
  }
  mapping$inference_consensus_count <- as_read_support(mapping$inference_consensus_count)
  eligible <- mapping$grouping_eligible
  representative <- mapping$is_inference_representative
  links <- match(mapping$inference_representative_sequence_id, mapping$sequence_id)
  backbone <- mapping$is_backbone_representative
  backbone_links <- match(mapping$backbone_representative_sequence_id, mapping$sequence_id)
  if (any(backbone & !eligible) || any(representative & !backbone) ||
      anyNA(backbone_links[eligible]) || anyNA(links[eligible]) ||
      any(!is_airr_null(mapping$inference_representative_sequence_id[!eligible])) ||
      any(!is_airr_null(mapping$backbone_representative_sequence_id[!eligible])) ||
      any(mapping$inference_consensus_count[!representative] != 0)) {
    stop("Invalid prepared representative links or support")
  }
  if (any(eligible)) {
    at <- which(eligible)
    to <- links[at]
    via <- backbone_links[at]
    if (any(!representative[to]) ||
        any(!backbone[via]) ||
        any(mapping$library_id[at] != mapping$library_id[to]) ||
        any(mapping$locus[at] != mapping$locus[to]) ||
        any(mapping$library_id[at] != mapping$library_id[via]) ||
        any(mapping$locus[at] != mapping$locus[via]) ||
        any(mapping$assembly_backbone_id[at] != mapping$assembly_backbone_id[via]) ||
        any(records$sequence[via] != records$sequence[to]) ||
        any(links[via] != to) || any(backbone_links[backbone] != which(backbone)) ||
        any(links[representative] != which(representative))) stop("Invalid prepared representative context")
    tcr <- !mapping$locus[at] %in% c("IGH", "IGK", "IGL")
    if (any(via[tcr] != at[tcr])) stop("TCR backbone variants must remain distinct")
    keys <- paste(mapping$library_id[representative], mapping$locus[representative],
                  records$sequence[representative], sep = "\x1f")
    if (anyDuplicated(keys)) stop("Inference representatives must have distinct full sequences")
    totals <- vapply(split(as_read_support(records$consensus_count, "consensus_count")[at], to), sum, numeric(1))
    if (any(mapping$inference_consensus_count[as.integer(names(totals))] != totals)) {
      stop("Prepared representatives did not preserve eligible support")
    }
  }
  mapping
}
