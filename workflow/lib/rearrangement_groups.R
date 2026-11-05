# Build non-IGH clonotype tables from records that passed common QC.
# build_rearrangement_group_table() defines membership and support below.
# These helpers use base R.

#' Give each distinct clonotype_key a label within one library and locus.
#'
#' Sort the distinct keys and number them as <library_id>_<locus>_ct<NNNNN>.
#' The library and locus make the label unique across a release. Reordering
#' input rows leaves the labels unchanged. Missing keys get NA.
#' set_pipeline_collation() fixes the string sorting order across runs.
#' Cross-library overlap compares clonotype_key, which records the junction and
#' V/J assignments, rather than this numbered label.
make_rearrangement_group_id <- function(library_id, locus, clonotype_key) {
  key  <- as.character(clonotype_key)
  uniq <- sort(unique(key[!is.na(key)]))
  idx  <- match(key, uniq)
  out  <- paste0(library_id, "_", locus, "_ct", sprintf("%05d", idx))
  out[is.na(key)] <- NA_character_
  out
}

#' Choose one common-QC group member for descriptive sequence annotations.
#'
#' Rank by decreasing read support, unambiguous junction and V/J calls, V identity
#' and alignment length, then increasing junction, full sequence, V, J and C
#' calls and sequence_id. Missing optional columns make that criterion neutral.
#' Content precedes sequence_id because TRUST4 numbers CDR3 variants differently
#' between runs. Return the row index and ranking reason; an empty group returns
#' integer(0) and a missing reason.
#'
#' Support leads because tied gene calls do not reduce the reads supporting a
#' record. All candidates already passed common QC, so productivity adds no rank.
select_rearrangement_group_representative <- function(df) {
  n <- nrow(df)
  if (n == 0L) return(list(index = integer(0), reason = NA_character_))
  # rep_len (not rep) so a length-n vector default (e.g. sequence_alignment -> junc) is
  # passed through unchanged rather than tiled to length n*n.
  col <- function(name, default) if (name %in% names(df)) df[[name]] else rep_len(default, n)
  junc       <- as.character(col("junction", NA_character_))
  unamb      <- !is.na(junc) & grepl("^[ACGTacgt]+$", junc) &
                vj_gene_status(col("v_call", NA_character_)) == "resolved" &
                vj_gene_status(col("j_call", NA_character_)) == "resolved"
  support <- if ("consensus_count" %in% names(df)) {
    as_read_support(df$consensus_count, "consensus_count")
  } else {
    rep(0, n)
  }
  identity   <- suppressWarnings(as.numeric(col("v_identity", NA)))
  identity[is.na(identity)] <- -Inf
  usable_len <- nchar(as.character(col("sequence_alignment", junc)))
  usable_len[is.na(usable_len)] <- 0L
  content    <- lapply(c("sequence", "v_call", "j_call", "c_call"),
                       function(name) as.character(col(name, NA_character_)))
  seqid      <- as.character(col("sequence_id", seq_len(n)))
  ord <- do.call(order, c(list(-support, -as.integer(unamb), -identity, -usable_len, junc),
                          content, list(seqid)))
  win <- ord[1]
  reason <- sprintf(
    "consensus_count=%s;tie=unambiguous>v_identity>len>content>seq_id", support[win])
  list(index = win, reason = reason)
}

# The released columns of a locus's clonotype membership, in order.
# member_consensus_count is the member's transferred support, the value the AIRR
# export carries in the same column.
CLONOTYPE_MEMBER_COLUMNS <- c(
  "library_id", "locus", "clonotype_id", "sequence_id", "member_consensus_count")

# The released columns of a locus's clonotype table, in order. v_gene_key and
# j_gene_key are the analytical V and J entities the key holds; n_contigs
# counts the TRUST4 contigs, called assembly backbones in the grouping code.
REARRANGEMENT_GROUP_COLUMNS <- c(
  "library_id", "locus", "clonotype_id", "clonotype_key",
  "v_gene_key", "j_gene_key", "junction", "junction_aa", "junction_length",
  "n_contigs", "n_members", "consensus_count", "support_fraction",
  "representative_sequence_id"
)

#' The rearrangement-group table of a locus that produced no eligible group.
#'
#' Typed rather than make_empty_frame()'s all-character, because the counts and
#' the support fraction are numbers everywhere else this table appears.
make_empty_rearrangement_group_table <- function() {
  data.frame(
    library_id = character(),
    locus = character(),
    clonotype_id = character(),
    clonotype_key = character(),
    v_gene_key = character(),
    j_gene_key = character(),
    junction = character(),
    junction_aa = character(),
    junction_length = integer(),
    n_contigs = integer(),
    n_members = integer(),
    consensus_count = numeric(),
    support_fraction = numeric(),
    representative_sequence_id = character(),
    stringsAsFactors = FALSE
  )
}

#' Group one library's common-QC non-IGH records by exact junction and V/J entity.
#'
#' records contains one AIRR row per reported CDR3 with original read support.
#' junction is the nucleotide sequence spanning the receptor's V-to-J junction;
#' junction_aa is its amino-acid translation. The group requires identical A/C/G/T
#' junction and resolved v_gene and j_gene assignments. Different nucleotide
#' sequences stay separate even when they encode the same amino-acid sequence.
#' The amino-acid column describes a representative; it does not define membership.
#' gene_calls.R owns V/J resolution; rearrangement_keys.R owns the shared key.
#' Conflicting entities remain separate even when their junctions match. Sharing
#' one possible gene does not establish identical assignments: two calls can each
#' overlap a third call without sharing any gene with each other.
#'
#' IGK and IGL choose the dominant common-QC backbone variant before eligibility;
#' an ineligible dominant variant has no minor substitute. TCR retains eligible
#' variants. All loci then merge identical full DNA before junction/V/J grouping.
#' prepare_assembly_mapping() in assembly_backbone.R owns this preparation and
#' support transfer. select_rearrangement_group_representative() then chooses
#' the group's descriptive annotations from the retained sequences.
#' Incompatible or incomplete present C calls also exclude records. Their support
#' stays in common QC, outside the grouped denominator.
#'
#' groups follows REARRANGEMENT_GROUP_COLUMNS and uses summed recovered read
#' support. Sequence annotations describe one representative; V/J tokens describe
#' the shared key. clonotype_id names one library and locus and is unique across
#' a release; clonotype_key compares libraries.
#' representative_sequence_id is the representative's library-local sequence_id,
#' which build_clonotypes.R qualifies with the library as it qualifies membership.
#' mapping records every input row's grouping decision. Counts report input,
#' dropped and eligible rows. n_dropped_ambiguous_vj overlaps other failures and
#' is a subset of n_dropped. Empty input returns typed empty groups and mapping.
build_rearrangement_group_table <- function(records, library_id, locus) {
  require_columns(records, c("sequence_id", "junction", "v_call", "j_call",
                             "v_gene_status", "j_gene_status",
                             "v_gene", "j_gene",
                             "consensus_count"),
                  "common QC-pass records")
  if (!"library_id" %in% names(records)) records$library_id <- rep(library_id, nrow(records))
  if (!"locus" %in% names(records)) records$locus <- rep(locus, nrow(records))
  if (any(records$locus != locus)) stop("Non-IGH grouping input has another locus")
  # Non-positive read support is a contract violation rather than a group to
  # drop: the common QC gate upstream has already removed those records.
  records$consensus_count <- as_read_support(
    records$consensus_count, "consensus_count")
  if (any(records$consensus_count <= 0)) {
    stop(sprintf("[%s] common QC input contains non-positive read support", locus))
  }

  records$clonotype_key <- exact_rearrangement_group_key(
    records$junction, records$v_gene, records$j_gene)
  junction_ok <- !is_airr_null(records$junction) & grepl("^[ACGTacgt]+$", records$junction)
  v_ok <- records$v_gene_status %in% c("resolved", "assigned_group") & !is_airr_null(records$v_gene)
  j_ok <- records$j_gene_status %in% c("resolved", "assigned_group") & !is_airr_null(records$j_gene)
  c_status <- grouping_constant_status(records)
  c_ok <- c_status == "eligible"
  eligible <- !is.na(records$clonotype_key) & v_ok & j_ok & c_ok
  status <- vapply(seq_len(nrow(records)), function(i) {
    failures <- c(if (!junction_ok[i]) "ineligible_junction",
      if (!v_ok[i]) paste0("v_call_", records$v_gene_status[i]),
      if (!j_ok[i]) paste0("j_call_", records$j_gene_status[i]),
      if (!c_ok[i]) c_status[i])
    if (length(failures)) paste(failures, collapse = ";") else "eligible"
  }, character(1))
  dropped <- sum(!eligible)
  # Which of the excluded records name several germline entities, so the caller
  # can report an unanswerable call apart from an absent one or an unreadable
  # junction.
  reducible_tie <- as.character(records$v_gene_status) == "ambiguous" |
    as.character(records$j_gene_status) == "ambiguous"
  dropped_ambiguous_vj <- sum(!eligible & reducible_tie)
  # Light-chain backbones use their dominant common-QC variant before eligibility.
  # Eligible sources transfer support only when that representative is retained.
  # TCR loci keep eligible variants. Both branches then merge identical retained DNA.
  mapping <- prepare_assembly_mapping(
    records, library_id, eligible,
    status,
    collapse = locus %in% c("IGK", "IGL"))
  n_eligible <- sum(mapping$grouping_eligible)
  dropped <- nrow(mapping) - n_eligible
  records <- records[mapping$is_inference_representative, , drop = FALSE]
  records$consensus_count <- mapping$inference_consensus_count[mapping$is_inference_representative]

  # Nothing survived, so there is no group to collapse into and the code below
  # would have no rows to bind. Every locus still owes its caller the released
  # header and the counts.
  if (n_eligible == 0L) {
    return(list(groups = make_empty_rearrangement_group_table(), mapping = mapping,
                membership = make_empty_frame(CLONOTYPE_MEMBER_COLUMNS),
                n_records = length(eligible), n_dropped = dropped,
                n_dropped_ambiguous_vj = dropped_ambiguous_vj,
                n_eligible = 0L))
  }

  records$clonotype_id <- make_rearrangement_group_id(
    library_id, locus, records$clonotype_key)

  # Unlike column_or_default(), which returns a whole column as long as the frame
  # has rows, this takes the one-row representative frame and returns a single
  # scalar for the group's descriptive fields.
  value_or_na <- function(row, column) {
    if (column %in% names(row)) as.character(row[[column]]) else NA_character_
  }

  contributors <- mapping[mapping$grouping_eligible, , drop = FALSE]
  backbone_groups <- unique(data.frame(
    clonotype_key = records$clonotype_key[match(
      contributors$inference_representative_sequence_id, records$sequence_id)],
    backbone = contributors$assembly_backbone_id))
  backbone_counts <- table(backbone_groups$clonotype_key)
  groups <- split(seq_len(nrow(records)), records$clonotype_key)
  rows <- lapply(groups, function(indices) {
    members <- records[indices, , drop = FALSE]
    representative <- select_rearrangement_group_representative(members)$index
    row <- members[representative, , drop = FALSE]

    data.frame(
      library_id = library_id,
      locus = locus,
      clonotype_id = row$clonotype_id[1],
      clonotype_key = row$clonotype_key[1],
      # The one germline entity the group is keyed on, which every member names
      # by construction, so the V/J incidence panel counts what the key counts.
      v_gene_key = as.character(row$v_gene[1]),
      j_gene_key = as.character(row$j_gene[1]),
      junction = as.character(row$junction[1]),
      junction_aa = value_or_na(row, "junction_aa"),
      # The key's junction, so every member has this nucleotide length.
      junction_length = nchar(as.character(row$junction[1])),
      # How many separately assembled receptors the group was built from, as a
      # count unit alongside read support. The two answer different questions: a
      # group carrying hundreds of reads may still come from one assembly.
      # Neither count is a cell or an independent molecule. Members of one
      # TRUST4 contig are alternative spellings of its CDR3, so contigs
      # (assembly backbones) are counted, not rows.
      n_contigs = as.integer(backbone_counts[[row$clonotype_key[1]]]),
      # The distinct retained sequences the group holds, one per row of its
      # clonotype_members table, as n_members counts IGH clone members.
      n_members = length(indices),
      consensus_count = sum(members$consensus_count),
      # The member whose record supplies junction_aa, so its AIRR annotation
      # (sequence, alignments, cdr3_aa, c_call) is one join away. Group members
      # share the nucleotide junction but not necessarily the flanking contig.
      representative_sequence_id = as.character(row$sequence_id[1]),
      stringsAsFactors = FALSE
    )
  })

  # Rows are released in group-label order, which is also key order.
  groups_out <- do.call(rbind, rows)
  groups_out <- groups_out[order(groups_out$clonotype_id), , drop = FALSE]
  groups_out$support_fraction <- groups_out$consensus_count /
    sum(groups_out$consensus_count)

  membership <- data.frame(
    library_id = records$library_id, locus = records$locus,
    clonotype_id = records$clonotype_id, sequence_id = records$sequence_id,
    member_consensus_count = records$consensus_count,
    stringsAsFactors = FALSE
  )
  list(groups = groups_out[, REARRANGEMENT_GROUP_COLUMNS], mapping = mapping,
       membership = membership,
       n_records = length(eligible), n_dropped = dropped,
       n_dropped_ambiguous_vj = dropped_ambiguous_vj,
       n_eligible = n_eligible)
}
