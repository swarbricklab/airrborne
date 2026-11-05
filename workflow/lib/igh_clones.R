# Decide which common-QC IGH junctions may enter clone inference, then build
# membership from assigned sequence representatives. These helpers use base R.
# Common QC lives in record_qc.R. prepare_igh_grouping() in igh_grouping.R calls
# clone_junction_eligible() before choosing representatives for fitting and
# assignment. A failed grouping record stays in the common-QC AIRR table.
#
# assign_scope_donor_clones.R calls SCOPer separately for each donor within each
# scope, using the shared scope cutoff chosen by scope_clone_threshold.R.
# Clone membership describes recovered sequences and read support, not cells.

#' Select junctions eligible for IGH clone calibration and assignment.
#'
#' Mirrors scoper::hierarchicalClones under IUPAC = FALSE: an eligible junction
#' is non-empty, contains only uppercase ATCG plus the accepted missing symbols
#' N and ?, and has no more than `max_n` missing symbols. The assignment call does
#' not set IUPAC and relies on the package default being FALSE; a scoper upgrade
#' that changed that default would widen the clusterer's alphabet beyond this gate.
clone_junction_eligible <- function(junction, max_n = 0L) {
  if (length(max_n) != 1L || is.na(max_n) || max_n < 0L) {
    stop("max_n must be one non-negative number")
  }
  sequence <- as.character(junction)
  n_non_atcg <- nchar(gsub("[ATCG]", "", sequence))
  accepted <- !is.na(sequence) & grepl("^[ATCGN?]+$", sequence)
  accepted & nzchar(sequence) & n_non_atcg <= max_n
}

#' Select records with usable junctions and resolved V/J reference-sequence groups.
#'
#' junction is the nucleotide junction; max_n is params.lineage.max_n.
#' v_gene_status and j_gene_status come from reference-based call resolution in
#' gene_calls.R. Both must be resolved or assigned_group. Several gene names may
#' qualify only when they name the same complete reference sequence. Calls spanning
#' different reference sequences could connect separate clones, so they fail.
#' Returns one logical per record. prepare_igh_grouping() applies these tests
#' separately to record each failure reason and also checks constant calls.
#' Review that function's V/J test alongside this helper when requirements change.
clone_record_eligible <- function(junction, v_gene_status, j_gene_status,
                                    max_n = 0L) {
  clone_junction_eligible(junction, max_n) &
    v_gene_status %in% c("resolved", "assigned_group") &
    j_gene_status %in% c("resolved", "assigned_group")
}

# ---------------------------------------------------------------- assignment ---

#' Give each eligible source row its sequence representative's assigned clone.
#'
#' representative gives row indices into df; clustered contains the sequence_id
#' and clone_id returned by SCOPer for representatives. Returns one label per df
#' row. Match by sequence_id because SCOPer may reorder its result. Stop if any
#' representative is absent or any assigned label is missing or blank.
#'
#' This restores labels to source rows for the membership builder. It does not
#' cluster again or transfer support: build_clone_membership_frame() applies the
#' prepared support totals once to retained representatives.
propagate_clone_id_to_siblings <- function(df, representative, clustered) {
  representative_id <- as.character(df$sequence_id[representative])
  assigned_index <- match(representative_id, as.character(clustered$sequence_id))
  if (anyNA(assigned_index)) {
    stop("hierarchicalClones omitted an inference representative")
  }
  clone_id <- clustered$clone_id[assigned_index]
  if (any(is.na(clone_id) | as.character(clone_id) == "")) {
    stop("an eligible IGH row did not receive a propagated clone_id")
  }
  clone_id
}

# ---------------------------------------------------------------- membership ---

#' Combine scope, donor and SCOPer's local cluster label into one clone ID.
#'
#' Labels apply only within this donor's run in this scope. Another donor or scope
#' can reuse the same SCOPer label without naming the same clone. The ID also does
#' not establish identity across changed inputs or methods. Scope and donor names
#' exclude ~, keeping the three parts distinct.
make_clone_id <- function(scope_id, donor_id, clone_id) {
  paste(scope_id, donor_id, clone_id, sep = "~")
}

# The membership table every downstream IGH metric reads. It stays narrow because
# a whole scope's worth of it is loaded at once.
# `v_gene` and `j_gene` sit beside the calls they come from:
# the single germline entity each call resolves to, empty where it resolves to
# none. Gene usage votes on the entity rather than on the gene name, so two names
# for one germline read as agreement instead of as a disagreement between members.
CLONE_MEMBER_COLUMNS <- c(
  "scope_id", "donor_id", "library_id", "sequence_id", "clone_id",
  "consensus_count", "v_call", "j_call", "v_gene",
  "j_gene", "junction", "junction_aa")

# The record shard adds what a cross-sample lookup or a lineage reconstruction
# needs and a metric does not: the assembled contig, a match key that does not
# depend on which library or scope the receptor came from, and the IMGT-gapped
# observed and germline alignments a downstream phylogenetic reconstruction reads.
# `tissue_id` and the QC-gated isotype are attached later, where the scope manifest
# and the per-library IGH feature tables are both in scope.
CLONE_RECORD_COLUMNS <- c(
  "scope_id", "donor_id", "library_id", "sequence_id", "clone_id",
  "clonotype_key", "v_call", "j_call", "v_gene",
  "j_gene", "junction", "junction_aa",
  "consensus_count", "sequence",
  "junction_length", "sequence_alignment", "germline_alignment_d_mask")

#' Build assigned membership for one donor within one scope.
#'
#' cp contains eligible source rows with clone_id labels. mapping identifies each
#' retained representative and the recovered read support assigned to it during
#' preparation. Keep one row per retained full sequence within each library and
#' apply that support total once. Support from ineligible records stays in common
#' QC only; it does not enter clone abundance or downstream clone summaries.
#'
#' Returns the columns needed by both membership and sequence-record tables,
#' ordered by sequence_id. scope_id and donor_id become part of every clone ID.
build_clone_membership_frame <- function(cp, mapping, scope_id, donor_id) {
  eligible <- parse_airr_boolean(mapping$grouping_eligible) == "true"
  mapping <- mapping[eligible, , drop = FALSE]
  mapping <- validate_grouping_mapping(cp, mapping)
  cp <- cp[mapping$is_inference_representative, , drop = FALSE]
  cp$consensus_count <- mapping$inference_consensus_count[mapping$is_inference_representative]

  support <- as_read_support(cp$consensus_count, "consensus_count")

  out_df <- data.frame(
    scope_id              = scope_id,
    donor_id              = donor_id,
    library_id                = as.character(cp$library_id),
    sequence_id          = as.character(cp$sequence_id),
    clone_id            = make_clone_id(scope_id, donor_id, cp$clone_id),
    consensus_count   = support,
    v_call                = as.character(cp$v_call),
    j_call                = as.character(cp$j_call),
    v_gene      = as.character(cp$v_gene),
    j_gene      = as.character(cp$j_gene),
    junction              = as.character(cp$junction),
    junction_aa           = as.character(cp$junction_aa),
    # Record-shard column only; the membership write selects CLONE_MEMBER_COLUMNS.
    # `clonotype_key` is the exact nucleotide junction plus the V and J
    # germline entity already used to match non-IGH rearrangement groups between
    # libraries, computed here so both arms of the pipeline offer an outside dataset
    # the same key. It is NOT how IGH clones are defined. Somatic hypermutation
    # splits one clone across many exact keys, so `clone_id` stays the clonal
    # unit. It is NA where the junction is not unambiguous A/C/G/T or the record's V
    # or J call names no single germline entity.
    clonotype_key = exact_rearrangement_group_key(
      cp$junction, cp$v_gene, cp$j_gene),
    sequence              = as.character(cp$sequence),
    # Record-shard columns only. These are the receptor as IgBLAST aligned it:
    # `sequence_alignment` and `germline_alignment_d_mask` are IMGT-gapped to the
    # same length, with D and junctional positions masked in the germline.
    # Germline alignments can differ among members of a clone. A downstream
    # method requiring one clone germline must reconcile them before inference.
    # The pipeline builds no lineage tree.
    junction_length           = as.character(cp$junction_length),
    sequence_alignment        = as.character(cp$sequence_alignment),
    germline_alignment_d_mask = as.character(cp$germline_alignment_d_mask),
    stringsAsFactors = FALSE
  )

  # Deterministic row order.
  out_df[order(out_df$sequence_id), , drop = FALSE]
}
