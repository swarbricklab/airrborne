# Count V/J usage over supplied receptor units using base R. IGH units are
# library observations of donor/scope clones; other loci use exact clonotypes.
# Each assigned unit counts once in n_groups and adds its read support to
# consensus_count.

#' Read the readable label of each composite V/J entity from reference metadata.
#'
#' `path` is the reference allele metadata that annotation also reads. Returns
#' a character vector of `identity_label` values named by `identity_id`. An
#' identifier must carry one label. Two exact-sequence classes can share a label,
#' so the label never replaces the identifier as a key.
read_vj_entity_labels <- function(path) {
  metadata <- read_pipeline_table(
    path, required_columns = c("identity_id", "identity_label"),
    context = "Reference allele metadata")
  usable <- !is_airr_null(metadata$identity_id)
  pairs <- unique(metadata[usable, c("identity_id", "identity_label"), drop = FALSE])
  if (anyDuplicated(pairs$identity_id)) {
    stop("Reference metadata gives one identity_id more than one identity_label: ", path)
  }
  stats::setNames(as.character(pairs$identity_label), pairs$identity_id)
}

#' Name each V/J analytical entity for a reader.
#'
#' A composite entity takes its reference label, the sorted member genes of its
#' exact sequence class. A single gene is its own label. The empty entity of an
#' unassigned row stays empty. This matches `vj_entity_label()` in the AIRR
#' export, so one entity reads the same in every released table.
vj_entity_display_label <- function(entities, entity_labels) {
  labels <- as.character(entities)
  matched <- unname(entity_labels[labels])
  labels[!is.na(matched)] <- matched[!is.na(matched)]
  labels
}

#' Read the reference facts gene_usage reports for each V/J entity.
#'
#' `path` is the reference allele metadata that annotation also reads, and only
#' its V-REGION and J-REGION rows count. A gene equivalence group joins every
#' gene that shares a complete reference sequence with another, directly or
#' through a third gene, and is written as its sorted genes. A gene that shares
#' no sequence is its own group.
#'
#' Returns a list of three named character vectors: `labels` from
#' read_vj_entity_labels(); `equivalence_groups`, the gene equivalence group of
#' each gene and of each composite entity; and `alleles`, the sorted alleles of
#' each composite entity, which share one complete reference sequence.
read_vj_entity_reference <- function(path) {
  metadata <- read_pipeline_table(
    path, required_columns = c("allele", "gene", "region", "identity_id"),
    context = "Reference allele metadata")
  metadata <- metadata[metadata$region %in% c("V-REGION", "J-REGION"), , drop = FALSE]
  gene <- as.character(metadata$gene)
  identity <- as.character(metadata$identity_id)
  classed <- !is_airr_null(identity) & nzchar(identity)
  class_genes <- lapply(split(gene[classed], identity[classed]), unique)
  composites <- class_genes[lengths(class_genes) > 1L]
  genes <- sort(unique(gene), method = "radix")
  root <- stats::setNames(genes, genes)
  # Each pass merges every group a composite touches; passes repeat until no
  # group changes, which joins genes linked only through a third gene.
  repeat {
    before <- root
    for (members in composites) {
      touched <- names(root)[root %in% root[members]]
      root[touched] <- sort(root[touched], method = "radix")[1L]
    }
    if (identical(before, root)) break
  }
  group_names <- vapply(split(names(root), root), function(members) {
    paste(sort(members, method = "radix"), collapse = ",")
  }, character(1))
  groups <- stats::setNames(unname(group_names[root]), names(root))
  composite_groups <- vapply(composites, function(members) groups[[members[1L]]],
                             character(1))
  alleles <- vapply(names(composites), function(id) {
    paste(sort(unique(as.character(metadata$allele[classed & identity == id])),
               method = "radix"), collapse = ",")
  }, character(1))
  list(labels = read_vj_entity_labels(path),
       equivalence_groups = c(groups, composite_groups),
       alleles = alleles)
}

#' Name the gene equivalence group of each V/J entity.
#'
#' `reference` comes from read_vj_entity_reference(). A gene absent from the
#' reference is its own group. A composite entity absent from it stops the
#' stage, since its genes and alleles cannot be named.
vj_entity_equivalence_group <- function(entities, reference) {
  entities <- as.character(entities)
  groups <- unname(reference$equivalence_groups[entities])
  composite <- grepl("^[VJ]SEQ_", entities)
  if (any(is.na(groups) & composite)) {
    stop("Reference metadata lacks the composite entities ",
         paste(unique(entities[is.na(groups) & composite]), collapse = ", "))
  }
  groups[is.na(groups)] <- entities[is.na(groups)]
  groups
}

# The columns count_gene_usage() returns, which the gene_usage.tsv producers
# write after their identifier columns.
GENE_USAGE_COUNT_COLUMNS <- c(
  "gene_key", "gene_names", "gene_equivalence_group", "identical_alleles",
  "n_groups", "group_fraction", "consensus_count", "support_fraction"
)

#' Count the receptor units and the read support of each V/J entity.
#'
#' `entities` is the one germline entity each unit was assigned, as a non-IGH
#' clonotype's token or the per-clone vote of summarise_clone_vj_call_sets()
#' writes it. `support` is each unit's summed consensus_count in the library,
#' and `reference` comes from read_vj_entity_reference(). Units with the empty
#' entity were assigned no gene and are excluded, so both fractions are shares
#' of the units that were assigned one.
#'
#' Returns one row per entity with the columns in GENE_USAGE_COUNT_COLUMNS.
#' `identical_alleles` is empty for a single gene. Equivalence groups come in
#' order of their summed n_groups, largest first, then by name, and the entities
#' of a group follow the same order, so one group's rows stay together. Returns
#' no rows when nothing was assigned. The caller prepends its identifier columns.
count_gene_usage <- function(entities, support, reference) {
  entities <- as.character(entities)
  assigned <- entities != ""
  support <- as_read_support(support[assigned], "consensus_count")
  entities <- entities[assigned]
  keys <- sort(unique(entities), method = "radix")
  units <- factor(entities, levels = keys)
  n_groups <- as.integer(table(units))
  consensus_count <- as.numeric(vapply(split(support, units), sum, numeric(1)))
  alleles <- unname(reference$alleles[keys])
  alleles[is.na(alleles)] <- ""
  usage <- data.frame(
    gene_key = keys,
    gene_names = vj_entity_display_label(keys, reference$labels),
    gene_equivalence_group = vj_entity_equivalence_group(keys, reference),
    identical_alleles = as.character(alleles),
    n_groups = n_groups,
    group_fraction = n_groups / sum(n_groups),
    consensus_count = consensus_count,
    support_fraction = consensus_count / sum(consensus_count),
    stringsAsFactors = FALSE
  )
  if (!nrow(usage)) return(usage)
  group_units <- stats::ave(usage$n_groups, usage$gene_equivalence_group, FUN = sum)
  usage <- usage[order(-group_units, usage$gene_equivalence_group, -usage$n_groups,
                       usage$gene_key, method = "radix"), , drop = FALSE]
  rownames(usage) <- NULL
  usage
}

# The columns count_vj_pair_usage() returns, which the vj_pair_usage.tsv
# producers write after their identifier columns.
VJ_PAIR_USAGE_COUNT_COLUMNS <- c(
  "v_gene_key", "v_gene_names", "j_gene_key", "j_gene_names",
  "n_groups", "group_fraction", "consensus_count", "support_fraction"
)

#' Count the receptor units and the read support of each observed V and J pair.
#'
#' `v_entities`, `j_entities` and `support` describe the units count_gene_usage()
#' reads, one position per unit, and `reference` comes from
#' read_vj_entity_reference(). A composite entity counts once under its own key
#' and is not split into its genes. A unit with an empty entity at either segment
#' has no pair and is excluded, so both fractions are shares of the units with a
#' pair. When every unit has both entities, summing the rows over one segment
#' gives count_gene_usage() for the other.
#'
#' Returns one row per observed pair with the columns in
#' VJ_PAIR_USAGE_COUNT_COLUMNS, by n_groups, largest first, then by V key and J
#' key. Returns no rows when no unit has a pair. The caller prepends its
#' identifier columns.
count_vj_pair_usage <- function(v_entities, j_entities, support, reference) {
  v_entities <- as.character(v_entities)
  j_entities <- as.character(j_entities)
  paired <- v_entities != "" & j_entities != ""
  support <- as_read_support(support[paired], "consensus_count")
  v_entities <- v_entities[paired]
  j_entities <- j_entities[paired]
  pairs <- unique(data.frame(v_gene_key = v_entities, j_gene_key = j_entities,
                             stringsAsFactors = FALSE))
  units <- match(paste(v_entities, j_entities, sep = "\t"),
                 paste(pairs$v_gene_key, pairs$j_gene_key, sep = "\t"))
  n_groups <- tabulate(units, nbins = nrow(pairs))
  consensus_count <- as.numeric(vapply(
    split(support, factor(units, levels = seq_len(nrow(pairs)))), sum, numeric(1)))
  usage <- data.frame(
    v_gene_key = pairs$v_gene_key,
    v_gene_names = vj_entity_display_label(pairs$v_gene_key, reference$labels),
    j_gene_key = pairs$j_gene_key,
    j_gene_names = vj_entity_display_label(pairs$j_gene_key, reference$labels),
    n_groups = n_groups,
    group_fraction = n_groups / sum(n_groups),
    consensus_count = consensus_count,
    support_fraction = consensus_count / sum(consensus_count),
    stringsAsFactors = FALSE
  )
  usage <- usage[order(-usage$n_groups, usage$v_gene_key, usage$j_gene_key,
                       method = "radix"), , drop = FALSE]
  rownames(usage) <- NULL
  usage
}

# What a clone's vote for one gene type can come to. `assigned` names a gene;
# the other three say why none could be named. Common QC and clone inference
# rule the three out, so calculate_scope_clone_diversity.R stops on any of them.
CLONE_VOTE_STATUSES <- c(
  "assigned", "no_callable_member", "members_ambiguous",
  "resolved_members_disagree"
)

# Vote ONE clone's category for one gene type. `calls` are the AIRR V or J
# calls of its members and `tokens` the single germline entity each of those
# calls resolves to, empty where it resolves to none. Returns the entity and the
# status, in that order.
#
# Only members that resolve to one entity vote. A member whose call spans two
# germlines is not evidence for either, and apportioning it between them would
# invent a precision the aligner explicitly declined to provide, so it abstains
# rather than being counted or broken up. Where the voters disagree the clone
# is reported as such and given no gene. Assigning the more numerous entity
# would hide the disagreement among the supplied member annotations.
#
# The vote is on the germline entity rather than the gene name, so a pair of IMGT
# names for one piece of sequence counts as agreement. That matters because a
# call naming such a pair is retained by clone inference on its shared token:
# comparing names would report those members as disagreeing about a germline they
# agree on, and inflate the one count this panel exists to expose.
#
# `calls` is still needed beside `tokens` to tell a clone whose members named
# no gene at all from one whose members named genes that resolve to no single
# entity. The pipeline's gates rule out both: common QC drops a record whose V
# or J call names no locus-compatible gene, and clone inference excludes a
# record whose call names several germline entities. Both are reported as
# categories that stay empty, so a nonzero count would be visible rather than
# silently absorbed.
.vote_clone_gene <- function(calls, tokens) {
  if (all(vj_gene_status(calls) == "missing")) {
    return(c("", "no_callable_member"))
  }
  entities <- unique(tokens[nzchar(tokens)])
  if (!length(entities)) return(c("", "members_ambiguous"))
  if (length(entities) > 1L) return(c("", "resolved_members_disagree"))
  c(entities, "assigned")
}

#' Give each donor-pooled clone one V and one J gene, or say why it has none.
#'
#' `membership` is a scope's whole clone-membership table, carrying
#' `clone_id`, `v_call`, `j_call` and the `v_gene` /
#' `j_gene` germline entities those calls resolve to. Each clone's
#' gene is voted over every
#' in-scope member, once across the scope, so the same clone reads the same way
#' in every library that carries it. SCOPer assigns clones within exact V/J
#' analytical entities. This summary checks the supplied member annotations.
#'
#' Returns one row per clone (`clone_id`, `v_gene`, `j_gene`, `v_vote_status`,
#' `j_vote_status`) and no rows for an empty membership table. A member with no
#' `clone_id` at all gets no row, because nothing can ever look it up: the
#' per-library aggregation drops unassigned members before joining these on.
#'
#' Calls that name several germline entities are already out of clone
#' inference, so every pipeline member names one entity. A nonzero
#' `resolved_members_disagree` count indicates inconsistent supplied membership.
summarise_clone_vj_call_sets <- function(membership) {
  clone_id <- as.character(membership$clone_id)
  # Clones come out in code-point order, fixed for the whole pipeline by
  # set_pipeline_collation().
  members <- split(
    seq_along(clone_id),
    factor(clone_id, levels = sort(unique(clone_id)))
  )
  vote <- function(calls, tokens) {
    voted <- vapply(
      members,
      function(rows) .vote_clone_gene(calls[rows], tokens[rows]),
      character(2), USE.NAMES = FALSE)
    list(gene = voted[1L, ], status = voted[2L, ])
  }
  v <- vote(as.character(membership$v_call),
            as.character(membership$v_gene))
  j <- vote(as.character(membership$j_call),
            as.character(membership$j_gene))
  data.frame(
    clone_id = as.character(names(members)),
    v_gene = v$gene,
    j_gene = j$gene,
    v_vote_status = v$status,
    j_vote_status = j$status,
    stringsAsFactors = FALSE
  )
}
