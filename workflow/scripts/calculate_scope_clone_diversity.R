#!/usr/bin/env Rscript
# Calculate library IGH diversity, dominance, curves, abundance and V/J usage
# from donor/scope memberships. aggregate_clone_support_per_library() sums
# transferred read support for each clone. build_repertoire_diversity() in
# diversity.R applies params.diversity.min_consensus_count to Shannon entropy,
# Hill q = 1 and q = 2, dominance and target-depth estimates. Observed richness,
# support and analytical rarefaction curves remain available below that minimum.

suppressPackageStartupMessages({
  library(optparse)
  library(dplyr)
})

source(file.path(dirname(sub("^--file=", "",
  grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1])),
  "..", "lib", "runtime.R"))
load_pipeline_modules(c(
  "table_io", "diversity", "gene_calls", "gene_usage"))
stage_log <- start_stage_log("calculate_scope_clone_diversity")

script_options <- parse_script_arguments(list(
  make_option("--manifest_file", type = "character",
              help = "Scope manifest naming each library's donor and tissue. Required."),
  make_option("--membership_files", type = "character",
              help = "Comma-separated donor clone membership tables for the scope. Required."),
  make_option("--scope_id", type = "character",
              help = "Scope the memberships were inferred within. Required."),
  make_option("--reference_metadata_file", type = "character",
              help = "Reference allele metadata naming each composite V/J entity. Required."),
  make_option("--rarefaction_targets", type = "character",
              help = "IGH support depth as an IGH=depth pair, or empty for no target. Required."),
  make_option("--minimum_consensus_count", type = "character",
              help = "Minimum IGH support for diversity values as an IGH=count pair, or empty for zero. Required."),
  make_option("--curves_file", type = "character",
              help = "Path to write the rarefaction curves to. Required."),
  make_option("--clone_diversity_file", type = "character",
              help = "Path to write one IGH diversity row per in-scope library to, without scope context. Required."),
  make_option("--group_abundance_file", type = "character",
              help = "Path to write the per-clone support vectors to. Required."),
  make_option("--gene_usage_file", type = "character",
              help = "Path to write clone V/J usage to. Required."),
  make_option("--vj_pair_usage_file", type = "character",
              help = "Path to write clone V and J pair usage to. Required.")
))

# The usage schema both arms of gene_usage.tsv share before
# collect_scope_gene_usage adds scope, tissue and donor IDs: the library and
# locus, the grouping unit and gene type, then the columns count_gene_usage()
# returns for the gene each clone voted for.
# calculate_library_clonotype_diversity.R writes the clonotype arm in the same
# columns.
IGH_GENE_USAGE_COLUMNS <- c(
  "library_id", "locus", "group_unit", "gene_type", GENE_USAGE_COUNT_COLUMNS
)

# The pair usage schema, which calculate_library_clonotype_diversity.R writes in
# the same columns before collect_scope_vj_pair_usage adds scope, tissue and
# donor IDs.
IGH_VJ_PAIR_USAGE_COLUMNS <- c(
  "library_id", "locus", "group_unit", VJ_PAIR_USAGE_COUNT_COLUMNS
)

manifest <- read_pipeline_table(
  script_options$manifest_file,
  required_columns = c("library_id", "donor_id", "tissue_id"),
  context = "Scope manifest")
# The outputs list libraries in library order, which collect_scope_group_abundance
# needs to merge these support vectors with the per-library non-IGH tables.
manifest <- manifest[order(manifest$library_id), , drop = FALSE]
membership_files <- split_path_list(script_options$membership_files)
scope_id <- script_options$scope_id
entity_reference <- read_vj_entity_reference(script_options$reference_metadata_file)

membership <- bind_rows(lapply(membership_files, read_pipeline_table))
require_columns(membership,
                c("library_id", "sequence_id", "clone_id",
                  "consensus_count", "v_call", "j_call",
                  "v_gene", "j_gene", "junction"),
                "Clone membership")
log_info("Scope %s: %d membership rows across %d libraries",
         scope_id, nrow(membership), nrow(manifest))

clone_calls <- summarise_clone_vj_call_sets(membership)

targets <- parse_named_values(script_options$rarefaction_targets)
minimums <- parse_named_values(script_options$minimum_consensus_count)
diversity_rows <- list()
curve_rows <- list()
abundance_rows <- list()
usage_rows <- list()
pair_usage_rows <- list()
for (row in seq_len(nrow(manifest))) {
  # The diversity, abundance, curve and gene usage collectors add scope, tissue
  # and donor context later.
  identifiers <- list(library_id = as.character(manifest$library_id[row]))
  lineages <- aggregate_clone_support_per_library(membership,
                                                   identifiers$library_id,
                                                   clone_calls)

  metrics <- build_repertoire_diversity(
    lineages, identifiers$library_id, "IGH", "clone",
    locus_setting(targets, "IGH"), locus_setting(minimums, "IGH", 0))
  curve_rows[[identifiers$library_id]] <- metrics$curve
  diversity_rows[[identifiers$library_id]] <- metrics$detail
  # The same aggregation the row above was summarised from, released per clone
  # by support rank. Equally supported clones follow their smallest member
  # junction in this library, then clone_id. The scope identifiers are attached
  # by collect_scope_group_abundance, which is where the non-IGH half of the
  # table gets them too.
  abundance_rows[[identifiers$library_id]] <- build_group_abundance_rows(
    lineages, identifiers$library_id, "IGH", "clone", "clone_id",
    order_columns = c("min_junction", "clone_id"))

  if (!nrow(lineages)) next
  for (gene_type in c("V", "J")) {
    # Common QC and clone eligibility leave every member a call that resolves to
    # one entity, and clone inference groups within one V and J entity, so every
    # clone's members vote for one gene. Any other vote means the supplied
    # membership is inconsistent, and the usage would silently lose clones.
    vote_status <- lineages[[paste0(tolower(gene_type), "_vote_status")]]
    unassigned <- table(vote_status[vote_status != "assigned"])
    if (length(unassigned)) {
      stop(sprintf(paste("%s has IGH clones assigned no %s gene (%s); clone",
                         "membership is inconsistent with common QC and clone inference"),
                   identifiers$library_id, gene_type,
                   paste(names(unassigned), unassigned, sep = ": ", collapse = ", ")))
    }
    # One row per gene the library's clones were assigned.
    usage_rows[[paste(identifiers$library_id, gene_type)]] <- data.frame(
      library_id = identifiers$library_id, locus = "IGH", group_unit = "clone",
      gene_type = gene_type,
      count_gene_usage(lineages[[paste0(tolower(gene_type), "_gene")]],
                       lineages$consensus_count, entity_reference),
      stringsAsFactors = FALSE
    )
  }
  # The check above leaves every clone a V and a J entity, so each one has a pair.
  pair_usage_rows[[identifiers$library_id]] <- data.frame(
    library_id = identifiers$library_id, locus = "IGH", group_unit = "clone",
    count_vj_pair_usage(lineages$v_gene, lineages$j_gene, lineages$consensus_count,
                        entity_reference),
    stringsAsFactors = FALSE
  )
}

diversity_out <- bind_rows(diversity_rows)
abundance_out <- bind_rows(abundance_rows)
if (nrow(abundance_out) == 0L) {
  abundance_out <- make_empty_group_abundance_table()
}
usage_out <- bind_rows(usage_rows)
if (nrow(usage_out) == 0L) {
  usage_out <- make_empty_frame(IGH_GENE_USAGE_COLUMNS)
}
pair_usage_out <- bind_rows(pair_usage_rows)
if (nrow(pair_usage_out) == 0L) {
  pair_usage_out <- make_empty_frame(IGH_VJ_PAIR_USAGE_COLUMNS)
}

write_pipeline_table(bind_rows(curve_rows),
                     script_options$curves_file,
                     columns = RAREFACTION_CURVE_COLUMNS)
write_pipeline_table(diversity_out,
                     script_options$clone_diversity_file,
                     columns = DIVERSITY_DETAIL_COLUMNS)
write_pipeline_table(abundance_out,
                     script_options$group_abundance_file,
                     columns = CLONE_ABUNDANCE_COLS)
write_pipeline_table(usage_out,
                     script_options$gene_usage_file,
                     columns = IGH_GENE_USAGE_COLUMNS)
write_pipeline_table(pair_usage_out,
                     script_options$vj_pair_usage_file,
                     columns = IGH_VJ_PAIR_USAGE_COLUMNS)
log_info("Wrote %d IGH clone-diversity rows, %d clone-abundance rows, %d IGH gene-usage rows and %d IGH V and J pair rows",
         nrow(diversity_out), nrow(abundance_out), nrow(usage_out), nrow(pair_usage_out))

finish_stage_log(stage_log)
