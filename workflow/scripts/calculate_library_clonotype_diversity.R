#!/usr/bin/env Rscript
# Calculate library diversity, dominance, curves and support vectors from supplied
# non-IGH clonotypes. Each input row is one exact nucleotide rearrangement group.
# build_repertoire_diversity() in diversity.R applies params.diversity.min_consensus_count
# to estimates while retaining observed richness, support and curves.
# Every non-IGH locus receives a diversity row and a support vector, even when
# its group table is empty. params.rarefaction.targets supplies target depths,
# and only a locus with a target receives a curve.
# count_gene_usage() in gene_usage.R counts clonotypes and their read support for
# every non-IGH locus, over the clonotypes assigned a V or J entity, and
# count_vj_pair_usage() counts them for each observed V and J pair.

suppressPackageStartupMessages(library(optparse))

source(file.path(dirname(sub("^--file=", "",
  grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1])),
  "..", "lib", "runtime.R"))
load_pipeline_modules(c(
  "table_io", "diversity", "gene_calls", "gene_usage"))
stage_log <- start_stage_log("calculate_library_clonotype_diversity")

# Outputs list loci in this text order, the order collect_scope_group_abundance
# merges each library's loci in.
NONIGH_LOCI <- c("IGK", "IGL", "TRA", "TRB", "TRD", "TRG")

script_options <- parse_script_arguments(c(
  lapply(NONIGH_LOCI, function(locus) {
    make_option(paste0("--", tolower(locus), "_file"), type = "character",
                help = sprintf("%s clonotype table for one library. Required.", locus))
  }),
  list(
    make_option("--library_id", type = "character",
                help = "Library the clonotype tables belong to. Required."),
    make_option("--reference_metadata_file", type = "character",
                help = "Reference allele metadata naming each composite V/J entity. Required."),
    make_option("--rarefaction_targets", type = "character",
                help = paste("Positive support depth per locus as LOCUS=depth pairs;",
                             "an omitted locus has no target. Required.")),
    make_option("--minimum_consensus_count", type = "character",
                help = paste("Minimum support for diversity values per locus as",
                             "LOCUS=count pairs; an omitted locus uses zero. Required.")),
    make_option("--diversity_file", type = "character",
                help = "Path to write one diversity row per locus to. Required."),
    make_option("--curves_file", type = "character",
                help = "Path to write the rarefaction curves to. Required."),
    make_option("--group_abundance_file", type = "character",
                help = "Path to write the released support vectors to. Required."),
    make_option("--gene_usage_file", type = "character",
                help = "Path to write V/J usage to. Required."),
    make_option("--vj_pair_usage_file", type = "character",
                help = "Path to write V and J pair usage to. Required."))
))
locus_files <- setNames(
  lapply(NONIGH_LOCI, function(locus) script_options[[paste0(tolower(locus), "_file")]]),
  NONIGH_LOCI)
sample_name <- script_options$library_id
targets <- parse_named_values(script_options$rarefaction_targets)
minimums <- parse_named_values(script_options$minimum_consensus_count)
entity_reference <- read_vj_entity_reference(script_options$reference_metadata_file)

# The usage schema both arms of gene_usage.tsv share before
# collect_scope_gene_usage adds scope, tissue and donor IDs: the identifiers this
# stage prepends, then the columns count_gene_usage() returns.
# calculate_scope_clone_diversity.R writes the IGH clone arm in the same
# columns.
NONIGH_GENE_USAGE_COLUMNS <- c(
  "library_id", "locus", "group_unit", "gene_type", GENE_USAGE_COUNT_COLUMNS
)

# The pair usage schema, which calculate_scope_clone_diversity.R writes in the
# same columns before collect_scope_vj_pair_usage adds scope, tissue and donor IDs.
NONIGH_VJ_PAIR_USAGE_COLUMNS <- c(
  "library_id", "locus", "group_unit", VJ_PAIR_USAGE_COUNT_COLUMNS
)

diversity_rows <- list()
curve_rows <- list()
abundance_rows <- list()
usage_rows <- list()
pair_usage_rows <- list()
for (locus in names(locus_files)) {
  df <- read_pipeline_table(
    locus_files[[locus]],
    required_columns = c("clonotype_key", "n_members",
                         "consensus_count", "v_gene_key",
                         "j_gene_key"),
    context = sprintf("%s rearrangement groups", locus))
  log_info("[%s] loaded %d nucleotide rearrangement groups", locus, nrow(df))

  # The clone-size vector is released beside each diversity row, so the two
  # tables describe the same repertoires. Its rows follow support rank, then
  # clonotype_key. A locus without a rarefaction target has no curve rows.
  metrics <- build_repertoire_diversity(
    df, sample_name, locus, "clonotype", locus_setting(targets, locus),
    locus_setting(minimums, locus, 0))
  diversity_rows[[locus]] <- metrics$detail
  curve_rows[[locus]] <- metrics$curve
  abundance_rows[[locus]] <- build_group_abundance_rows(
    df, sample_name, locus, "clonotype", "clonotype_key")
  for (gene_type in c("V", "J")) {
    usage <- count_gene_usage(
      df[[paste0(tolower(gene_type), "_gene_key")]], df$consensus_count,
      entity_reference
    )
    if (!nrow(usage)) next
    usage_rows[[paste(locus, gene_type)]] <- cbind(
      library_id = sample_name, locus = locus, group_unit = "clonotype",
      gene_type = gene_type, usage
    )
  }
  pair_usage <- count_vj_pair_usage(df$v_gene_key, df$j_gene_key, df$consensus_count,
                                    entity_reference)
  if (nrow(pair_usage)) {
    pair_usage_rows[[locus]] <- cbind(
      library_id = sample_name, locus = locus, group_unit = "clonotype", pair_usage
    )
  }
}

diversity_out <- do.call(rbind, diversity_rows)
abundance_out <- do.call(rbind, abundance_rows)
usage_out <- if (length(usage_rows)) {
  do.call(rbind, usage_rows)
} else {
  make_empty_frame(NONIGH_GENE_USAGE_COLUMNS)
}
pair_usage_out <- if (length(pair_usage_rows)) {
  do.call(rbind, pair_usage_rows)
} else {
  make_empty_frame(NONIGH_VJ_PAIR_USAGE_COLUMNS)
}

write_pipeline_table(do.call(rbind, curve_rows),
                     script_options$curves_file,
                     columns = RAREFACTION_CURVE_COLUMNS)
write_pipeline_table(diversity_out, script_options$diversity_file,
                     columns = DIVERSITY_DETAIL_COLUMNS)
write_pipeline_table(abundance_out,
                     script_options$group_abundance_file,
                     columns = CLONE_ABUNDANCE_COLS)
write_pipeline_table(usage_out,
                     script_options$gene_usage_file,
                     columns = NONIGH_GENE_USAGE_COLUMNS)
write_pipeline_table(pair_usage_out,
                     script_options$vj_pair_usage_file,
                     columns = NONIGH_VJ_PAIR_USAGE_COLUMNS)
log_info("Wrote %d diversity rows, %d clone-abundance rows, %d gene-usage rows and %d V and J pair rows",
         nrow(diversity_out), nrow(abundance_out), nrow(usage_out), nrow(pair_usage_out))

finish_stage_log(stage_log)
