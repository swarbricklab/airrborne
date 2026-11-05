#!/usr/bin/env Rscript
# Calculate pairwise sharing, donor summaries and shared-clone support from
# supplied IGH memberships and non-IGH keys. overlap_metrics.R owns the pair
# metrics and donor occupancy summaries. Pair support minima flag shallow pairs
# without masking their values. Empty endpoints remain in the pair table with
# unavailable comparisons.

suppressPackageStartupMessages(library(optparse))

source(file.path(dirname(sub("^--file=", "",
  grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1])),
  "..", "lib", "runtime.R"))
load_pipeline_modules(c(
  "table_io", "gene_calls", "diversity", "gene_usage", "observed_representatives", "overlap_metrics",
  "pooled_records"))
stage_log <- start_stage_log("calculate_donor_overlap")

NONIGH_LOCI <- c("IGK", "IGL", "TRA", "TRB", "TRD", "TRG")

script_options <- parse_script_arguments(c(
  list(
    make_option("--membership_file", type = "character",
                help = "IGH clone membership for one donor in one scope. Required."),
    make_option("--records_file", type = "character",
                help = "IGH lineage records for the same donor and scope. Required."),
    make_option("--scope_id", type = "character",
                help = "Scope the clones were inferred within. Required."),
    make_option("--donor_id", type = "character",
                help = "Donor whose libraries are compared. Required."),
    make_option("--library_ids", type = "character",
                help = "Comma-separated libraries of this donor in the scope. Required."),
    make_option("--tissue_ids", type = "character",
                help = "Comma-separated tissues, aligned with --library_ids. Required."),
    make_option("--loci", type = "character",
                help = "Comma-separated loci to compare. Required."),
    make_option("--minimum_consensus_count", type = "character",
                help = paste("Minimum pair support per locus as LOCUS=count pairs;",
                             "an omitted locus uses zero. Required.")),
    make_option("--donor_overlap_pairs_file", type = "character",
                help = "Path to write the pairwise overlap rows to. Required."),
    make_option("--donor_overlap_file", type = "character",
                help = "Path to write the donor and locus summaries to. Required."),
    make_option("--shared_groups_file", type = "character",
                help = "Path to write the shared-clone rows to. Required."),
    make_option("--reference_metadata_file", type = "character",
                help = "Reference allele metadata naming composite V/J entities. Required."),
    make_option("--airr_files", type = "character",
                help = paste("Comma-separated common QC-pass AIRR tables, one per library in",
                             "--library_ids order. Required."))),
  lapply(NONIGH_LOCI, function(locus) {
    make_option(paste0("--", tolower(locus), "_files"), type = "character", default = "",
                help = sprintf(paste("Comma-separated %s clonotype tables, one per library",
                                     "in --library_ids order. Required when --loci names %s."),
                               locus, locus))
  })
))

scope_id <- script_options$scope_id
donor_id <- script_options$donor_id
library_ids <- split_path_list(script_options$library_ids)
tissue_ids <- split_path_list(script_options$tissue_ids)
loci <- split_path_list(script_options$loci)
entity_labels <- read_vj_entity_labels(script_options$reference_metadata_file)
pair_minima <- parse_named_values(script_options$minimum_consensus_count)

# The rule selects donors with at least two libraries in the scope. Check that
# library and tissue identifiers align so a tissue cannot be assigned to the
# wrong library.
if (length(library_ids) != length(tissue_ids) || anyDuplicated(library_ids)) {
  stop("Overlap library and tissue identifiers must be unique and aligned")
}
tissue_by_sample <- setNames(tissue_ids, library_ids)
airr_files <- split_path_list(script_options$airr_files)

# One library's common-QC sequences by library-qualified ID. A clonotype table
# names its representative record, and that record's sequence decides between
# equally supported observations in summarise_shared_clone_descriptors(). The
# file path names its library, as for the clonotype tables below. A library with
# no common-QC records has no sequences, and paste() would still build one ID.
read_library_sequences <- function(library_id, path) {
  if (!library_id %in% strsplit(path, "/", fixed = TRUE)[[1]]) {
    stop(sprintf("Common QC-pass AIRR input %s is not a %s file", path, library_id))
  }
  records <- read_pipeline_table(
    path, required_columns = c("sequence_id", "sequence"),
    context = sprintf("Common QC-pass AIRR table for %s", library_id))
  if (!nrow(records)) return(character(0))
  setNames(records$sequence, make_pooled_sequence_id(library_id, records$sequence_id))
}

membership <- read_pipeline_table(
  script_options$membership_file,
  required_columns = c("library_id", "clone_id", "consensus_count"),
  context = sprintf("IGH clone membership for scope %s donor %s",
                    scope_id, donor_id)
)
log_info("Scope %s donor %s: %d IGH clone membership rows over %d library(s)",
         scope_id, donor_id, nrow(membership), length(library_ids))

# Membership supplies clone support. Lineage records supply the junction and gene
# calls that describe shared clones, and the sequence that decides between tied
# observations.
igh_records <- read_pipeline_table(
  script_options$records_file,
  required_columns = c("library_id", "clone_id", "sequence_id",
                       "consensus_count", "junction", "junction_aa",
                       "v_gene", "j_gene", "sequence"),
  context = sprintf("IGH lineage records for scope %s donor %s",
                    scope_id, donor_id)
)


all_pairs <- list()
all_summaries <- list()
all_shared_groups <- list()
library_sequences <- list()
for (locus in loci) {
  paths <- NULL
  if (locus != "IGH") {
    supplied <- script_options[[paste0(tolower(locus), "_files")]]
    if (is.null(supplied) || !nzchar(supplied)) {
      stop(sprintf("--loci names %s but --%s_files was not supplied", locus, tolower(locus)))
    }
    paths <- split_path_list(supplied)
  }
  group_unit <- if (locus == "IGH") "clone" else "clonotype"

  repertoires <- list()
  descriptor_records <- list()
  for (index in seq_along(library_ids)) {
    library_id <- library_ids[index]
    repertoire <- if (locus == "IGH") {
      rows <- membership[as.character(membership$library_id) == library_id, , drop = FALSE]
      descriptor_records[[library_id]] <- as_descriptor_records(
        igh_records[as.character(igh_records$library_id) == library_id, , drop = FALSE],
        "clone_id", "sequence_id", entity_labels = entity_labels)
      build_sharing_units(
        rows$clone_id,
        as_read_support(rows$consensus_count, "consensus_count")
      )
    } else {
      # The group tables arrive as one positional list per locus, so this pairs
      # paths[index] with library_ids[index]. Every per-library output lives under
      # libraries/<library_id>/, which is what makes the pairing checkable: a
      # reordered or short input list would otherwise attribute one library's
      # receptors to another and the released table would look ordinary.
      if (!library_id %in% strsplit(paths[index], "/", fixed = TRUE)[[1]]) {
        stop(sprintf("%s overlap input %d is %s, which is not a %s file",
                     locus, index, paths[index], library_id))
      }
      groups <- read_pipeline_table(
        paths[index],
        required_columns = c("clonotype_id", "clonotype_key",
                             "consensus_count", "junction", "junction_aa",
                             "v_gene_key", "j_gene_key", "representative_sequence_id"),
        context = sprintf("%s rearrangement groups for %s", locus, library_id)
      )
      if (is.null(library_sequences[[library_id]])) {
        library_sequences[[library_id]] <- read_library_sequences(library_id, airr_files[index])
      }
      groups$representative_sequence <- unname(
        library_sequences[[library_id]][groups$representative_sequence_id])
      if (anyNA(groups$representative_sequence)) {
        stop(sprintf("%s clonotypes of %s name representatives absent from %s",
                     locus, library_id, airr_files[index]))
      }
      descriptor_records[[library_id]] <- as_descriptor_records(
        groups, "clonotype_key", "clonotype_id", "representative_sequence",
        "representative_sequence_id", c("v_gene_key", "j_gene_key"),
        entity_labels = entity_labels)
      build_sharing_units(
        groups$clonotype_key,
        as_read_support(groups$consensus_count, "consensus_count")
      )
    }
    if (nrow(repertoire)) {
      repertoires[[library_id]] <- repertoire
    }
  }

  if (length(repertoires) < 2L) {
    log_info("[%s] fewer than two libraries contributed; each pair is reported with an empty endpoint", locus)
  }
  pair_df <- build_pairwise_overlap_rows(repertoires, donor_id, locus, group_unit,
    tissue_by_sample, libraries = library_ids,
    minimum_consensus_count = locus_setting(pair_minima, locus, 0))
  if (nrow(pair_df)) {
    all_pairs[[locus]] <- pair_df
  }
  all_summaries[[locus]] <- summarise_donor_locus_occupancy(
    repertoires, pair_df, donor_id, locus, group_unit, length(library_ids))

  descriptors <- summarise_shared_clone_descriptors(
    do.call(rbind, unname(descriptor_records)))
  shared_df <- build_donor_locus_shared_clone_rows(
    repertoires, descriptors, donor_id, locus, group_unit, library_ids,
    tissue_by_sample)
  if (nrow(shared_df)) {
    all_shared_groups[[locus]] <- shared_df
  }
  log_info("[%s] %d clone(s) detected in two or more of the donor's libraries",
           locus, length(unique(shared_df$group_key)))
}

empty <- make_empty_overlap_tables()
pair_rows <- if (length(all_pairs)) do.call(rbind, all_pairs) else empty$donor_overlap_pairs
summary_rows <- if (length(all_summaries)) do.call(rbind, all_summaries) else empty$donor_overlap_summary
# The shared-clone rows run by locus in text order; within a locus,
# build_donor_locus_shared_clone_rows() has ordered them by clone support.
shared_clone_rows <- if (length(all_shared_groups)) {
  do.call(rbind, unname(all_shared_groups[sort(names(all_shared_groups))]))
} else empty$shared_groups

outputs <- script_options[c("donor_overlap_pairs_file", "donor_overlap_file",
                            "shared_groups_file")]
ensure_output_dirs(outputs)
write_pipeline_table(pair_rows, outputs$donor_overlap_pairs_file,
                     columns = OVERLAP_SCHEMAS$donor_overlap_pairs)
write_pipeline_table(summary_rows,
                     outputs$donor_overlap_file,
                     columns = OVERLAP_SCHEMAS$donor_overlap_summary)
write_pipeline_table(shared_clone_rows,
                     outputs$shared_groups_file,
                     columns = OVERLAP_SCHEMAS$shared_groups)
log_info("Wrote %d pairwise, %d donor/locus overlap and %d shared-clone rows",
         nrow(pair_rows), nrow(summary_rows), nrow(shared_clone_rows))

finish_stage_log(stage_log)
