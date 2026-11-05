#!/usr/bin/env Rscript
# Pool one donor's unchanged common-QC IGH rows and prepared grouping across
# that donor's libraries in the scope manifest, read in sorted library_id order.
# Namespace sequence, backbone and representative IDs as <library_id>~<id> with
# make_pooled_sequence_id(), revalidate the pooled mapping and write both tables
# in pooled sequence_id order. Grouping eligibility stays in the mapping; this
# script does not discard ineligible source rows or add their support again.

source(file.path(dirname(sub("^--file=", "",
  grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1])),
  "..", "lib", "runtime.R"))
load_pipeline_modules(c(
  "table_io", "assembly_backbone", "pooled_records"))
stage_log <- start_stage_log("pool_scope_donor_igh")
suppressPackageStartupMessages(library(optparse))
opts <- parse_script_arguments(list(
  make_option("--out_file", type = "character", help = "Path to write the donor's pooled common-QC IGH table to. Required."),
  make_option("--grouping_file", type = "character", help = "Path to write the pooled prepared grouping to. Required."),
  make_option("--reannotated_files", type = "character", help = "Comma-separated common-QC AIRR files, one per donor library in sorted library_id order. Required."),
  make_option("--grouping_files", type = "character", help = "Comma-separated IGH grouping files in the same order. Required."),
  make_option("--manifest_file", type = "character", help = "Scope manifest. Required."),
  make_option("--scope_id", type = "character", help = "Scope identifier. Required."),
  make_option("--donor_id", type = "character", help = "Donor identifier. Required.")
))
manifest <- read_pipeline_table(opts$manifest_file,
  c("scope_id", "donor_id", "library_id", "tissue_id"))
if (any(manifest$scope_id != opts$scope_id)) stop("Manifest belongs to another scope")
libraries <- sort(manifest$library_id[manifest$donor_id == opts$donor_id])
files <- split_path_list(opts$reannotated_files)
grouping_files <- split_path_list(opts$grouping_files)
if (!length(libraries) || anyDuplicated(libraries) ||
    length(files) != length(libraries) || length(grouping_files) != length(libraries) ||
    anyDuplicated(files) || anyDuplicated(grouping_files)) {
  stop("Pool inputs must cover every manifest library exactly once")
}
frames <- mappings <- list()
columns <- character()
for (i in seq_along(libraries)) {
  source <- read_pipeline_table(files[i], c("library_id", "sequence_id", "locus", "consensus_count"))
  if (any(source$library_id != libraries[i])) stop("Pool input belongs to another library")
  columns <- union(columns, names(source))
  source <- source[source$locus == "IGH", , drop = FALSE]
  mapping <- validate_grouping_mapping(source, read_pipeline_table(grouping_files[i]))
  source$sequence_id <- make_pooled_sequence_id(source$library_id, source$sequence_id)
  for (field in c("sequence_id", "assembly_backbone_id", "backbone_representative_sequence_id",
                  "inference_representative_sequence_id")) {
    present <- !is_airr_null(mapping[[field]])
    mapping[[field]][present] <- make_pooled_sequence_id(mapping$library_id[present], mapping[[field]][present])
  }
  frames[[i]] <- source
  mappings[[i]] <- mapping
}
pooled <- bind_frames_on_column_union(frames, columns)
pooled <- pooled[order(pooled$sequence_id), , drop = FALSE]
mapping <- do.call(rbind, mappings)
mapping <- mapping[match(pooled$sequence_id, mapping$sequence_id), , drop = FALSE]
mapping <- validate_grouping_mapping(pooled, mapping)
write_pipeline_table(pooled, opts$out_file)
write_pipeline_table(mapping, opts$grouping_file, GROUPING_COLUMNS)
log_info("Pooled %d common-QC IGH records from %d libraries", nrow(pooled), length(libraries))
finish_stage_log(stage_log)
