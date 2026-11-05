# Summarise one donor's scope-specific IGH clones and regional SHM evidence.
# write_donor_igh_summaries() joins membership, feature and source tables and
# writes six sorted summary tables plus the SHM population audit. Its call to
# choose_shm_evidence() in shm_summaries.R selects one qualifying backbone source
# per member and region. It ranks qualifying sources by supported positions,
# original read support and sequence ID.
# Missing qualifying evidence leaves regional SHM unavailable. Selection has its
# own alignment eligibility in calculate_donor_selection.R.

source(file.path(dirname(sub("^--file=", "",
  grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1])),
  "..", "lib", "runtime.R"))
load_pipeline_modules(c(
  "table_io", "observed_representatives", "igh_summaries", "isotype_resolution", "isotype_metrics", "shm_summaries",
  "gene_calls", "igh_features", "igh_clones", "pooled_records", "qc_audit",
  "shm_regions", "assembly_backbone", "gene_usage"))
stage_log <- start_stage_log("calculate_donor_igh_summaries")

suppressPackageStartupMessages(library(optparse))
script_options <- parse_script_arguments(list(
  make_option("--qc_audit_file", type = "character", help = "Regional SHM population counts. Required."),
  make_option("--pooled_file", type = "character", help = "Original donor common-QC IGH records. Required."),
  make_option("--grouping_file", type = "character", help = "Prepared backbone and exact-sequence links. Required."),
  make_option("--out_clone_shm", type = "character",
              help = "Clone support and scored mutation states. Required."),
  make_option("--out_shm", type = "character",
              help = "Per-library SHM output for this donor. Required."),
  make_option("--out_clones", type = "character",
              help = "Clone output for this donor. Required."),
  make_option("--out_lineage_isotype", type = "character",
              help = "Library-local clone isotype output for this donor. Required."),
  make_option("--out_members", type = "character",
              help = "Enriched clone-member output for this donor. Required."),
  make_option("--out_clone_member_shm_coverage", type = "character",
              help = "Complete member-region coverage audit. Required."),
  make_option("--out_clone_member_shm", type = "character",
              help = "Member and V-region SHM rows for this donor. Required."),
  make_option("--scope_id", type = "character",
              help = "Analysis scope of these clone assignments. Required."),
  make_option("--donor_id", type = "character",
              help = "Donor whose clone tables are being summarised. Required."),
  make_option("--manifest_file", type = "character",
              help = "Scope manifest identifying the donor's libraries. Required."),
  make_option("--feature_files", type = "character",
              help = "Comma-separated IGH feature tables for the donor's libraries. Required."),
  make_option("--membership_file", type = "character",
              help = "Scope-specific clone membership for this donor. Required."),
  make_option("--records_file", type = "character",
              help = "Matching clone records with sequences and alignments. Required."),
  make_option("--min_subregion_fraction", type = "double",
              help = "Required fraction of read-supported sub-region positions. Required."),
  make_option("--min_germline_occupied_v", type = "integer",
              help = "Minimum germline-occupied whole-V positions. Required."),
  make_option("--reference_metadata_file", type = "character",
              help = "Reference allele metadata naming composite V/J entities. Required.")
))

write_donor_igh_summaries(
  input_paths = list(records = script_options$records_file,
                     pooled = script_options$pooled_file,
                     grouping = script_options$grouping_file,
                     membership = script_options$membership_file,
                     features = split_path_list(script_options$feature_files),
                     manifest = script_options$manifest_file),
  output_paths = list(clone_shm = script_options$out_clone_shm, shm = script_options$out_shm,
                      clones = script_options$out_clones,
                      clone_isotype = script_options$out_lineage_isotype,
                      members = script_options$out_members,
                      clone_member_shm = script_options$out_clone_member_shm,
                      clone_member_shm_coverage = script_options$out_clone_member_shm_coverage,
                      qc_audit = script_options$qc_audit_file),
  scope_id = script_options$scope_id,
  donor_id = script_options$donor_id,
  min_subregion_fraction = script_options$min_subregion_fraction,
  min_germline_occupied_v = script_options$min_germline_occupied_v,
  region_levels = V_REGION_LEVELS,
  feature_columns = igh_feature_columns(V_REGION_LEVELS),
  entity_labels = read_vj_entity_labels(script_options$reference_metadata_file)
)
finish_stage_log(stage_log)
