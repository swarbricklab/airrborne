# Run selection while preserving clone populations, unavailable regions and Sigma.
# Run: python test/run_tests.py --component selection --level integration

load_test_modules(c(
  "table_io", "gene_calls", "isotype_resolution", "baseline_selection",
  "assembly_backbone", "qc_audit"))

work_dir <- tempfile("selection-regression-")
dir.create(work_dir)

pooled_file <- file.path(work_dir, "pooled.tsv")
membership_file <- file.path(work_dir, "membership.tsv")
manifest_file <- file.path(work_dir, "manifest.tsv")
feature_file <- file.path(work_dir, "features.tsv")
output_file <- file.path(work_dir, "selection.tsv")

replace_base <- function(sequence, position, base) {
  chars <- strsplit(sequence, "", fixed = TRUE)[[1]]
  chars[position] <- base
  paste(chars, collapse = "")
}

germline <- paste(rep("A", shazam::IMGT_V@seqLength), collapse = "")
regions <- as.character(shazam::IMGT_V@boundaries)
stopifnot(regions[80] == "cdr", regions[6] == "fwr")

# AAA -> ACA at the middle base of a CDR codon is a replacement mutation.
cdr_only <- replace_base(germline, 80L, "C")

# AAA -> AAG at the third base of an FWR codon is silent. All CDR positions are
# dot-padded to reproduce the unobserved-sequence case.
dot_cdr_chars <- strsplit(germline, "", fixed = TRUE)[[1]]
dot_cdr_chars[regions == "cdr"] <- "."
dot_cdr_chars[6] <- "G"
dot_cdr <- paste(dot_cdr_chars, collapse = "")

pooled <- data.frame(
  sequence_id = c(
    "cdr_only~cdr_only_0", "cdr_only~cdr_only_1", "dot_cdr~dot_cdr_0",
    "no_lineage~no_lineage_0", "no_dmask~no_dmask_0"
  ),
  library_id = c("cdr_only", "cdr_only", "dot_cdr", "no_lineage", "no_dmask"),
  sequence_alignment = c(cdr_only, cdr_only, dot_cdr, germline, germline),
  germline_alignment_d_mask = c(germline, germline, germline, germline, ""),
  v_call = "IGHV1-1*01",
  j_call = "IGHJ4*01",
  consensus_count = 1L,
  stringsAsFactors = FALSE
)
membership <- data.frame(
  library_id = pooled$library_id[c(1, 2, 3, 5)],
  sequence_id = pooled$sequence_id[c(1, 2, 3, 5)],
  # One donor-pooled lineage may be observed in several libraries. Reusing the
  # identity here proves that collapseClones still constructs one effective
  # sequence per library because library_id is an explicit grouping field.
  clone_id = "scope_donor_1",
  consensus_count = c(5L, 1L, 3L, 2L),
  stringsAsFactors = FALSE
)
manifest <- data.frame(
  library_id = c("cdr_only", "dot_cdr", "empty", "no_dmask", "no_lineage"),
  donor_id = "donor",
  tissue_id = paste0("tissue_", seq_len(5)),
  stringsAsFactors = FALSE
)
features <- data.frame(
  library_id = pooled$library_id,
  sequence_id = c(
    "cdr_only_0", "cdr_only_1", "dot_cdr_0", "no_lineage_0", "no_dmask_0"
  ),
  # One lineage in cdr_only is represented independently in IgG and IgA.
  qualified_c_call = c("IGHG1*01", "IGHA1*01", "IGHM*01", "IGHG1*01", "IGHG1*01"),
  shm_frequency = NA_real_,
  stringsAsFactors = FALSE
)

write.table(pooled, pooled_file, sep = "\t", row.names = FALSE, quote = FALSE)
write.table(
  membership, membership_file, sep = "\t", row.names = FALSE, quote = FALSE
)
write.table(manifest, manifest_file, sep = "\t", row.names = FALSE, quote = FALSE)
write.table(features, feature_file, sep = "\t", row.names = FALSE, quote = FALSE,
            na = "")

run <- run_pipeline_script("calculate_donor_selection.R", c(
  paste0("--qc_audit_file=", audit_file <- tempfile(fileext = ".tsv")),
  paste0("--pooled_file=", pooled_file),
  paste0("--membership_file=", membership_file),
  paste0("--manifest_file=", manifest_file),
  paste0("--out_file=", output_file),
  "--scope_id=scope", "--donor_id=donor",
  "--threads=1", "--collapse_min_freq=0.6",
  paste0("--isotype_call_files=", feature_file)
))
expect_script_ok(run)

result <- read.delim(output_file, stringsAsFactors = FALSE, check.names = FALSE)
audit <- read_qc_audit(audit_file)
# Membership supplies transferred support (5 + 1 for cdr_only), while the
# pooled representatives each carry one original read. The no_lineage record is
# outside membership; no_dmask enters the audit but fails selection eligibility.
# Both cdr_only representatives collapse to one effective clone in all_IGH.
test_that("selection records distinct input populations even when effective clone counts agree", {
  counts <- audit[audit$isotype_group == "all_IGH" & audit$stage == "selection_input", ]
  counts <- counts[match(c("cdr_only", "dot_cdr", "empty", "no_dmask", "no_lineage"), counts$library_id), ]
  expect_equal(counts$n_units_in, c(2, 1, 0, 1, 0))
  expect_equal(counts$n_units_out, c(2, 1, 0, 0, 0))
  expect_equal(counts$consensus_count_in, c(6, 3, 0, 2, 0))
  expect_equal(counts$consensus_count_out, c(6, 3, 0, 0, 0))
  expect_true(all(counts$input_unit == "inference_representative"))
  expect_true(all(counts$output_unit == "selection_input_record"))
  source_counts <- audit[audit$isotype_group == "all_IGH" &
                          audit$stage == "selection_input_source_support", ]
  source_counts <- source_counts[match(counts$library_id, source_counts$library_id), ]
  expect_equal(source_counts$n_units_in, counts$n_units_in)
  expect_equal(source_counts$n_units_out, counts$n_units_out)
  expect_equal(source_counts$consensus_count_in, c(2, 1, 0, 1, 0))
  expect_equal(source_counts$consensus_count_out, c(2, 1, 0, 0, 0))
  estimates <- result[result$isotype_group == "all_IGH", ]
  expect_equal(estimates$n_clone_consensus_sequences[match(c("cdr_only", "dot_cdr"), estimates$library_id)], c(1, 1))
  expect_equal(estimates$sigma_cdr_status[estimates$library_id == "dot_cdr"], "no_usable_mutation_pdfs")
})

expected_cols <- c(
  "scope_id", "library_id", "donor_id", "tissue_id", "isotype_group", "n_clone_consensus_sequences",
  "n_clone_consensus_sequences_with_positions_cdr", "n_clone_consensus_sequences_used_cdr", "sigma_cdr_status",
  "sigma_cdr", "sigma_cdr_ci_low", "sigma_cdr_ci_high",
  "n_clone_consensus_sequences_with_positions_fwr", "n_clone_consensus_sequences_used_fwr", "sigma_fwr_status",
  "sigma_fwr", "sigma_fwr_ci_low", "sigma_fwr_ci_high"
)
test_that("every library gets one row per isotype group", {
  expect_identical(names(result), expected_cols)
  expect_equal(nrow(result), 20L)
  expect_identical(unique(result$isotype_group), c("all_IGH", "IgM", "IgG", "IgA"))
  expect_true(!anyDuplicated(result[c("library_id", "isotype_group")]))
})

# Masking the CDR the record never observed is what leaves a representative with
# no CDR mutation PDF, so the two masked representatives are the only two rows
# that hold a lineage yet report no CDR estimate: dot_cdr's all-IGH row and its
# IgM row.
test_that("the all-IGH and IgM representatives of dot_cdr are masked", {
  masked_rows <- result[
    result$sigma_cdr_status == "no_usable_mutation_pdfs", , drop = FALSE
  ]
  expect_equal(nrow(masked_rows), 2L)
  expect_identical(unique(masked_rows$library_id), "dot_cdr")
  expect_setequal(masked_rows$isotype_group, c("all_IGH", "IgM"))
})

# Three of the five pooled records carry both a lineage and a D-masked germline,
# and each is the only member of its library x isotype x lineage group, so each
# becomes one representative in exactly one resolved isotype group. Summing
# n_clone_consensus_sequences over those groups therefore counts the records selection admitted.
test_that("only the records with a lineage and a germline are admitted", {
  expect_equal(sum(result$n_clone_consensus_sequences[result$isotype_group != "all_IGH"]), 3L)
})

cdr_row <- result[
  result$library_id == "cdr_only" & result$isotype_group == "all_IGH", , drop = FALSE
]
test_that("a mutated CDR gives a CDR estimate and no FWR estimate", {
  expect_equal(cdr_row$n_clone_consensus_sequences, 1L)
  expect_equal(cdr_row$n_clone_consensus_sequences_used_cdr, 1L)
  expect_equal(cdr_row$sigma_cdr_status, "estimated")
  expect_true(is.finite(cdr_row$sigma_cdr))
  expect_true(cdr_row$sigma_cdr_ci_low <= cdr_row$sigma_cdr)
  expect_true(cdr_row$sigma_cdr <= cdr_row$sigma_cdr_ci_high)
  # The value the whole stage produces, so a change in how calcBaseline is called
  # cannot pass unnoticed behind a finite number and an ordered interval. The sign
  # is positive because this record's one mutation replaces an amino acid in the
  # CDR; that is a property of this fixture, not of B cells, and a real repertoire
  # is not expected to select in any particular direction.
  expect_gt(cdr_row$sigma_cdr, 0)
  expect_equal(cdr_row$sigma_cdr, 2.2661557844, tolerance = 1e-6)
  expect_equal(cdr_row$n_clone_consensus_sequences_used_fwr, 0L)
  expect_equal(cdr_row$sigma_fwr_status, "no_usable_mutation_pdfs")
  expect_true(all(is.na(cdr_row[c(
    "sigma_fwr", "sigma_fwr_ci_low", "sigma_fwr_ci_high"
  )])))
  expect_equal(cdr_row$n_clone_consensus_sequences_with_positions_cdr, 1L)
})

dot_row <- result[
  result$library_id == "dot_cdr" & result$isotype_group == "all_IGH", , drop = FALSE
]
test_that("a dot-padded CDR gives an FWR estimate and no CDR estimate", {
  expect_equal(dot_row$n_clone_consensus_sequences, 1L)
  expect_equal(dot_row$n_clone_consensus_sequences_used_cdr, 0L)
  expect_equal(dot_row$sigma_cdr_status, "no_usable_mutation_pdfs")
  expect_true(all(is.na(dot_row[c(
    "sigma_cdr", "sigma_cdr_ci_low", "sigma_cdr_ci_high"
  )])))
  expect_equal(dot_row$n_clone_consensus_sequences_used_fwr, 1L)
  expect_equal(dot_row$sigma_fwr_status, "estimated")
  expect_true(is.finite(dot_row$sigma_fwr))
  expect_true(dot_row$sigma_fwr_ci_low <= dot_row$sigma_fwr)
  expect_true(dot_row$sigma_fwr <= dot_row$sigma_fwr_ci_high)
  # Negative here because this record's one framework mutation is silent, which
  # is again a property of the fixture rather than an expectation about a real
  # repertoire.
  expect_lt(dot_row$sigma_fwr, 0)
  expect_equal(dot_row$sigma_fwr, -4.1732686130, tolerance = 1e-6)
  # The dot-padded CDR has no callable position, so it cannot give a CDR density.
  expect_equal(dot_row$n_clone_consensus_sequences_with_positions_cdr, 0L)
  expect_equal(dot_row$n_clone_consensus_sequences_with_positions_fwr, 1L)
})

excluded_rows <- result[
  result$library_id %in% c("empty", "no_dmask", "no_lineage") &
    result$isotype_group == "all_IGH", , drop = FALSE
]
test_that("a library with no eligible lineage reports why", {
  expect_true(all(excluded_rows$n_clone_consensus_sequences == 0L))
  expect_true(all(excluded_rows$n_clone_consensus_sequences_used_cdr == 0L))
  expect_true(all(excluded_rows$n_clone_consensus_sequences_used_fwr == 0L))
  expect_true(all(excluded_rows$sigma_cdr_status == "no_eligible_clones"))
  expect_true(all(excluded_rows$sigma_fwr_status == "no_eligible_clones"))
  expect_true(all(is.na(excluded_rows[c(
    "sigma_cdr", "sigma_cdr_ci_low", "sigma_cdr_ci_high",
    "sigma_fwr", "sigma_fwr_ci_low", "sigma_fwr_ci_high"
  )])))
})

cdr_isotypes <- result[
  result$library_id == "cdr_only" & result$isotype_group %in% c("IgG", "IgA"),
  , drop = FALSE
]
cdr_igm <- result[
  result$library_id == "cdr_only" & result$isotype_group == "IgM", , drop = FALSE
]
dot_igm <- result[
  result$library_id == "dot_cdr" & result$isotype_group == "IgM", , drop = FALSE
]
test_that("one lineage contributes to every isotype it contains", {
  expect_equal(nrow(cdr_isotypes), 2L)
  expect_true(all(cdr_isotypes$n_clone_consensus_sequences == 1L))
  expect_true(all(cdr_isotypes$n_clone_consensus_sequences_used_cdr == 1L))
  expect_true(all(cdr_isotypes$sigma_cdr_status == "estimated"))
  expect_equal(cdr_igm$n_clone_consensus_sequences, 0L)
  expect_equal(cdr_igm$sigma_cdr_status, "no_eligible_clones")
  expect_equal(dot_igm$n_clone_consensus_sequences, 1L)
  expect_equal(dot_igm$n_clone_consensus_sequences_used_fwr, 1L)
  expect_equal(dot_igm$sigma_fwr_status, "estimated")
})

# The lineage membership table is what admits a sequence to selection, so the
# alternative CDR3s TRUST4 reports for an assembled contig cannot each cast a
# vote in the lineage consensus: they are reduced to one row before membership is
# written and therefore never reach this stage. Adding such rows to the pooled
# substrate without adding them to membership must leave every value unchanged.
sibling_pooled <- rbind(pooled, transform(
  pooled[1, , drop = FALSE],
  sequence_id = "cdr_only~cdr_only_2"
))
sibling_features <- rbind(features, transform(
  features[1, , drop = FALSE],
  sequence_id = "cdr_only_2"
))
sibling_pooled_file <- file.path(work_dir, "pooled_siblings.tsv")
sibling_feature_file <- file.path(work_dir, "features_siblings.tsv")
sibling_output_file <- file.path(work_dir, "selection_siblings.tsv")
write.table(sibling_pooled, sibling_pooled_file, sep = "\t", row.names = FALSE,
            quote = FALSE)
write.table(sibling_features, sibling_feature_file, sep = "\t",
            row.names = FALSE, quote = FALSE, na = "")
sibling_run <- run_pipeline_script("calculate_donor_selection.R", c(
  paste0("--qc_audit_file=", sibling_audit_file <- tempfile(fileext = ".tsv")),
  paste0("--pooled_file=", sibling_pooled_file),
  paste0("--membership_file=", membership_file),
  paste0("--manifest_file=", manifest_file),
  paste0("--out_file=", sibling_output_file),
  "--scope_id=scope", "--donor_id=donor",
  "--threads=1", "--collapse_min_freq=0.6",
  paste0("--isotype_call_files=", sibling_feature_file)
))
expect_script_ok(sibling_run)

test_that("pooled rows outside membership change nothing", {
  expect_identical(
    readLines(sibling_output_file, warn = FALSE),
    readLines(output_file, warn = FALSE)
  )
  expect_identical(readLines(sibling_audit_file), readLines(audit_file))
})

# A donor can reach this stage with nothing to estimate at all: an empty pool, or
# one whose records never joined a lineage. Its libraries still owe the scope a
# row each per isotype group, because an absent row and a row saying nothing was
# estimated read very differently downstream.
empty_pooled_file <- file.path(work_dir, "pooled_empty.tsv")
empty_membership_file <- file.path(work_dir, "membership_empty.tsv")
empty_feature_file <- file.path(work_dir, "features_empty.tsv")
empty_manifest_file <- file.path(work_dir, "manifest_empty.tsv")
empty_output_file <- file.path(work_dir, "selection_empty.tsv")

empty_manifest <- data.frame(
  library_id = c("empty_one", "empty_two"),
  donor_id = "donor",
  tissue_id = c("tissue_one", "tissue_two"),
  stringsAsFactors = FALSE
)
write.table(pooled[0, , drop = FALSE], empty_pooled_file, sep = "\t",
            row.names = FALSE, quote = FALSE)
write.table(membership[0, , drop = FALSE], empty_membership_file, sep = "\t",
            row.names = FALSE, quote = FALSE)
write.table(features[0, , drop = FALSE], empty_feature_file, sep = "\t",
            row.names = FALSE, quote = FALSE, na = "")
write.table(empty_manifest, empty_manifest_file, sep = "\t", row.names = FALSE,
            quote = FALSE)

empty_run <- run_pipeline_script("calculate_donor_selection.R", c(
  paste0("--qc_audit_file=", tempfile(fileext = ".tsv")),
  paste0("--pooled_file=", empty_pooled_file),
  paste0("--membership_file=", empty_membership_file),
  paste0("--manifest_file=", empty_manifest_file),
  paste0("--out_file=", empty_output_file),
  "--scope_id=scope", "--donor_id=donor",
  "--threads=1", "--collapse_min_freq=0.6",
  paste0("--isotype_call_files=", empty_feature_file)
))
expect_script_ok(empty_run)

empty_result <- read.delim(empty_output_file, stringsAsFactors = FALSE,
                           check.names = FALSE)
test_that("the donor stage rejects invalid frequency for populated and empty donors", {
  for (empty in c(FALSE, TRUE)) {
    for (frequency in c("0", "-0.1", "1.01", "NaN", "Inf")) {
      rejected_output <- tempfile(tmpdir = work_dir, fileext = ".tsv")
      rejected_audit <- tempfile(tmpdir = work_dir, fileext = ".tsv")
      rejected <- run_pipeline_script("calculate_donor_selection.R", c(
        paste0("--qc_audit_file=", rejected_audit),
        paste0("--pooled_file=", if (empty) empty_pooled_file else pooled_file),
        paste0("--membership_file=", if (empty) empty_membership_file else membership_file),
        paste0("--manifest_file=", if (empty) empty_manifest_file else manifest_file),
        paste0("--out_file=", rejected_output),
        "--scope_id=scope", "--donor_id=donor", "--threads=1",
        paste0("--collapse_min_freq=", frequency),
        paste0("--isotype_call_files=", if (empty) empty_feature_file else feature_file)))
      expect_false(identical(rejected$status, 0L))
      expect_true(any(grepl("collapse_min_freq must be one finite number in (0, 1]",
                           rejected$output, fixed = TRUE)))
      expect_false(file.exists(rejected_output))
      expect_false(file.exists(rejected_audit))
    }
  }
})

test_that("a donor with no eligible record still gets the full grid", {
  expect_identical(names(empty_result), SELECTION_COLUMNS)
  expect_equal(nrow(empty_result),
               nrow(empty_manifest) * length(ANALYSIS_IGH_ISOTYPE_GROUPS))
  expect_setequal(empty_result$library_id, empty_manifest$library_id)
  expect_identical(unique(empty_result$isotype_group), ANALYSIS_IGH_ISOTYPE_GROUPS)
  expect_true(all(empty_result$n_clone_consensus_sequences == 0L))
  expect_true(all(empty_result$sigma_cdr_status == "no_eligible_clones"))
  expect_true(all(empty_result$sigma_fwr_status == "no_eligible_clones"))
  expect_true(all(empty_result$n_clone_consensus_sequences_with_positions_cdr == 0))
  expect_true(all(empty_result$n_clone_consensus_sequences_with_positions_fwr == 0))
  expect_true(all(is.na(empty_result[c(
    "sigma_cdr", "sigma_cdr_ci_low", "sigma_cdr_ci_high",
    "sigma_fwr", "sigma_fwr_ci_low", "sigma_fwr_ci_high"
  )])))
})
