# Run donor fitting, scope threshold choice and IGH clone assignment.
# Run: python test/run_tests.py --component igh_clones --level integration

load_test_modules(c(
  "table_io", "gene_calls", "assembly_backbone", "donor_distance_threshold",
  "rearrangement_keys", "igh_clones", "igh_grouping", "qc_audit"))

prepare_test_grouping <- function(path) {
  frame <- read_pipeline_table(path)
  if (!"library_id" %in% names(frame)) frame$library_id <- rep("test_library", nrow(frame))
  # Give distinct sources different full-DNA flanks while retaining the junction
  # distances this calibration fixture tests. Fixed lengths preserve ranking ties.
  flanks <- vapply(seq_len(nrow(frame)), function(i)
    paste(ifelse(as.logical(intToBits(i)[1:16]), "C", "A"), collapse = ""), character(1))
  defaults <- list(sequence_alignment = frame$junction, sequence = paste0(frame$junction, flanks),
    germline_alignment_d_mask = frame$junction, junction_length = nchar(frame$junction),
    junction_aa = rep("fixture", nrow(frame)),
    v_gene = sub("\\*.*$", "", frame$v_call), j_gene = sub("\\*.*$", "", frame$j_call))
  for (segment in c("v", "j")) {
    if (!paste0(segment, "_gene") %in% names(frame) &&
        any(grepl("[,;|]", frame[[paste0(segment, "_call")]]))) {
      stop("Constructed multi-gene calls require explicit analytical entities")
    }
  }
  for (name in names(defaults)) {
    if (!name %in% names(frame)) frame[[name]] <- defaults[[name]]
  }
  write_pipeline_table(frame, path)
  mapping <- prepare_igh_grouping(frame, unique(frame$library_id), 0)
  grouping <- tempfile(fileext = ".tsv")
  write_pipeline_table(mapping, grouping, GROUPING_COLUMNS)
  grouping
}

# The combined entry point persists representatives only. Exercise its propagation
# helper separately against the assigned labels to retain the source-row assertions.
assigned_source <- function(input, members) {
  frame <- read_pipeline_table(input)
  mapping <- prepare_igh_grouping(frame, unique(frame$library_id), 0)
  db <- frame[mapping$grouping_eligible, , drop = FALSE]
  eligible <- mapping[mapping$grouping_eligible, , drop = FALSE]
  representative <- match(eligible$inference_representative_sequence_id, db$sequence_id)
  db$clone_id <- propagate_clone_id_to_siblings(db, representative, read_pipeline_table(members))
  db
}

write_scope_threshold <- function(scope_id, threshold) {
  path <- tempfile(fileext = ".tsv")
  write.table(data.frame(
    scope_id = scope_id, clone_threshold = threshold,
    threshold_source = "configured", n_donors = 0L, n_donors_fitted = 0L
  ), path, sep = "\t", row.names = FALSE, quote = FALSE)
  path
}

write_diagnostics <- function(scope_id, donor_id) {
  path <- tempfile(fileext = ".tsv")
  breaks <- seq(0, 1, by = 0.01)
  write.table(data.frame(
    scope_id = scope_id, donor_id = donor_id,
    distance_bin_start = head(breaks, -1L),
    distance_bin_end = tail(breaks, -1L),
    distance_bin_midpoint = head(breaks, -1L) + 0.005,
    n_distances_in_bin = 0L, observed_density = 0,
    near_component_density = NA_real_, far_component_density = NA_real_
  ), path, sep = "\t", row.names = FALSE, quote = FALSE, na = "NA")
  path
}

write_donor <- function(scope_id, donor_id, threshold = NA_real_,
                        n_distances = 0L, fit_status = "too_few_distances") {
  path <- tempfile(fileext = ".tsv")
  row <- data.frame(
    scope_id = scope_id, donor_id = donor_id, n_libraries = 1L,
    n_members = n_distances, n_distances = n_distances,
    n_distances_fitted = n_distances, fit_status = fit_status, threshold = threshold
  )
  for (column in names(DONOR_FIT_PARAMETERS)) row[[column]] <- NA_real_
  write.table(row, path, sep = "\t", row.names = FALSE, quote = FALSE, na = "NA")
  path
}

# Exact-entity DTN fitting retains a fixed 100-bin review surface even when
# the donor is too sparse to fit.
patterns <- data.frame(
  sequence_id = sprintf("seq%03d", seq_len(121L)),
  junction = c(rep("AAA", 30L), rep("AAT", 30L), rep(c("AAAA", "AAAT"), each = 15L),
               rep("AAAT", 30L), "AAN"),
  v_call = c(rep(c(
    "IGHV1-1*01", "IGHV2-1*01", "IGHV1-1*01,IGHV2-1*01", "IGHV1-1*01"
  ), each = 30L), "IGHV1-1*01"),
  j_call = "IGHJ4*01", locus = "IGH", consensus_count = 1L,
  v_gene = c(rep(c("IGHV1-1", "IGHV2-1", "VSEQ_complete_a", "IGHV1-1"),
                 each = 30L), "IGHV1-1"),
  j_gene = "IGHJ4",
  # The tie kind the reannotation stage wrote against the germline reference.
  # The multi-gene rows here stand for a tie the reference reduces to one entity,
  # which lineage inference keeps. Two junctions within that composite supply
  # distances without joining any single-gene entity.
  v_gene_status = c(rep(c("resolved", "resolved", "assigned_group",
                            "resolved"), each = 30L), "resolved"),
  j_gene_status = "resolved",
  stringsAsFactors = FALSE
)
donor_input <- tempfile(fileext = ".tsv")
donor_output <- tempfile(fileext = ".tsv")
donor_diagnostics <- tempfile(fileext = ".tsv")
write.table(patterns, donor_input, sep = "\t", row.names = FALSE, quote = FALSE)
donor_run <- run_pipeline_script("fit_donor_threshold.R", c(
  paste0("--grouping_file=", prepare_test_grouping(donor_input)),
  paste0("--input_file=", donor_input),
  paste0("--output_file=", donor_output),
  paste0("--diagnostic_file=", donor_diagnostics),
  "--scope_id=test_scope", "--donor_id=test_donor",
  "--minimum_n=500",
    "--fit_edge=0.9", "--fit_specificity=0.995", "--fit_subsample=10000"
))
expect_script_ok(donor_run)

test_that("a donor too sparse to fit still leaves a full review surface", {
  donor_result <- read.delim(donor_output, stringsAsFactors = FALSE)
  diagnostics <- read.delim(donor_diagnostics, stringsAsFactors = FALSE)
  expect_equal(nrow(donor_result), 1L)
  expect_true(donor_result$n_distances > 0L)
  expect_identical(donor_result$fit_status, "too_few_distances")
  expect_identical(donor_result$n_libraries, 1L)
  expect_true(donor_result$n_members >= donor_result$n_distances)
  expect_equal(nrow(diagnostics), 100L)
  expect_identical(names(diagnostics), c(
    "scope_id", "donor_id", "distance_bin_start", "distance_bin_end",
    "distance_bin_midpoint", "n_distances_in_bin", "observed_density",
    "near_component_density", "far_component_density"
  ))
})

# The same separated mixture exercises the file handoff into fitting and assignment.
fit_fixture <- separated_distance_fixture()
fit_fixture$v_gene_status <- "resolved"
fit_fixture$j_gene_status <- "resolved"
n_pairs <- nrow(fit_fixture) / 2L
pair_index <- rep(seq_len(n_pairs), each = 2L)
fit_input <- tempfile(fileext = ".tsv")
fit_output <- tempfile(fileext = ".tsv")
fit_diagnostics <- tempfile(fileext = ".tsv")
write.table(fit_fixture, fit_input, sep = "\t", row.names = FALSE, quote = FALSE)
fit_run <- run_pipeline_script("fit_donor_threshold.R", c(
  paste0("--grouping_file=", prepare_test_grouping(fit_input)),
  paste0("--input_file=", fit_input),
  paste0("--output_file=", fit_output),
  paste0("--diagnostic_file=", fit_diagnostics),
  "--scope_id=fit_scope", "--donor_id=fit_donor",
  "--minimum_n=500",
    "--fit_edge=0.9", "--fit_specificity=0.995", "--fit_subsample=10000"
))
expect_script_ok(fit_run)
fit_result <- read.delim(fit_output, stringsAsFactors = FALSE)

test_that("a separated mixture is fitted", {
  expect_equal(fit_result$n_distances, 600L)
  expect_true(is.finite(fit_result$threshold))
  expect_identical(fit_result$fit_status, "fitted")
  # Every row of the fixture is its own retained representative in one library.
  expect_identical(c(fit_result$n_libraries, fit_result$n_members), c(1L, 600L))
})

fit_output_repeat <- tempfile(fileext = ".tsv")
fit_diagnostics_repeat <- tempfile(fileext = ".tsv")
repeat_run <- run_pipeline_script("fit_donor_threshold.R", c(
  paste0("--grouping_file=", prepare_test_grouping(fit_input)),
  paste0("--input_file=", fit_input),
  paste0("--output_file=", fit_output_repeat),
  paste0("--diagnostic_file=", fit_diagnostics_repeat),
  "--scope_id=fit_scope", "--donor_id=fit_donor",
  "--minimum_n=500",
    "--fit_edge=0.9", "--fit_specificity=0.995", "--fit_subsample=10000"
))
expect_script_ok(repeat_run)

test_that("the same input gives the same fit and the same review surface", {
  fit_result_repeat <- read.delim(fit_output_repeat, stringsAsFactors = FALSE)
  expect_identical(fit_result_repeat$fit_status, "fitted")
  expect_identical(fit_result_repeat$threshold, fit_result$threshold)
  expect_identical(
    readLines(fit_diagnostics_repeat, warn = FALSE),
    readLines(fit_diagnostics, warn = FALSE)
  )
})

# Resolve the fitted donor through the same scope handoff used by assignment.
# Its retained calibration must describe the fit that supplied the cutoff.
fit_scope_output <- tempfile(fileext = ".tsv")
fit_calibration_output <- tempfile(fileext = ".tsv")
fit_scope_diagnostics <- tempfile(fileext = ".tsv")
fit_scope_run <- run_pipeline_script("resolve_scope_clone_threshold.R", c(
  paste0("--output_file=", fit_scope_output),
  paste0("--donor_output_file=", fit_calibration_output),
  paste0("--diagnostics_output_file=", fit_scope_diagnostics),
  "--scope_id=fit_scope", "--configured_threshold=auto",
  paste0("--donor_files=", fit_output),
  paste0("--diagnostic_files=", fit_diagnostics)
))
expect_script_ok(fit_scope_run)

test_that("the successful fit supplies the scope cutoff and retained calibration", {
  threshold <- read.delim(fit_scope_output, stringsAsFactors = FALSE)
  calibration <- read.delim(fit_calibration_output, stringsAsFactors = FALSE)
  expect_identical(threshold$threshold_source, "automatic_fit")
  expect_equal(threshold$clone_threshold, fit_result$threshold)
  expect_equal(c(threshold$n_donors, threshold$n_donors_fitted), c(1L, 1L))
  expect_identical(names(calibration), c(
    "scope_id", "donor_id", "n_libraries", "n_members", "n_distances", "n_distances_fitted",
    "fit_status", "fitted_threshold", "near_weight", "near_shape", "near_scale",
    "far_weight", "far_mean", "far_sd", "absolute_log_likelihood", "fit_sensitivity"))
  expect_identical(calibration$n_distances_fitted, 600L)
  expect_identical(calibration$donor_id, "fit_donor")
  expect_equal(c(calibration$n_libraries, calibration$n_members, calibration$n_distances),
               c(1L, 600L, 600L))
  expect_identical(calibration$fit_status, "fitted")
  expect_equal(calibration$fitted_threshold, fit_result$threshold)
  expect_equal(read.delim(fit_scope_diagnostics), read.delim(fit_diagnostics))
})

# The written parameters keep enough digits to rebuild the written component
# densities. Values below the smallest normal double hold too few bits to compare.
test_that("the written mixture reproduces the written component densities", {
  fit <- read.delim(fit_calibration_output, stringsAsFactors = FALSE)
  bins <- read.delim(fit_scope_diagnostics, stringsAsFactors = FALSE)
  agrees <- function(rebuilt, written) {
    all(abs(rebuilt - written) <= pmax(1e-9 * pmax(abs(rebuilt), abs(written)),
                                       .Machine$double.xmin))
  }
  expect_true(agrees(fit$near_weight * dgamma(bins$distance_bin_midpoint, shape = fit$near_shape,
                                              scale = fit$near_scale),
                     bins$near_component_density))
  expect_true(agrees(fit$far_weight * dnorm(bins$distance_bin_midpoint, mean = fit$far_mean,
                                            sd = fit$far_sd),
                     bins$far_component_density))
})

test_that("a donor without a fit leaves every fit column empty", {
  calibration_file <- tempfile(fileext = ".tsv")
  run <- run_pipeline_script("resolve_scope_clone_threshold.R", c(
    paste0("--output_file=", tempfile(fileext = ".tsv")),
    paste0("--donor_output_file=", calibration_file),
    paste0("--diagnostics_output_file=", tempfile(fileext = ".tsv")),
    "--scope_id=test_scope", "--configured_threshold=auto",
    paste0("--donor_files=", donor_output),
    paste0("--diagnostic_files=", donor_diagnostics)
  ))
  expect_script_ok(run)
  row <- read_pipeline_table(calibration_file)
  expect_identical(row$fit_status, "too_few_distances")
  expect_true(all(unlist(row[c("n_distances_fitted", "fitted_threshold",
                               names(DONOR_FIT_PARAMETERS))]) == ""))
})

# The fitted cutoff is fed into assignment and the partition it produces is
# checked. Every other assignment case here supplies a hand-written threshold,
# which cannot notice the fit moving. Each pair carries its own V gene, so a pair
# can only ever merge with itself: the near pairs, 1 to 4 mismatches apart, must
# each become one lineage, and the far pairs, 14 to 18 apart, must each stay two.
partition_fixture <- fit_fixture
partition_fixture$sequence_alignment <- partition_fixture$junction
partition_input <- tempfile(fileext = ".tsv")
partition_output <- tempfile(fileext = ".tsv")
write.table(partition_fixture, partition_input, sep = "\t", row.names = FALSE,
            quote = FALSE)
partition_run <- run_pipeline_script("assign_scope_donor_clones.R", c(
  paste0("--grouping_file=", prepare_test_grouping(partition_input)),
  paste0("--input_file=", partition_input),
  paste0("--threshold_file=", fit_scope_output),
  paste0("--membership_file=", partition_output),
  paste0("--records_file=", partition_records <- tempfile(fileext = ".tsv")),
  paste0("--qc_audit_file=", partition_audit <- tempfile(fileext = ".tsv")),
  "--donor_id=fit_donor",
  "--scope_id=fit_scope", "--max_n=0"
))
expect_script_ok(partition_run)

test_that("the fitted cutoff joins the near pairs and separates the far ones", {
  partition <- assigned_source(partition_input, partition_output)
  expect_equal(nrow(partition), 2L * n_pairs)
  clones_per_pair <- tapply(partition$clone_id, pair_index,
                            function(ids) length(unique(ids)))
  near <- seq_len(n_pairs / 2L)
  far <- setdiff(seq_len(n_pairs), near)
  expect_true(all(clones_per_pair[near] == 1L))
  expect_true(all(clones_per_pair[far] == 2L))
  # 150 joined pairs and 150 split ones.
  expect_equal(length(unique(partition$clone_id)), 450L)
})

# Hand-written class tokens stand for distinct complete reference sequences.
# This fixture tests the inference-tool boundary; annotation tests establish
# whether a named-allele call can receive such a token. An identical junction
# cannot merge single genes, duplicated-location genes or overlapping composites.
ambiguous <- data.frame(
  sequence_id = paste0("ambiguous~seq", seq_len(8), "_0"),
  locus = "IGH",
  junction = c(rep("AAACCCGGG", 3), "AAACCCGGA", rep("AAACCCGGG", 3), "AAACCCGGA"),
  v_call = c(
    "IGHV3-23*01", "IGHV3-23D*01", "IGHV3-23*01,IGHV3-23D*01",
    "IGHV3-23D*01,IGHV3-23*01", "IGHV3-23D*02,IGHV4-59*01",
    rep("IGHV3-23*01,IGHV3-23D*01", 3)
  ),
  j_call = c(
    rep("IGHJ4*01", 5), "IGHJ4*01,IGHJ5*01", "IGHJ5*02,IGHJ6*01",
    "IGHJ5*01,IGHJ4*01"
  ),
  v_gene = c("IGHV3-23", "IGHV3-23D", "VSEQ_complete_a", "VSEQ_complete_a",
             "VSEQ_complete_b", rep("VSEQ_complete_a", 3)),
  j_gene = c(rep("IGHJ4", 5), "JSEQ_complete_a", "JSEQ_complete_b", "JSEQ_complete_a"),
  v_gene_status = c(rep("resolved", 2), rep("assigned_group", 6)),
  j_gene_status = c(rep("resolved", 5), rep("assigned_group", 3)),
  consensus_count = 1L, stringsAsFactors = FALSE
)
ambiguous_input <- tempfile(fileext = ".tsv")
ambiguous_output <- tempfile(fileext = ".tsv")
write.table(ambiguous, ambiguous_input, sep = "\t", row.names = FALSE,
            quote = FALSE)
ambiguous_run <- run_pipeline_script("assign_scope_donor_clones.R", c(
  paste0("--grouping_file=", prepare_test_grouping(ambiguous_input)),
  paste0("--input_file=", ambiguous_input),
  paste0("--threshold_file=", write_scope_threshold("ambiguous_scope", 0.2)),
  paste0("--membership_file=", ambiguous_output),
  paste0("--records_file=", ambiguous_records <- tempfile(fileext = ".tsv")),
  paste0("--qc_audit_file=", ambiguous_audit <- tempfile(fileext = ".tsv")),
  "--donor_id=test_donor",
  "--scope_id=ambiguous_scope", "--max_n=0"
))
expect_script_ok(ambiguous_run)

test_that("assignment joins only the same complete-reference V and J entities", {
  ambiguous_result <- assigned_source(ambiguous_input, ambiguous_output)
  ambiguous_result <- ambiguous_result[match(ambiguous$sequence_id, ambiguous_result$sequence_id), ]
  expected <- c(1, 2, 3, 3, 4, 5, 6, 5)
  expect_equal(nrow(ambiguous_result), 8L)
  expect_identical(outer(ambiguous_result$clone_id, ambiguous_result$clone_id, `==`),
                   outer(expected, expected, `==`))
  expect_identical(ambiguous_result$v_call, ambiguous$v_call)
  expect_identical(ambiguous_result$j_call, ambiguous$j_call)
})

test_that("distance fitting uses the same exact entity boundaries as assignment", {
  result <- fit_donor_distance_threshold(ambiguous, 500L, 0.9, 0.995, 10000)
  # SHazaM 1.3.2 rounds its nearest distances to four decimal places.
  expect_equal(result$distances, rep(round(1 / 9, 4), 4), tolerance = 1e-12)
  expect_identical(result$fit_status, "too_few_distances")
})

# The alternative CDR3s of one assembly must inherit their representative's
# lineage instead of acting as extra vertices in a single-link chain. At 0.2 the
# minor row is within reach of both outer rows (2/12 = 0.167) while the two
# representatives are not within reach of each other (4/12 = 0.333), so
# representative-first clustering yields two lineages rather than one chain.
bridge <- data.frame(
  sequence_id = c("bridge~assemble1_0", "bridge~assemble1_1", "bridge~assemble2_0"),
  library_id = "bridge", locus = "IGH",
  junction = c("AAAAAAAAAAAA", "CCAAAAAAAAAA", "CCCCAAAAAAAA"),
  junction_aa = c("KKKK", "PKKK", "PPKK"),
  v_call = "IGHV1-1*01", j_call = "IGHJ4*01",
  # The germline entity each call resolves to, written by the reannotation stage
  # and carried through pooling; the membership stage builds its match key from
  # these rather than re-reading the calls. The tie kinds come from the same
  # stage and are what the lineage gate reads.
  v_gene = "IGHV1-1", j_gene = "IGHJ4",
  v_gene_status = "resolved", j_gene_status = "resolved",
  consensus_count = c(10L, 1L, 10L),
  stringsAsFactors = FALSE
)
bridge$sequence_alignment <- bridge$junction
bridge$junction_length <- nchar(bridge$junction)
# These rows exercise which records cluster together, not how mutated they are,
# so each carries its own observed sequence as its germline and no mutation is
# implied. The membership table needs the column present because a lineage's
# members are handed to a tree in this form.
bridge$germline_alignment_d_mask <- bridge$sequence_alignment
# Distinct per row, so the record shard's contig can be traced to the row it came
# from rather than to any sibling of the same backbone.
bridge$sequence <- paste0("GGGG", bridge$junction, "TTTT")
bridge_input <- tempfile(fileext = ".tsv")
bridge_output <- tempfile(fileext = ".tsv")
write.table(bridge, bridge_input, sep = "\t", row.names = FALSE, quote = FALSE)
bridge_run <- run_pipeline_script("assign_scope_donor_clones.R", c(
  paste0("--grouping_file=", prepare_test_grouping(bridge_input)),
  paste0("--input_file=", bridge_input),
  paste0("--threshold_file=", write_scope_threshold("bridge_scope", 0.2)),
  paste0("--membership_file=", bridge_output),
  paste0("--records_file=", bridge_records <- tempfile(fileext = ".tsv")),
  paste0("--qc_audit_file=", bridge_audit <- tempfile(fileext = ".tsv")),
  "--donor_id=test_donor",
  "--scope_id=bridge_scope", "--max_n=0"
))
expect_script_ok(bridge_run)

test_that("alternative CDR3s inherit their representative's lineage", {
  bridge_result <- assigned_source(bridge_input, bridge_output)
  bridge_clone <- setNames(as.character(bridge_result$clone_id),
                           as.character(bridge_result$sequence_id))
  expect_equal(nrow(bridge_result), 3L)
  expect_equal(length(unique(bridge_clone)), 2L)
  expect_equal(bridge_clone[["bridge~assemble1_0"]],
               bridge_clone[["bridge~assemble1_1"]])
  expect_true(bridge_clone[["bridge~assemble1_0"]] !=
                bridge_clone[["bridge~assemble2_0"]])
  # Assignment propagates labels and must not rewrite read support; the
  # membership stage is where siblings are reduced to one row.
  expect_identical(as.integer(bridge_result$consensus_count), c(10L, 1L, 10L))
})

test_that("an ineligible dominant variant has no minor replacement in assignment", {
  blocked <- bridge
  blocked$v_call[1] <- "IGHV1-1*01,IGHV2-1*01"
  blocked$v_gene[1] <- ""
  blocked$v_gene_status[1] <- "ambiguous"
  input <- tempfile(fileext = ".tsv")
  members <- tempfile(fileext = ".tsv")
  write_pipeline_table(blocked, input)
  grouping <- prepare_test_grouping(input)
  mapping <- read_pipeline_table(grouping)
  expect_identical(as.logical(mapping$grouping_eligible), c(FALSE, FALSE, TRUE))
  expect_identical(mapping$grouping_status[2], "ineligible_backbone_representative")
  run <- run_pipeline_script("assign_scope_donor_clones.R", c(
    paste0("--input_file=", input), paste0("--grouping_file=", grouping),
    paste0("--threshold_file=", write_scope_threshold("blocked_scope", 0.2)),
    paste0("--membership_file=", members),
    paste0("--records_file=", tempfile(fileext = ".tsv")),
    paste0("--qc_audit_file=", tempfile(fileext = ".tsv")),
    "--donor_id=test_donor", "--scope_id=blocked_scope", "--max_n=0"))
  expect_script_ok(run)
  result <- read_pipeline_table(members)
  expect_identical(result$sequence_id, "bridge~assemble2_0")
  expect_identical(result$consensus_count, "10")
})

# Control: the same three junctions on three separate assemblies really do chain
# into one lineage at this threshold, so the two-lineage result above is the
# backbone unit at work rather than a property of the distances.
distinct_contigs <- bridge
distinct_contigs$sequence_id <- c("bridge~assemble1_0", "bridge~assemble2_0",
                                  "bridge~assemble3_0")
distinct_input <- tempfile(fileext = ".tsv")
distinct_output <- tempfile(fileext = ".tsv")
write.table(distinct_contigs, distinct_input, sep = "\t", row.names = FALSE,
            quote = FALSE)
distinct_run <- run_pipeline_script("assign_scope_donor_clones.R", c(
  paste0("--grouping_file=", prepare_test_grouping(distinct_input)),
  paste0("--input_file=", distinct_input),
  paste0("--threshold_file=", write_scope_threshold("bridge_scope", 0.2)),
  paste0("--membership_file=", distinct_output),
  paste0("--records_file=", distinct_records <- tempfile(fileext = ".tsv")),
  paste0("--qc_audit_file=", distinct_audit <- tempfile(fileext = ".tsv")),
  "--donor_id=test_donor",
  "--scope_id=bridge_scope", "--max_n=0"
))
expect_script_ok(distinct_run)

test_that("three separate assemblies do chain into one lineage", {
  distinct_result <- assigned_source(distinct_input, distinct_output)
  expect_equal(nrow(distinct_result), 3L)
  expect_equal(length(unique(distinct_result$clone_id)), 1L)
})

# A record whose V or J call names several germline entities takes no part in a
# lineage. The gate is applied before clustering, so the record is not there
# to be assigned one, and the rest of the pool is unaffected.
reducible <- bridge
reducible$v_gene_status <- c("resolved", "resolved", "ambiguous")
reducible_input <- tempfile(fileext = ".tsv")
reducible_output <- tempfile(fileext = ".tsv")
write.table(reducible, reducible_input, sep = "\t", row.names = FALSE,
            quote = FALSE)
reducible_run <- run_pipeline_script("assign_scope_donor_clones.R", c(
  paste0("--grouping_file=", prepare_test_grouping(reducible_input)),
  paste0("--input_file=", reducible_input),
  paste0("--threshold_file=", write_scope_threshold("bridge_scope", 0.2)),
  paste0("--membership_file=", reducible_output),
  paste0("--records_file=", reducible_records <- tempfile(fileext = ".tsv")),
  paste0("--qc_audit_file=", reducible_audit <- tempfile(fileext = ".tsv")),
  "--donor_id=test_donor",
  "--scope_id=bridge_scope", "--max_n=0"
))
expect_script_ok(reducible_run)

test_that("a record naming several germline entities joins no lineage", {
  reducible_result <- assigned_source(reducible_input, reducible_output)
  expect_identical(as.character(reducible_result$sequence_id),
                   c("bridge~assemble1_0", "bridge~assemble1_1"))
  expect_equal(length(unique(reducible_result$clone_id)), 1L)
})

# The combined entry point writes both projections from the assigned representatives.
membership_output <- bridge_output
records_output <- bridge_records
membership <- read.delim(membership_output, stringsAsFactors = FALSE)
expected_uid <- c("bridge~assemble1_0", "bridge~assemble2_0")

test_that("one backbone becomes one membership row carrying its siblings' support", {
  expect_equal(nrow(membership), 2L)
  # The row kept as the membership unit is the row that was clustered.
  expect_identical(sort(as.character(membership$sequence_id)), sort(expected_uid))
  expect_equal(sum(as.numeric(membership$consensus_count)),
               sum(as.numeric(bridge$consensus_count)))
  expect_true(membership$consensus_count[
    membership$sequence_id == "bridge~assemble1_0"] == 11)
  expect_equal(length(unique(membership$clone_id)), 2L)
})

# The record shard describes exactly the same receptors, and adds the contig and
# the library-independent match key that a cross-sample lookup needs.
test_that("the record shard describes the same receptors", {
  records <- read.delim(records_output, stringsAsFactors = FALSE)
  expect_identical(as.character(records$sequence_id),
                   as.character(membership$sequence_id))
  expect_identical(as.character(records$clone_id),
                   as.character(membership$clone_id))
  expect_identical(as.numeric(records$consensus_count),
                   as.numeric(membership$consensus_count))
  # The contig kept is the representative's, not a sibling's.
  expect_identical(as.character(records$sequence),
                   as.character(bridge$sequence[match(records$sequence_id,
                                                      bridge$sequence_id)]))
  # Every row of the fixture resolves to the same V and J entity, so the key is
  # the junction with those two.
  expect_identical(as.character(records$clonotype_key),
                   paste(c("AAAAAAAAAAAA", "CCCCAAAAAAAA"), "IGHV1-1", "IGHJ4",
                         sep = "~"))
})

# A V/J-discordant backbone is still one assembly, so it becomes one member: the
# rows describe one contig and the annotations disagree because the substituted
# CDR3 moved the boundary the aligner reads.
discordant <- bridge[1:2, ]
discordant$v_call <- c("IGHV1-1*01", "IGHV2-1*01")
discordant_input <- tempfile(fileext = ".tsv")
discordant_output <- tempfile(fileext = ".tsv")
discordant_membership <- tempfile(fileext = ".tsv")
write.table(discordant, discordant_input, sep = "\t", row.names = FALSE,
            quote = FALSE)
discordant_assign <- run_pipeline_script("assign_scope_donor_clones.R", c(
  paste0("--grouping_file=", prepare_test_grouping(discordant_input)),
  paste0("--input_file=", discordant_input),
  paste0("--threshold_file=", write_scope_threshold("bridge_scope", 0.2)),
  paste0("--membership_file=", discordant_output),
  paste0("--records_file=", discordant_records <- tempfile(fileext = ".tsv")),
  paste0("--qc_audit_file=", discordant_audit <- tempfile(fileext = ".tsv")),
  "--donor_id=test_donor",
  "--scope_id=bridge_scope", "--max_n=0"
))
expect_script_ok(discordant_assign)
discordant_membership <- discordant_output

test_that("a V/J-discordant backbone becomes one member holding both supports", {
  discordant_result <- read.delim(discordant_membership, stringsAsFactors = FALSE)
  expect_equal(nrow(discordant_result), 1L)
  expect_equal(as.numeric(discordant_result$consensus_count), 11)
})

# Calibration measures distances between receptors, so adding alternative CDR3s
# to every assembly must not add distances.
backbone_junctions <- vapply(seq_len(20L), function(k) {
  paste0(strrep("A", 20L - k), strrep("C", k))
}, character(1))
one_per_contig <- data.frame(
  sequence_id = sprintf("assemble%d_0", seq_len(20L)),
  junction = backbone_junctions, v_call = "IGHV1-1*01", j_call = "IGHJ4*01",
  locus = "IGH", consensus_count = 10L,
  v_gene_status = "resolved", j_gene_status = "resolved",
  stringsAsFactors = FALSE
)
with_siblings <- do.call(rbind, lapply(seq_len(20L), function(k) {
  data.frame(
    sequence_id = sprintf("assemble%d_%d", k, 0:2),
    junction = c(backbone_junctions[k],
                 sub("^.", "G", backbone_junctions[k]),
                 sub("^.", "T", backbone_junctions[k])),
    v_call = "IGHV1-1*01", j_call = "IGHJ4*01", locus = "IGH",
    consensus_count = c(10L, 1L, 1L),
    v_gene_status = "resolved", j_gene_status = "resolved",
    stringsAsFactors = FALSE
  )
}))
count_distances <- function(frame) {
  input <- tempfile(fileext = ".tsv")
  out <- tempfile(fileext = ".tsv")
  diag <- tempfile(fileext = ".tsv")
  write.table(frame, input, sep = "\t", row.names = FALSE, quote = FALSE)
  run <- run_pipeline_script("fit_donor_threshold.R", c(
    paste0("--grouping_file=", prepare_test_grouping(input)),
  paste0("--input_file=", input),
    paste0("--output_file=", out),
    paste0("--diagnostic_file=", diag),
    "--scope_id=bb_scope", "--donor_id=bb_donor",
    "--minimum_n=500",
    "--fit_edge=0.9", "--fit_specificity=0.995", "--fit_subsample=10000"
  ))
  expect_script_ok(run)
  read.delim(out, stringsAsFactors = FALSE)$n_distances
}

test_that("alternative CDR3s add no calibration distances", {
  expect_identical(count_distances(one_per_contig),
                   count_distances(with_siblings))
})

# Calibration and clustering share one gate, so marking a record's call
# ambiguous leaves calibration exactly where deleting the record would.
test_that("a record naming several germline entities is not calibrated on", {
  marked <- one_per_contig
  marked$j_gene_status[1:5] <- "ambiguous"
  expect_identical(count_distances(marked),
                   count_distances(one_per_contig[6:20, , drop = FALSE]))
})

# One fitted donor is sufficient in automatic mode. Manual mode uses its numeric
# cutoff, and automatic mode with no fits uses the fixed human-BCR fallback 0.1.
single_donor <- write_donor("single_scope", "d1", 0.17, 600L, "fitted")
single_diag <- write_diagnostics("single_scope", "d1")
scope_output <- tempfile(fileext = ".tsv")
calibration_output <- tempfile(fileext = ".tsv")
scope_diagnostics <- tempfile(fileext = ".tsv")
scope_run <- run_pipeline_script("resolve_scope_clone_threshold.R", c(
  paste0("--output_file=", scope_output),
  paste0("--donor_output_file=", calibration_output),
  paste0("--diagnostics_output_file=", scope_diagnostics),
  "--scope_id=single_scope", "--configured_threshold=auto",
  paste0("--donor_files=", single_donor),
  paste0("--diagnostic_files=", single_diag)
))
expect_script_ok(scope_run)

test_that("one fitted donor is sufficient in automatic mode", {
  scope_result <- read.delim(scope_output, stringsAsFactors = FALSE)
  expect_equal(scope_result$threshold_source, "automatic_fit")
  expect_true(isTRUE(all.equal(scope_result$clone_threshold, 0.17)))
  expect_equal(c(scope_result$n_donors, scope_result$n_donors_fitted), c(1L, 1L))
})

manual_output <- tempfile(fileext = ".tsv")
manual_run <- run_pipeline_script("resolve_scope_clone_threshold.R", c(
  paste0("--output_file=", manual_output),
  paste0("--donor_output_file=", tempfile(fileext = ".tsv")),
  paste0("--diagnostics_output_file=", tempfile(fileext = ".tsv")),
  "--scope_id=single_scope", "--configured_threshold=0.08",
  paste0("--donor_files=", single_donor),
  paste0("--diagnostic_files=", single_diag)
))
expect_script_ok(manual_run)

test_that("manual mode uses the configured cutoff", {
  manual_result <- read.delim(manual_output, stringsAsFactors = FALSE)
  expect_equal(manual_result$threshold_source, "configured")
  expect_true(isTRUE(all.equal(manual_result$clone_threshold, 0.08)))
})

sparse_donor <- write_donor("sparse_scope", "d1")
sparse_diag <- write_diagnostics("sparse_scope", "d1")
sparse_output <- tempfile(fileext = ".tsv")
sparse_run <- run_pipeline_script("resolve_scope_clone_threshold.R", c(
  paste0("--output_file=", sparse_output),
  paste0("--donor_output_file=", tempfile(fileext = ".tsv")),
  paste0("--diagnostics_output_file=", tempfile(fileext = ".tsv")),
  "--scope_id=sparse_scope", "--configured_threshold=auto",
  paste0("--donor_files=", sparse_donor),
  paste0("--diagnostic_files=", sparse_diag)
))
expect_script_ok(sparse_run)

test_that("automatic mode with no fits uses the human-BCR fallback", {
  sparse_result <- read.delim(sparse_output, stringsAsFactors = FALSE)
  expect_equal(sparse_result$threshold_source, "automatic_fallback")
  expect_true(isTRUE(all.equal(sparse_result$clone_threshold, 0.1)))
  expect_equal(c(sparse_result$n_donors, sparse_result$n_donors_fitted), c(1L, 0L))
})

test_that("a failed donor fit still assigns all eligible representatives using the scope fit", {
  fitted <- write_donor("mixed_scope", "fitted", fit_result$threshold, fit_result$n_distances, "fitted")
  failed <- write_donor("mixed_scope", "failed", NA_real_, read.delim(donor_output)$n_distances, "fit_failed")
  resolved <- tempfile(fileext = ".tsv")
  run <- run_pipeline_script("resolve_scope_clone_threshold.R", c(
    paste0("--output_file=", resolved),
    paste0("--donor_output_file=", tempfile(fileext = ".tsv")),
    paste0("--diagnostics_output_file=", tempfile(fileext = ".tsv")),
    "--scope_id=mixed_scope", "--configured_threshold=auto",
    paste0("--donor_files=", fitted, ",", failed),
    paste0("--diagnostic_files=", write_diagnostics("mixed_scope", "fitted"), ",",
           write_diagnostics("mixed_scope", "failed"))
  ))
  expect_script_ok(run)
  threshold <- read.delim(resolved)
  expect_equal(threshold$clone_threshold, fit_result$threshold)
  expect_equal(c(threshold$n_donors, threshold$n_donors_fitted), c(2L, 1L))
  grouping <- prepare_test_grouping(donor_input)
  mapping <- read_pipeline_table(grouping)
  members <- tempfile(fileext = ".tsv")
  audit <- tempfile(fileext = ".tsv")
  run <- run_pipeline_script("assign_scope_donor_clones.R", c(
    paste0("--input_file=", donor_input), paste0("--grouping_file=", grouping),
    paste0("--threshold_file=", resolved), "--scope_id=mixed_scope", "--donor_id=failed",
    "--max_n=0", paste0("--membership_file=", members),
    paste0("--records_file=", tempfile(fileext = ".tsv")), paste0("--qc_audit_file=", audit)
  ))
  expect_script_ok(run)
  output <- read_pipeline_table(members)
  expect_equal(nrow(output), sum(as.logical(mapping$is_inference_representative)))
  expect_gt(nrow(output), 0L)
  expect_setequal(output$sequence_id, mapping$sequence_id[as.logical(mapping$is_inference_representative)])
  counts <- read_qc_audit(audit)
  expect_equal(counts$n_units_out[counts$stage == "donor_clones"], length(unique(output$clone_id)))
})

# Distinct valid context tuples must stay distinct through the assignment stage.
test_that("assignment preserves partitions while names distinguish scope and donor boundaries", {
  outputs <- lapply(1:2, function(i) {
    scope <- c("study_A", "study")[i]
    donor <- c("B", "A_B")[i]
    library <- c("LIB_A", "LIB_B")[i]
    db <- bridge[1:2, ]
    db$library_id <- library
    db$sequence_id <- paste0(library, "~assemble1_", 0:1)
    input <- tempfile(fileext = ".tsv")
    members_file <- tempfile(fileext = ".tsv")
    records_file <- tempfile(fileext = ".tsv")
    write_pipeline_table(db, input)
    run <- run_pipeline_script("assign_scope_donor_clones.R", c(
      paste0("--input_file=", input),
      paste0("--grouping_file=", prepare_test_grouping(input)),
      paste0("--threshold_file=", write_scope_threshold(scope, 0.1)),
      paste0("--membership_file=", members_file),
      paste0("--records_file=", records_file),
      paste0("--qc_audit_file=", tempfile(fileext = ".tsv")),
      paste0("--scope_id=", scope), paste0("--donor_id=", donor), "--max_n=0"))
    expect_script_ok(run)
    members <- read_pipeline_table(members_file)
    records <- read_pipeline_table(records_file)
    expect_equal(nrow(members), 1L)
    expect_identical(members$sequence_id, paste0(library, "~assemble1_0"))
    expect_equal(as.numeric(members$consensus_count), sum(as.numeric(db$consensus_count)))
    expect_identical(records$clone_id, members$clone_id)
    expect_identical(records$sequence_id, members$sequence_id)
    members$clone_id
  })
  expect_identical(unlist(outputs), c("study_A~B~1", "study~A_B~1"))
})

test_that("combined assignment writes both empty projections for empty or ineligible pools", {
  source <- data.frame(library_id = "0460", sequence_id = "0460~assemble1_0",
    sequence = "ACGT", sequence_alignment = "ACGT", locus = "IGH", consensus_count = 5,
    junction = "", v_call = "IGHV1-1*01", j_call = "IGHJ4*01",
    v_gene = "IGHV1-1", j_gene = "IGHJ4", v_gene_status = "resolved",
    j_gene_status = "resolved", stringsAsFactors = FALSE)
  threshold <- tempfile(fileext = ".tsv")
  write_pipeline_table(data.frame(scope_id = "scope", clone_threshold = 0.1,
    threshold_source = "configured", n_donors = 0, n_donors_fitted = 0), threshold)
  for (frame in list(source[0, ], source)) {
    input <- tempfile(fileext = ".tsv")
    grouping <- tempfile(fileext = ".tsv")
    members <- tempfile(fileext = ".tsv")
    records <- tempfile(fileext = ".tsv")
    audit <- tempfile(fileext = ".tsv")
    write_pipeline_table(frame, input)
    write_pipeline_table(prepare_igh_grouping(frame, "0460", 0), grouping)
    run <- run_pipeline_script("assign_scope_donor_clones.R", c(
      paste0("--input_file=", input), paste0("--grouping_file=", grouping),
      paste0("--threshold_file=", threshold), "--scope_id=scope", "--donor_id=001",
      "--max_n=0", paste0("--membership_file=", members),
      paste0("--records_file=", records), paste0("--qc_audit_file=", audit)))
    expect_script_ok(run)
    expect_identical(names(read_pipeline_table(members)), CLONE_MEMBER_COLUMNS)
    expect_identical(names(read_pipeline_table(records)), CLONE_RECORD_COLUMNS)
    expect_equal(nrow(read_pipeline_table(members)), 0L)
    expect_equal(nrow(read_pipeline_table(records)), 0L)
    counts <- read_qc_audit(audit)
    expect_true(all(counts$n_units_out == 0 & counts$consensus_count_out == 0))
    expect_true(all(counts$donor_id == "001"))
  }
})
