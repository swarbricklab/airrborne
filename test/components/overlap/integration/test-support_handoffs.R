# Run pairwise sharing readers with exact endpoint support and file/library pairing.
# Run: python test/run_tests.py --component overlap --level integration

load_test_modules(c("table_io"))

test_that("overlap readers reject malformed counts and mispaired library files", {
  directory <- tempfile("overlap-support-")
  dir.create(directory)
  on.exit(unlink(directory, recursive = TRUE), add = TRUE)
  libraries <- c("0460", "NA", "empty")
  input <- list(membership = file.path(directory, "members.tsv"),
                records = file.path(directory, "records.tsv"),
                trb = file.path(directory, "libraries", libraries, "TRB.tsv"),
                airr = file.path(directory, "libraries", libraries, "airr.tsv"))
  products <- c("donor_overlap_pairs", "donor_overlap_summary", "shared_groups")
  output <- setNames(as.list(file.path(directory, paste0(products, ".tsv"))), products)
  params <- list(scope_id = "scope", donor_id = "donor", library_ids = libraries,
    tissue_ids = c("T1", "T2", "T3"), loci = c("IGH", "TRB"),
    minimum_consensus_count = list())
  members <- data.frame(library_id = c("0460", "0460", "NA"),
    clone_id = c("shared", "private", "shared"), consensus_count = c(9, 1, 5))
  records <- members
  records$sequence_id <- c("a", "b", "c")
  records$junction <- records$sequence <- c("ATG", "AAA", "ATG")
  records$junction_aa <- c("M", "K", "M")
  records$v_gene <- "VSEQ_x"
  records$j_gene <- "IGHJ4"
  groups <- data.frame(clonotype_id = c("a", "b"), clonotype_key = c("shared", "private"),
    consensus_count = c(9, 1), junction = c("ATG", "AAA"), junction_aa = c("M", "K"),
    v_gene_key = "TRBV1", j_gene_key = "TRBJ1-1", representative_sequence_id = c("0460~t1", "0460~t2"))
  reference <- file.path(directory, "reference_metadata.tsv")
  writeLines(c("allele\tidentity_id\tidentity_label",
               "IGHV1-2*01\tVSEQ_x\tIGHV1-2,IGHV1-2D"), reference)
  write_pipeline_table(members, input$membership)
  write_pipeline_table(records, input$records)
  for (index in seq_along(libraries)) {
    frame <- if (index == 1L) groups else if (index == 2L) groups[1, ] else groups[FALSE, ]
    if (index == 2L) {
      frame$consensus_count <- 5
      frame$representative_sequence_id <- "NA~t1"
    }
    dir.create(dirname(input$trb[index]), recursive = TRUE)
    write_pipeline_table(frame, input$trb[index])
    # The empty library has no common-QC records, so its AIRR table is a header.
    airr <- data.frame(sequence_id = c("t1", "t2"), sequence = c("ATG", "AAA"))
    write_pipeline_table(if (index == 3L) airr[FALSE, ] else airr, input$airr[index])
  }
  execute <- function(ok = TRUE) run_stage_with_options("calculate_donor_overlap.R",
    list(membership_file = input$membership, records_file = input$records,
         scope_id = params$scope_id, donor_id = params$donor_id,
         library_ids = params$library_ids, tissue_ids = params$tissue_ids,
         loci = params$loci, minimum_consensus_count = params$minimum_consensus_count,
         trb_files = input$trb, airr_files = input$airr,
         donor_overlap_pairs_file = output$donor_overlap_pairs,
         donor_overlap_file = output$donor_overlap_summary,
         shared_groups_file = output$shared_groups, reference_metadata_file = reference),
    expect_success = ok)
  execute()
  pairs <- read_pipeline_table(output$donor_overlap_pairs)
  expect_equal(nrow(pairs), 6L) # Three unordered library pairs at each of two loci.
  for (locus in params$loci) {
    pair <- pairs[pairs$locus == locus & pairs$library_id_2 == "NA", ]
    expect_identical(pair$library_id_1, "0460")
    expect_equal(as.numeric(pair$grouped_consensus_count_1), 10)
    expect_equal(as.numeric(pair$grouped_consensus_count_2), 5)
    expect_equal(as.numeric(pair$shared_support_fraction_1), 0.9)
    expect_equal(as.numeric(pair$shared_support_fraction_2), 1)
  }
  expect_true(all(pairs$pair_status[pairs$library_id_2 == "empty"] == "empty_second"))
  expect_identical(unique(pairs$group_unit[pairs$locus == "IGH"]), "clone")
  expect_identical(unique(pairs$group_unit[pairs$locus == "TRB"]), "clonotype")
  summary <- read_pipeline_table(output$donor_overlap_summary)
  expect_identical(summary$n_libraries_in_donor, c("3", "3"))
  expect_identical(summary$n_libraries_nonempty, c("2", "2"))
  expect_identical(summary$group_unit, c("clone", "clonotype"))
  expect_identical(summary$n_pairs_total, c("3", "3"))
  expect_identical(summary$n_pairs_available, c("1", "1"))
  # Both loci describe a shared clone by its entity, and a composite reads as
  # its member genes.
  shared <- read_pipeline_table(output$shared_groups)
  expect_false("representative_sequence" %in% names(shared))
  expect_identical(unique(shared$representative_junction[shared$locus == "TRB"]), "ATG")
  expect_identical(unique(shared$v_gene_key[shared$locus == "IGH"]), "VSEQ_x")
  expect_identical(unique(shared$v_gene_names[shared$locus == "IGH"]), "IGHV1-2,IGHV1-2D")
  expect_identical(unique(shared$v_gene_names[shared$locus == "TRB"]), "TRBV1")
  # The IGH representative is the best-supported lineage record; the TRB one is
  # the best-supported clonotype's representative.
  expect_identical(unique(shared$representative_sequence_id[shared$locus == "IGH"]), "a")
  expect_identical(unique(shared$representative_sequence_id[shared$locus == "TRB"]), "0460~t1")
  valid_output <- readLines(output$donor_overlap_pairs)
  for (branch in c("IGH", "TRB")) {
    invalid <- if (branch == "IGH") members else groups
    invalid$consensus_count[2] <- NA_real_
    path <- if (branch == "IGH") input$membership else input$trb[1]
    write_pipeline_table(invalid, path)
    run <- execute(FALSE)
    expect_gt(run$status, 0L)
    expect_match(paste(run$output, collapse = "\n"), "consensus_count.*invalid")
    expect_identical(readLines(output$donor_overlap_pairs), valid_output)
    write_pipeline_table(if (branch == "IGH") members else groups, path)
  }
  paths <- input$trb
  input$trb <- paths[c(2, 1, 3)]
  run <- execute(FALSE)
  expect_gt(run$status, 0L)
  expect_match(paste(run$output, collapse = "\n"), "not a 0460 file")
  input$trb <- paths[-1]
  expect_gt(execute(FALSE)$status, 0L)
})

# A clonotype table names its representative, and the common-QC AIRR table holds
# the sequence that decides between equally supported observations.
test_that("a tied clonotype takes the representative with the smaller full sequence", {
  directory <- tempfile("overlap-tie-")
  dir.create(directory)
  on.exit(unlink(directory, recursive = TRUE), add = TRUE)
  libraries <- c("A", "B")
  group_files <- file.path(directory, "libraries", libraries, "TRB.tsv")
  airr_files <- file.path(directory, "libraries", libraries, "airr.tsv")
  for (index in seq_along(libraries)) {
    dir.create(dirname(group_files[index]), recursive = TRUE)
    write_pipeline_table(data.frame(
      clonotype_id = paste0(libraries[index], "_TRB_ct00001"), clonotype_key = "shared",
      consensus_count = 4, junction = "ATG", junction_aa = "M", v_gene_key = "TRBV1",
      j_gene_key = "TRBJ1-1", representative_sequence_id = paste0(libraries[index], "~t1")),
      group_files[index])
  }
  write_pipeline_table(data.frame(sequence_id = c("t1", "t9"), sequence = c("TTTT", "CCCC")), airr_files[1])
  write_pipeline_table(data.frame(sequence_id = "t1", sequence = "AAAA"), airr_files[2])
  membership <- file.path(directory, "members.tsv")
  writeLines("library_id\tclone_id\tconsensus_count", membership)
  records <- file.path(directory, "records.tsv")
  writeLines(paste(c("library_id", "clone_id", "sequence_id", "consensus_count", "junction",
                     "junction_aa", "v_gene", "j_gene", "sequence"), collapse = "\t"), records)
  reference <- file.path(directory, "reference_metadata.tsv")
  writeLines("allele\tidentity_id\tidentity_label", reference)
  output <- file.path(directory, c("pairs.tsv", "summary.tsv", "shared.tsv"))
  options <- list(membership_file = membership, records_file = records, scope_id = "scope",
                  donor_id = "donor", library_ids = libraries, tissue_ids = c("T1", "T2"),
                  loci = "TRB", minimum_consensus_count = list(), trb_files = group_files,
                  airr_files = airr_files, donor_overlap_pairs_file = output[1],
                  donor_overlap_file = output[2], shared_groups_file = output[3],
                  reference_metadata_file = reference)
  run_stage_with_options("calculate_donor_overlap.R", options)
  shared <- read_pipeline_table(output[3])
  expect_identical(shared$library_id, libraries)
  expect_identical(unique(shared$representative_sequence_id), "B~t1")
  # A representative absent from its AIRR table is a broken link.
  write_pipeline_table(data.frame(sequence_id = "t9", sequence = "CCCC"), airr_files[1])
  run <- run_stage_with_options("calculate_donor_overlap.R", options, expect_success = FALSE)
  expect_gt(run$status, 0L)
  expect_match(paste(run$output, collapse = "\n"), "name representatives absent")
})
