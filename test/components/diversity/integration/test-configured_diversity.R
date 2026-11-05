# Run configured loci and minima through diversity, curves and support vectors.
# Run: python test/run_tests.py --component diversity --level integration

load_test_modules(c("table_io"))

test_that("configured loci and minima reach diversity, curves and support vectors", {
  directory <- tempfile("clonotype-diversity-")
  dir.create(directory)
  on.exit(unlink(directory, recursive = TRUE), add = TRUE)
  loci <- c("IGK", "IGL", "TRA", "TRB", "TRD", "TRG")
  inputs <- setNames(as.list(file.path(directory, paste0(loci, ".tsv"))), tolower(loci))
  for (locus in loci) {
    counts <- switch(locus, IGK = c(5, 3), IGL = 7, TRB = c(1, 5), numeric())
    records <- data.frame(clonotype_key = if (length(counts)) paste0(locus, seq_along(counts)) else character(),
      n_members = seq_along(counts), consensus_count = counts,
      v_gene_key = if (locus == "IGK") c("IGKV1", "VSEQ_example") else rep(paste0(locus, "V1"), length(counts)),
      j_gene_key = rep(paste0(locus, "J1"), length(counts)))
    write_pipeline_table(records, inputs[[tolower(locus)]])
  }
  metadata <- file.path(directory, "reference_metadata.tsv")
  write_pipeline_table(data.frame(
    allele = c("IGKV1-39*01", "IGKV1D-39*01"), gene = c("IGKV1-39", "IGKV1D-39"),
    region = "V-REGION", identity_id = "VSEQ_example", identity_label = "IGKV1-39,IGKV1D-39"),
    metadata)
  execute <- function(name, targets, minimums) {
    paths <- setNames(as.list(file.path(directory, paste0(name, "_",
      c("diversity", "curves", "group_abundance", "gene_usage", "vj_pair_usage"), ".tsv"))),
      c("diversity", "curves", "group_abundance", "gene_usage", "vj_pair_usage"))
    run_stage_with_options("calculate_library_clonotype_diversity.R", c(
      setNames(inputs, paste0(names(inputs), "_file")),
      list(library_id = "0040", reference_metadata_file = metadata,
           rarefaction_targets = targets,
           minimum_consensus_count = minimums),
      setNames(paths, paste0(names(paths), "_file"))))
    lapply(paths, read_pipeline_table)
  }
  ordinary <- execute("ordinary", list(TRB = 4), list())
  configured <- execute("configured", list(IGK = 4, IGL = 8, TRB = 4), list(TRB = 7))
  # Every locus has a diversity row and a support vector, and a target adds a
  # curve. Loci follow text order, the order scope collection merges them in.
  expect_identical(names(ordinary$diversity), c(
    "library_id", "locus", "group_unit", "grouped_consensus_count",
    "min_consensus_count_for_diversity", "diversity_status", "richness",
    "shannon_entropy", "hill_shannon", "inverse_simpson", "clonality",
    "top_1_support_fraction", "top_5_support_fraction", "top_10_support_fraction",
    "rarefaction_target", "rarefaction_status", "rarefied_richness",
    "rarefied_hill_shannon", "rarefied_inverse_simpson",
    "inverse_simpson_asymptote", "asymptote_status"))
  expect_identical(ordinary$diversity$locus, loci)
  expect_identical(configured$diversity$locus, loci)
  expect_identical(unique(ordinary$group_abundance$locus), c("IGK", "IGL", "TRB"))
  expect_false(is.unsorted(ordinary$group_abundance$locus))
  expect_identical(configured$group_abundance, ordinary$group_abundance)
  # Each support vector carries its clonotype's member count unchanged, in
  # support-rank order, so TRB2 with 5 support units precedes TRB1 with 1.
  expect_identical(ordinary$group_abundance$group_key,
                   c("IGK1", "IGK2", "IGL1", "TRB2", "TRB1"))
  expect_identical(ordinary$group_abundance$group_rank, c("1", "2", "1", "1", "2"))
  expect_identical(ordinary$group_abundance$n_members, c("1", "2", "1", "2", "1"))
  expect_setequal(ordinary$curves$locus, "TRB")
  expect_setequal(configured$curves$locus, c("IGK", "IGL", "TRB"))
  untargeted <- ordinary$diversity[ordinary$diversity$locus %in% c("IGK", "IGL"), ]
  expect_identical(untargeted$rarefaction_status, c("not_configured", "not_configured"))
  expect_identical(untargeted$rarefaction_target, c("", ""))
  expect_identical(untargeted$rarefied_inverse_simpson, c("", ""))
  expect_identical(untargeted$diversity_status, c("available", "available"))
  expect_equal(as.numeric(untargeted$inverse_simpson), c(64 / 34, 1))
  expect_identical(configured$gene_usage, ordinary$gene_usage)
  expect_setequal(ordinary$gene_usage$locus, c("IGK", "IGL", "TRB"))
  names <- setNames(ordinary$gene_usage$gene_names, ordinary$gene_usage$gene_key)
  expect_identical(unname(names[c("VSEQ_example", "IGKV1")]), c("IGKV1-39,IGKV1D-39", "IGKV1"))
  groups <- setNames(ordinary$gene_usage$gene_equivalence_group, ordinary$gene_usage$gene_key)
  expect_identical(unname(groups[c("VSEQ_example", "IGKV1")]), c("IGKV1-39,IGKV1D-39", "IGKV1"))
  kappa_v <- ordinary$gene_usage[ordinary$gene_usage$locus == "IGK" &
                                     ordinary$gene_usage$gene_type == "V", ]
  expect_equal(as.numeric(kappa_v$support_fraction[kappa_v$gene_key == "IGKV1"]), 5 / 8)
  # Pair usage has one row per observed V and J pair, and a composite counts once
  # under its own key.
  expect_identical(configured$vj_pair_usage, ordinary$vj_pair_usage)
  expect_identical(names(ordinary$vj_pair_usage), c(
    "library_id", "locus", "group_unit", "v_gene_key", "v_gene_names", "j_gene_key", "j_gene_names",
    "n_groups", "group_fraction", "consensus_count", "support_fraction"))
  expect_setequal(ordinary$vj_pair_usage$locus, c("IGK", "IGL", "TRB"))
  kappa_pairs <- ordinary$vj_pair_usage[ordinary$vj_pair_usage$locus == "IGK", ]
  expect_identical(kappa_pairs$v_gene_key, c("IGKV1", "VSEQ_example"))
  expect_identical(kappa_pairs$v_gene_names, c("IGKV1", "IGKV1-39,IGKV1D-39"))
  expect_identical(kappa_pairs$j_gene_key, c("IGKJ1", "IGKJ1"))
  expect_equal(as.numeric(kappa_pairs$group_fraction), c(.5, .5))
  expect_equal(as.numeric(kappa_pairs$support_fraction), c(5, 3) / 8)
  expect_true(all(ordinary$vj_pair_usage$group_unit == "clonotype"))
  expect_true(all(configured$diversity$library_id == "0040"))
  kappa <- configured$diversity[configured$diversity$locus == "IGK", ]
  expect_equal(as.numeric(kappa$grouped_consensus_count), 8)
  expect_equal(as.numeric(kappa$inverse_simpson), 64 / 34)
  expect_identical(kappa$rarefaction_status, "available")
  lambda <- configured$diversity[configured$diversity$locus == "IGL", ]
  expect_identical(lambda$rarefaction_status, "below_target_support")
  expect_equal(as.numeric(lambda$inverse_simpson), 1)
  expect_identical(lambda$rarefied_inverse_simpson, "")
  trb <- configured$diversity[configured$diversity$locus == "TRB", ]
  expect_identical(trb$diversity_status, "below_minimum_support")
  expect_identical(trb$rarefaction_status, "below_minimum_support")
  expect_identical(trb$inverse_simpson, "")
  expect_identical(trb$top_1_support_fraction, "")
  expect_equal(as.numeric(trb$grouped_consensus_count), 6)
  expect_equal(as.numeric(trb$richness), 2)
  gated_curve <- configured$curves[configured$curves$locus == "TRB", ]
  original_curve <- ordinary$curves[ordinary$curves$locus == "TRB", ]
  rownames(gated_curve) <- rownames(original_curve) <- NULL
  expect_identical(gated_curve, original_curve)
  for (locus in c("IGK", "IGL", "TRB")) {
    expect_equal(sum(as.numeric(configured$group_abundance$consensus_count[
      configured$group_abundance$locus == locus])),
      as.numeric(configured$diversity$grouped_consensus_count[configured$diversity$locus == locus]))
  }
})

test_that("the IGH arm writes gene and pair usage from each clone's vote", {
  directory <- tempfile("clone-diversity-")
  dir.create(directory)
  on.exit(unlink(directory, recursive = TRUE), add = TRUE)
  manifest <- file.path(directory, "manifest.tsv")
  write_pipeline_table(data.frame(library_id = c("0460", "NA"), donor_id = "D1",
                                  tissue_id = c("T1", "T2")), manifest)
  membership <- file.path(directory, "membership.tsv")
  write_pipeline_table(data.frame(
    library_id = c("0460", "0460", "0460", "NA"), sequence_id = paste0("s", 1:4),
    clone_id = c("c1", "c1", "c2", "c1"), consensus_count = c(4, 2, 3, 5),
    v_call = c("IGHV1-2*01", "IGHV1-2*01", "IGHV3-23*01", "IGHV1-2*01"),
    j_call = c("IGHJ4*02", "IGHJ4*02", "IGHJ6*01", "IGHJ4*02"),
    v_gene = c("VSEQ_x", "VSEQ_x", "IGHV3-23", "VSEQ_x"),
    j_gene = c("IGHJ4", "IGHJ4", "IGHJ6", "IGHJ4"),
    junction = c("TGTAAA", "TGTAAG", "TGTCCC", "TGTAAA")), membership)
  metadata <- file.path(directory, "reference_metadata.tsv")
  write_pipeline_table(data.frame(
    allele = c("IGHV1-2*01", "IGHV1-2D*01"), gene = c("IGHV1-2", "IGHV1-2D"),
    region = "V-REGION", identity_id = "VSEQ_x", identity_label = "IGHV1-2,IGHV1-2D"), metadata)
  products <- c("curves", "clone_diversity", "group_abundance", "gene_usage", "vj_pair_usage")
  paths <- setNames(as.list(file.path(directory, paste0(products, ".tsv"))), products)
  run_stage_with_options("calculate_scope_clone_diversity.R", c(
    list(manifest_file = manifest, membership_files = membership, scope_id = "scope",
         reference_metadata_file = metadata, rarefaction_targets = list(IGH = 3),
         minimum_consensus_count = list()),
    setNames(paths, paste0(names(paths), "_file"))))
  pairs <- read_pipeline_table(paths$vj_pair_usage)
  expect_identical(pairs$library_id, c("0460", "0460", "NA"))
  expect_true(all(pairs$locus == "IGH" & pairs$group_unit == "clone"))
  expect_identical(pairs$v_gene_key, c("IGHV3-23", "VSEQ_x", "VSEQ_x"))
  expect_identical(pairs$v_gene_names, c("IGHV3-23", "IGHV1-2,IGHV1-2D", "IGHV1-2,IGHV1-2D"))
  expect_identical(pairs$j_gene_key, c("IGHJ6", "IGHJ4", "IGHJ4"))
  # A clone counts once in each library that observed it, with that library's support.
  expect_identical(pairs$n_groups, c("1", "1", "1"))
  expect_identical(pairs$consensus_count, c("3", "6", "5"))
  expect_equal(as.numeric(pairs$group_fraction), c(.5, .5, 1))
  expect_equal(as.numeric(pairs$support_fraction), c(3, 6, 9) / 9)
  usage <- read_pipeline_table(paths$gene_usage)
  v_usage <- usage[usage$gene_type == "V", ]
  expect_identical(v_usage$gene_key, c("VSEQ_x", "IGHV3-23", "VSEQ_x"))
  expect_identical(v_usage$consensus_count, c("6", "3", "5"))
})
