# Check effective sequences, regional opportunities and pinned BASELINe estimates.
# Run: python test/run_tests.py --component selection --level unit

load_test_modules(c(
  "table_io", "gene_calls", "isotype_resolution", "baseline_selection", "shm_regions"))

test_that("mask_unobserved_alignment_gaps keeps shared spacers and masks uncovered bases", {
  expect_identical(mask_unobserved_alignment_gaps("...AC.-T", "...ACGGT"),
                   "...ACNNT")
})

describe("mask_unobserved_alignment_gaps over a vector of alignments", {
  cases <- list(
    list(sequence = "A.C", germline = "AGC", expected = "ANC",
         why = "masks a dot that covers a germline base"),
    list(sequence = "A-C", germline = "AGC", expected = "ANC",
         why = "masks a dash that covers a germline base"),
    list(sequence = "...", germline = "...", expected = "...",
         why = "leaves a gap the germline shares untouched")
  )
  for (case in cases) {
    it(case$why, expect_identical(
      mask_unobserved_alignment_gaps(case$sequence, case$germline),
      case$expected))
  }
})


test_that("selection frequency is one finite number in (0, 1]", {
  for (value in list(1e-12, 0.5, 0.6, 1, 1L)) {
    expect_identical(validate_selection_collapse_frequency(value), value)
  }
  for (value in list(0, -0.1, 1.01, NA_real_, NaN, Inf, -Inf,
                    NULL, numeric(), c(0.5, 0.6), "0.6", TRUE, list(0.6), 0.6 + 0i)) {
    expect_error(validate_selection_collapse_frequency(value),
                 "collapse_min_freq must be one finite number in (0, 1]", fixed = TRUE)
  }
})

test_that("the collapse adapter validates frequency before using its inputs", {
  expect_error(collapse_clones_to_effective_sequences(data.frame(), 0, 1),
               "collapse_min_freq must be one finite number in (0, 1]", fixed = TRUE)
})

test_that("selection masks incomplete codons without shifting the IMGT frame", {
  expect_identical(mask_incomplete_selection_codons("GCTANNGCT", "GCTGCTGCT"),
                   "GCTNNNGCT")
  expect_identical(mask_incomplete_selection_codons("GCT...GCT", "GCT...GCT"),
                   "GCT...GCT")
  expect_identical(mask_incomplete_selection_codons("GC.GCT", "GC.GCT"), "NN.GCT")
  expect_identical(mask_incomplete_selection_codons("GCTAC", "GCTAC"), "GCTNN")
  expect_identical(mask_incomplete_selection_codons("GCTGCT", "GCTGCN"), "GCTNNN")
  expect_identical(mask_incomplete_selection_codons("GCTG-T", "GCTGCT"), "GCTNNN")
  expect_identical(mask_incomplete_selection_codons(character(), character()), character())
  expect_error(mask_incomplete_selection_codons("GCT", "GCTGCT"), "equal string lengths")
})

test_that("the released column order is fixed", {
  expect_identical(SELECTION_COLUMNS, c(
    "scope_id", "library_id", "donor_id", "tissue_id", "isotype_group", "n_clone_consensus_sequences",
    "n_clone_consensus_sequences_with_positions_cdr", "n_clone_consensus_sequences_used_cdr", "sigma_cdr_status",
    "sigma_cdr", "sigma_cdr_ci_low", "sigma_cdr_ci_high",
    "n_clone_consensus_sequences_with_positions_fwr", "n_clone_consensus_sequences_used_fwr", "sigma_fwr_status",
    "sigma_fwr", "sigma_fwr_ci_low", "sigma_fwr_ci_high"))
})

test_that("selection removes the partial-codon replacement counted by pinned SHazaM", {
  germline <- strrep("GCT", 104)
  partial <- germline
  substr(partial, 103, 105) <- "ANN"
  observed <- function(sequence) shazam::calcObservedMutations(
    sequence, germline, regionDefinition = shazam::IMGT_V, returnRaw = TRUE)
  unmasked <- observed(partial)
  expect_identical(as.integer(unmasked$pos$position), 103L)
  expect_equal(unmasked$pos$r, 1)
  masked <- mask_incomplete_selection_codons(partial, germline)
  expect_identical(substr(masked, 103, 105), "NNN")
  expect_true(all(is.na(observed(masked)$pos)))
  targeting <- shazam:::calculateTargeting(germlineSeq = germline, inputSeq = partial,
    targetingModel = shazam::HH_S5F, regionDefinition = shazam::IMGT_V)
  paths <- shazam:::calculateMutationalPaths(
    germlineSeq = paste(colnames(targeting), collapse = ""),
    regionDefinition = shazam::IMGT_V, codonTable = shazam:::CODON_TABLE)
  expect_false(any(paths[, 103:105] %in% c("r", "s")))
  contigs <- data.frame(library_id = "s1", isotype_group = "all_IGH", clone_id = "c1",
    sequence_alignment = partial, germline_alignment_d_mask = germline)
  effective <- collapse_clones_to_effective_sequences(contigs, 0.6, 1)$effective_sequences
  expect_identical(effective$clonal_sequence, masked)
  result <- suppressWarnings(estimate_baseline_selection(effective, "all_IGH", 1))
  expect_identical(result$pdf_counts$n_pdf_cdr, 0L)
  expect_identical(result$pdf_counts$n_pdf_fwr, 0L)
  # The SHM endpoint retains its separate incomplete-sequence policy.
  expect_equal(summarise_observed_mutations_by_region(partial, germline)$mutations["cdr1"],
               c(cdr1 = 1L))
})

test_that("selection context counts the effective sequences with a callable position", {
  germline <- strrep("GCT", 104)
  half_cdr <- germline
  substr(half_cdr, 79, 111) <- strrep("N", 33)
  substr(half_cdr, 112, 112) <- "A"
  no_cdr <- germline
  positions <- which(as.character(shazam::IMGT_V@boundaries) == "cdr")
  chars <- strsplit(no_cdr, "", fixed = TRUE)[[1]]
  chars[positions] <- "N"
  no_cdr <- paste(chars, collapse = "")
  representatives <- data.frame(clonal_sequence = c(germline, half_cdr, no_cdr),
                                 clonal_germline = germline)
  # The masked CDR has no callable position; the half-masked one keeps some.
  row <- summarise_selection_context(representatives, "cdr")
  expect_identical(names(row), "n_clone_consensus_sequences_with_positions_cdr")
  expect_equal(row$n_clone_consensus_sequences_with_positions_cdr, 2L)
  expect_equal(summarise_selection_context(representatives, "fwr")$n_clone_consensus_sequences_with_positions_fwr, 3L)
  missing <- summarise_selection_context(NULL, "cdr")
  expect_equal(missing$n_clone_consensus_sequences_with_positions_cdr, 0L)
})

# A focused test needs exactly the two regions; calcBaseline otherwise falls back
# to a local test without saying so.
test_that("the region partition is the CDR/FWR split of IMGT V", {
  expect_setequal(as.character(selection_region_definition()@regions),
                  c("cdr", "fwr"))
  expect_identical(selection_region_definition()@seqLength,
                   shazam::IMGT_V@seqLength)
})

# -------------------------------------------------------------- eligibility ---

describe("selection_record_eligible admits what an estimate needs", {
  it("admits a record with a lineage and a germline", {
    expect_true(selection_record_eligible("cohort_D1_1", "ACGT", "ACGT"))
  })
  it("refuses a record with no lineage", {
    expect_false(selection_record_eligible(NA_character_, "ACGT", "ACGT"))
    expect_false(selection_record_eligible("", "ACGT", "ACGT"))
  })
  it("refuses a record with no D-masked germline", {
    expect_false(selection_record_eligible("cohort_D1_1", "ACGT", ""))
  })
  it("refuses a record with no alignment", {
    expect_false(selection_record_eligible("cohort_D1_1", NA_character_, "ACGT"))
  })
  # Table writers use empty cells for missing values. The literal identifier
  # NA remains a string when selection reads the pooled IGH table.
  it("treats a cell holding the literal NA as the string it is", {
    expect_true(selection_record_eligible("NA", "ACGT", "ACGT"))
  })
  it("is vectorised", {
    expect_identical(
      selection_record_eligible(c("L1", "", "L3"), "ACGT", c("ACGT", "ACGT", "")),
      c(TRUE, FALSE, FALSE))
  })
})

# ----------------------------------------------------------- stratification ---

RECORDS <- data.frame(
  sequence_id = paste0("s", 1:4),
  isotype_group = c("IgG", "IgM", NA_character_, "IgG"),
  stringsAsFactors = FALSE)

describe("build_isotype_grouped_selection_input stratifies the repertoire", {
  grouped <- build_isotype_grouped_selection_input(RECORDS)

  # Every record answers the C-independent question, and a resolved one answers
  # its isotype's question as well: the isotype rows are separate readings of the
  # same repertoire, not a partition of the all-IGH row.
  it("gives every record an all_IGH copy",
     expect_identical(sum(grouped$isotype_group == "all_IGH"), 4L))
  it("gives a resolved record a second copy in its own group", {
    expect_identical(sum(grouped$isotype_group == "IgG"), 2L)
    expect_identical(sum(grouped$isotype_group == "IgM"), 1L)
  })
  it("keeps an unresolved C call in all_IGH only", {
    expect_identical(
      grouped$isotype_group[grouped$sequence_id == "s3"], "all_IGH")
  })
  it("counts one record twice at most",
     expect_identical(nrow(grouped), 7L))
  it("handles a pool with nothing resolved", {
    unresolved <- RECORDS
    unresolved$isotype_group <- NA_character_
    expect_identical(nrow(build_isotype_grouped_selection_input(unresolved)), 4L)
  })
  # The helper requires a populated input. calculate_donor_selection.R skips
  # stratification for an empty eligible pool and writes the full library grid;
  # test-selection.R checks that entry-point behaviour.
  it("stops on a pool with no eligible record", {
    expect_error(build_isotype_grouped_selection_input(RECORDS[0, ]),
                 "no eligible IGH record to stratify")
  })
})

# --------------------------------------------------- reading off an estimate ---

SUMMARY <- data.frame(
  library_id = c("s1", "s1", "s2"),
  region = c("cdr", "fwr", "cdr"),
  baseline_sigma = c(0.4, -0.2, NaN),
  baseline_ci_lower = c(0.1, -0.5, NA_real_),
  baseline_ci_upper = c(0.7, 0.1, NA_real_),
  stringsAsFactors = FALSE)

PDF_ROW <- data.frame(library_id = "s1", n_pdf_cdr = 3L, n_pdf_fwr = 0L,
                      stringsAsFactors = FALSE)
RESULT <- list(summary = SUMMARY, pdf_counts = PDF_ROW)

describe("summarise_region_for_library reads one library's estimate", {
  it("names the released columns for its region", {
    expect_identical(names(summarise_region_for_library(RESULT, "s1", "cdr", 5L)),
                     setdiff(selection_region_columns("cdr"), "n_clone_consensus_sequences_with_positions_cdr"))
    expect_identical(selection_region_columns("fwr"),
                     c("n_clone_consensus_sequences_with_positions_fwr", "n_clone_consensus_sequences_used_fwr", "sigma_fwr_status",
                       "sigma_fwr", "sigma_fwr_ci_low", "sigma_fwr_ci_high"))
  })
  it("takes the sigma, its interval and the density count", {
    cdr <- summarise_region_for_library(RESULT, "s1", "cdr", 5L)
    expect_equal(cdr$sigma_cdr, 0.4)
    expect_equal(cdr$sigma_cdr_ci_low, 0.1)
    expect_equal(cdr$sigma_cdr_ci_high, 0.7)
    expect_identical(cdr$n_clone_consensus_sequences_used_cdr, 3L)
    expect_identical(cdr$sigma_cdr_status, "estimated")
  })
  # A sigma with no density retained behind it is not an estimate either.
  it("refuses a sigma with no density behind it", {
    fwr <- summarise_region_for_library(RESULT, "s1", "fwr", 5L)
    expect_identical(fwr$sigma_fwr_status, "no_usable_mutation_pdfs")
    expect_identical(fwr$n_clone_consensus_sequences_used_fwr, 0L)
    expect_identical(fwr$sigma_fwr, NA_real_)
  })
  # A non-finite sigma is not an estimate, and must not be reported as one.
  it("refuses a non-finite sigma", {
    row <- summarise_region_for_library(RESULT, "s2", "cdr", 5L)
    expect_identical(row$sigma_cdr_status, "no_usable_mutation_pdfs")
    expect_identical(row$sigma_cdr, NA_real_)
  })
  # A library with no lineages in the group really used zero densities, which is a
  # known count; a library that had lineages and is absent from the estimate is a
  # case nothing was reported about, and zero would read as an answer.
  it("counts zero densities for a library that contributed nothing", {
    row <- summarise_region_for_library(RESULT, "s9", "cdr", 0L)
    expect_identical(row$sigma_cdr_status, "no_eligible_clones")
    expect_identical(row$n_clone_consensus_sequences_used_cdr, 0L)
  })
  it("reports no count for a library that contributed but was not summarised", {
    row <- summarise_region_for_library(RESULT, "s9", "cdr", 2L)
    expect_identical(row$sigma_cdr_status, "no_usable_mutation_pdfs")
    expect_identical(row$n_clone_consensus_sequences_used_cdr, NA_integer_)
  })
  it("handles a group that produced no estimate at all", {
    expect_identical(summarise_region_for_library(NULL, "s1", "cdr", 0L)$sigma_cdr_status,
                     "no_eligible_clones")
    expect_identical(summarise_region_for_library(NULL, "s1", "cdr", 4L)$sigma_cdr_status,
                     "no_usable_mutation_pdfs")
  })
})

# ------------------------------------------------------------- the released rows ---

SAMPLES <- data.frame(library_id = c("s1", "s2"), donor_id = "D1",
                      tissue_id = c("tumour", "normal"), stringsAsFactors = FALSE)
REPRESENTATIVES <- data.frame(
  library_id = c("s1", "s1", "s1"),
  isotype_group = c("all_IGH", "all_IGH", "IgG"),
  clone_id = c("cohort_D1_1", "cohort_D1_2", "cohort_D1_1"),
  clonal_sequence = strrep("GCT", 104), clonal_germline = strrep("GCT", 104),
  stringsAsFactors = FALSE)
BASELINE_RESULTS <- list(
  all_IGH = list(summary = SUMMARY, pdf_counts = PDF_ROW),
  IgM = NULL,
  IgG = NULL,
  IgA = NULL)

describe("build_selection_result_rows reports every library in every group", {
  rows <- build_selection_result_rows(BASELINE_RESULTS, REPRESENTATIVES,
                                      SAMPLES, "cohort")

  it("writes one row per library and isotype group",
     expect_identical(nrow(rows), 8L))
  it("declares the released columns in order",
     expect_identical(names(rows), SELECTION_COLUMNS))
  it("reports the groups in their released order",
     expect_identical(unique(rows$isotype_group), ANALYSIS_IGH_ISOTYPE_GROUPS))
  it("stamps the scope and carries the manifest's donor and tissue", {
    expect_identical(unique(rows$scope_id), "cohort")
    expect_identical(rows$tissue_id[rows$library_id == "s2"], rep("normal", 4L))
  })
  it("counts the effective sequences each library contributed", {
    expect_identical(rows$n_clone_consensus_sequences[rows$library_id == "s1"], c(2L, 0L, 1L, 0L))
    expect_true(all(rows$n_clone_consensus_sequences[rows$library_id == "s2"] == 0L))
  })
  it("reports the estimate where there is one", {
    row <- rows[rows$library_id == "s1" & rows$isotype_group == "all_IGH", ]
    expect_equal(row$sigma_cdr, 0.4)
    expect_equal(row$sigma_cdr_ci_low, 0.1)
    expect_identical(row$sigma_cdr_status, "estimated")
    expect_identical(row$n_clone_consensus_sequences_used_cdr, 3L)
    # The FWR estimate exists but no density was retained for it.
    expect_identical(row$sigma_fwr_status, "no_usable_mutation_pdfs")
    expect_identical(row$n_clone_consensus_sequences_used_fwr, 0L)
  })
  it("says why a group with lineages has no estimate", {
    row <- rows[rows$library_id == "s1" & rows$isotype_group == "IgG", ]
    expect_identical(row$sigma_cdr_status, "no_usable_mutation_pdfs")
    expect_identical(row$n_clone_consensus_sequences_used_cdr, NA_integer_)
    expect_true(is.na(row$sigma_cdr))
  })
  it("says why a library with nothing has no estimate", {
    row <- rows[rows$library_id == "s2" & rows$isotype_group == "all_IGH", ]
    expect_identical(row$sigma_cdr_status, "no_eligible_clones")
    expect_identical(row$sigma_fwr_status, "no_eligible_clones")
    expect_identical(row$n_clone_consensus_sequences_used_cdr, 0L)
  })
})

test_that("a donor whose libraries produced nothing still reports every row", {
  rows <- build_selection_result_rows(
    stats::setNames(vector("list", 4L), ANALYSIS_IGH_ISOTYPE_GROUPS),
    NULL, SAMPLES, "cohort")
  expect_identical(nrow(rows), 8L)
  expect_true(all(rows$n_clone_consensus_sequences == 0L))
  expect_true(all(rows$sigma_cdr_status == "no_eligible_clones"))
  expect_true(all(is.na(rows$sigma_cdr)))
})

# ------------------------------------------- the collapse and the estimate ---
#
# One real IMGT V region: an all-A germline, one sequence with a replacement
# mutation in a CDR codon, and one whose CDR was never observed and is dot-padded.

GERMLINE <- paste(rep("A", shazam::IMGT_V@seqLength), collapse = "")
V_REGIONS <- as.character(shazam::IMGT_V@boundaries)
stopifnot(V_REGIONS[80] == "cdr", V_REGIONS[6] == "fwr")
replace_base <- function(sequence, position, base) {
  chars <- strsplit(sequence, "", fixed = TRUE)[[1]]
  chars[position] <- base
  paste(chars, collapse = "")
}
# AAA -> ACA at the middle base of a CDR codon is a replacement mutation.
CDR_MUTATED <- replace_base(GERMLINE, 80L, "C")
# Every CDR position unobserved, with one silent FWR mutation, so this sequence
# can only ever speak about the framework.
dot_chars <- strsplit(GERMLINE, "", fixed = TRUE)[[1]]
dot_chars[V_REGIONS == "cdr"] <- "."
dot_chars[6] <- "G"
DOT_CDR <- paste(dot_chars, collapse = "")

CONTIGS <- data.frame(
  library_id = c("s1", "s1", "s2"),
  isotype_group = "all_IGH",
  clone_id = c("cohort_D1_1", "cohort_D1_2", "cohort_D1_1"),
  sequence_alignment = c(CDR_MUTATED, DOT_CDR, CDR_MUTATED),
  germline_alignment_d_mask = GERMLINE,
  stringsAsFactors = FALSE)

describe("collapse_clones_to_effective_sequences reduces a lineage to one sequence", {
  collapsed <- collapse_clones_to_effective_sequences(CONTIGS, 0.6, 1)
  effective <- collapsed$effective_sequences

  # One lineage shared by two libraries yields one effective sequence in each,
  # because the library is an explicit grouping field.
  it("keeps a donor-pooled lineage separate by library",
     expect_identical(nrow(effective), 3L))
  it("returns the collapsed sequence and its germline",
     expect_identical(names(effective), c("library_id", "isotype_group", "clone_id",
                                          "clonal_sequence", "clonal_germline")))
  # An unobserved CDR opposite a germline base becomes N, so SHazaM does not
  # count it as mutation opportunity that failed to mutate.
  it("masks the sequence that was never observed",
     expect_identical(collapsed$n_masked, 1L))
  it("masks it with N rather than dropping it", {
    dotted <- effective$clonal_sequence[effective$clone_id == "cohort_D1_2"]
    expect_false(grepl(".", dotted, fixed = TRUE))
    expect_true(grepl("N", dotted, fixed = TRUE))
    expect_identical(nchar(dotted), nchar(GERMLINE))
  })
  it("leaves an observed sequence alone", {
    expect_identical(
      effective$clonal_sequence[effective$clone_id == "cohort_D1_1" &
                                  effective$library_id == "s1"],
      CDR_MUTATED)
  })
})

test_that("positive frequencies preserve agreed bases and unobserved positions", {
  germline <- strrep("GCT", 104)
  variant <- germline
  substr(variant, 79, 79) <- "A"
  cases <- list(identical = c(germline, germline), variable = c(germline, variant),
    clipped = paste0("...", substring(c(germline, variant), 4)),
    outside_v = paste0(germline, c("GCT", "ACT")))
  for (name in names(cases)) {
    contigs <- data.frame(library_id = "library", isotype_group = "all_IGH",
      clone_id = "clone", sequence_alignment = cases[[name]],
      germline_alignment_d_mask = if (name == "outside_v") paste0(germline, "GCT") else germline)
    # Even identical members must reject zero before SHazaM's early return.
    expect_error(collapse_clones_to_effective_sequences(contigs, 0, 1),
                 "collapse_min_freq must be one finite number in (0, 1]", fixed = TRUE)
    for (frequency in c(1e-12, 0.4, 0.5, 0.6, 1)) {
      expected <- germline
      if (name %in% c("variable", "clipped")) {
        # Both A and G have frequency 1/2; A wins the fixed-order tie.
        if (frequency <= 0.5) substr(expected, 79, 79) <- "A" else
          substr(expected, 79, 81) <- "NNN"
      }
      if (name == "clipped") substr(expected, 1, 3) <- "NNN"
      effective <- collapse_clones_to_effective_sequences(contigs, frequency, 1)$effective_sequences
      expect_identical(effective$clonal_sequence, expected, info = paste(name, frequency))
      expect_identical(effective$clonal_germline, germline, info = paste(name, frequency))
    }
  }
})

# Truncated members pad the positions they did not observe. At the default 0.6
# frequency, a base that one of three members covers does not qualify, so the
# effective sequence marks the position N rather than taking that base. Input
# order does not change the result.
test_that("a base covered by a minority of members is not called at 0.6", {
  truncated <- paste0(strrep(".", 90), substring(GERMLINE, 91))
  contigs <- data.frame(library_id = "s1", isotype_group = "all_IGH",
    clone_id = "cohort_D1_1", sequence_alignment = c(CDR_MUTATED, truncated, truncated),
    germline_alignment_d_mask = GERMLINE, stringsAsFactors = FALSE)
  alone <- collapse_clones_to_effective_sequences(contigs[1, ], 0.6, 1)$effective_sequences
  mixed <- collapse_clones_to_effective_sequences(contigs, 0.6, 1)$effective_sequences
  expect_identical(substr(alone$clonal_sequence, 80, 80), "C")
  expect_identical(mixed$clonal_sequence, paste0(strrep("N", 90), substring(GERMLINE, 91)))
  expect_identical(mixed$clonal_germline, GERMLINE)
  expect_identical(
    collapse_clones_to_effective_sequences(contigs[3:1, ], 0.6, 1)$effective_sequences, mixed)
})

test_that("effective bases follow the configured frequency and deterministic ties", {
  contigs <- CONTIGS[rep(1, 5), , drop = FALSE]
  # The first three positions have largest base frequencies 2/5, 3/5 and 4/5.
  # SHazaM calls N below the floor and includes a base exactly at the floor.
  # Selection then masks the whole incomplete codon.
  contigs$sequence_alignment <- paste0(c("AAA", "AAA", "TAA", "TTA", "CCC"),
                                       substring(GERMLINE, 4))
  effective <- collapse_clones_to_effective_sequences(contigs, 0.6, 1)$effective_sequences
  expect_equal(nrow(effective), 1L)
  expect_identical(effective$clonal_sequence, paste0("NNN", substring(GERMLINE, 4)))
  expect_identical(effective$clonal_germline, GERMLINE)

  # At a 0.4 floor A and T tie at the first position. With both ambiguity and
  # stochastic ties disabled, the documented package ordering selects A.
  tied <- collapse_clones_to_effective_sequences(contigs, 0.4, 1)$effective_sequences
  expect_identical(tied$clonal_sequence, GERMLINE)
})

describe("estimate_baseline_selection convolves one group's densities", {
  effective <- collapse_clones_to_effective_sequences(
    CONTIGS, 0.6, 1)$effective_sequences
  # SHazaM notices that this synthetic germline carries no IMGT gap. A real
  # germline does; the fixture is a flat run of one base so that every mutation
  # in it is countable by eye. The notice is about the fixture, not the call.
  estimate <- suppressWarnings(
    estimate_baseline_selection(effective, "all_IGH", 1))

  it("summarises every library it saw",
     expect_setequal(estimate$summary$library_id, c("s1", "s2")))
  # The value, not just its finiteness. Sigma is signed and a call configured with
  # the wrong targeting model, test statistic or sequence/germline columns still
  # returns a finite number inside ordered confidence intervals, so nothing but
  # the number itself notices. This one comes from SHazaM on this fixture; the
  # tolerance is loose enough to survive a numerical change in the package and
  # tight enough that a different statistic cannot pass.
  it("reports the sigma SHazaM computes for the region that was mutated", {
    sigma <- summarise_region_for_library(estimate, "s1", "cdr", 2L)
    expect_identical(sigma$sigma_cdr_status, "estimated")
    expect_true(sigma$sigma_cdr_ci_low <= sigma$sigma_cdr &&
                  sigma$sigma_cdr <= sigma$sigma_cdr_ci_high)
    expect_equal(sigma$sigma_cdr, 2.2661557844, tolerance = 1e-6)
  })
  # The densities retained are not the number of lineages. s1 has two effective
  # sequences and each speaks about one region only: the mutated one about the
  # CDR, the masked one about the framework it did observe.
  it("counts the densities actually retained per region", {
    counts <- estimate$pdf_counts[estimate$pdf_counts$library_id == "s1", ]
    expect_identical(counts$n_pdf_cdr, 1L)
    expect_identical(counts$n_pdf_fwr, 1L)
  })
  it("leaves a library whose only sequence speaks about one region", {
    counts <- estimate$pdf_counts[estimate$pdf_counts$library_id == "s2", ]
    expect_identical(counts$n_pdf_cdr, 1L)
    expect_identical(counts$n_pdf_fwr, 0L)
  })
  it("reports nothing for a group holding no effective sequence",
     expect_null(estimate_baseline_selection(effective, "IgG", 1)))
  it("reports nothing when there are no effective sequences at all", {
    expect_null(estimate_baseline_selection(NULL, "all_IGH", 1))
    expect_null(estimate_baseline_selection(effective[0, ], "all_IGH", 1))
  })
})

# Library identifiers that readers can take for missing values or numbers must
# stay separate partitions through both SHazaM calls. Each library holds the same
# four clones: an unmutated one, one with a CDR replacement, one with an FWR
# replacement and one with both. Each region therefore has two usable densities
# in every library.
describe("selection keeps each library a separate partition", {
  libraries <- c("NA", "None", "0460", "nan", "1.0", "1")
  germline <- strrep("GCT", 104)
  variants <- vapply(list(integer(0), 79L, 4L, c(4L, 79L)), function(positions) {
    sequence <- germline
    for (position in positions) substr(sequence, position, position) <- "A"
    sequence
  }, character(1))
  contigs <- data.frame(library_id = rep(libraries, each = 4L), isotype_group = "all_IGH",
    clone_id = paste0("scope|donor|", seq_along(variants)), sequence_alignment = variants,
    germline_alignment_d_mask = germline, stringsAsFactors = FALSE)
  effective <- collapse_clones_to_effective_sequences(contigs, 0.6, 1)$effective_sequences
  pooled <- suppressWarnings(estimate_baseline_selection(effective, "all_IGH", 1))

  it("keeps literal identifiers apart", {
    expect_identical(nrow(effective), 24L)
    expect_setequal(effective$library_id, libraries)
    expect_setequal(pooled$pdf_counts$library_id, libraries)
    expect_true(all(pooled$pdf_counts$n_pdf_cdr == 2L & pooled$pdf_counts$n_pdf_fwr == 2L))
  })
  it("estimates each library as it would be estimated alone", {
    for (library in libraries) {
      alone <- suppressWarnings(estimate_baseline_selection(
        effective[effective$library_id == library, ], "all_IGH", 1))
      for (region in SELECTION_REGIONS) {
        row <- summarise_region_for_library(pooled, library, region, 4L)
        expect_identical(row[[paste0("sigma_", region, "_status")]], "estimated")
        expect_equal(row, summarise_region_for_library(alone, library, region, 4L),
                     tolerance = 1e-12, info = paste(library, region))
      }
    }
  })
  it("does not depend on input row order", {
    reversed <- collapse_clones_to_effective_sequences(
      contigs[rev(seq_len(nrow(contigs))), ], 0.6, 1)$effective_sequences
    expect_identical(reversed, effective)
    expect_identical(suppressWarnings(estimate_baseline_selection(
      effective[rev(seq_len(nrow(effective))), ], "all_IGH", 1)), pooled)
  })
})
