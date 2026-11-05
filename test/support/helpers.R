# Load explicit production modules and construct scientific test handoffs.
# Run: python test/run_tests.py

MODULE_ROOT <- normalizePath(Sys.getenv("AIRRBORNE_ROOT"), mustWork = TRUE)

load_test_modules <- function(modules) {
  target <- environment(load_test_modules)
  source(file.path(MODULE_ROOT, "workflow", "lib", "runtime.R"), local = target)
  target$set_pipeline_collation()
  for (module in modules) {
    source(file.path(MODULE_ROOT, "workflow", "lib", paste0(module, ".R")),
           local = target)
  }
  invisible(NULL)
}

run_pipeline_script <- function(script_name, arguments) {
  script <- file.path(MODULE_ROOT, "workflow", "scripts", script_name)
  rscript <- file.path(R.home("bin"), "Rscript")
  output <- suppressWarnings(system2(
    rscript, shQuote(c(script, as.character(arguments))),
    stdout = TRUE, stderr = TRUE
  ))
  status <- attr(output, "status")
  list(
    script = script_name,
    status = if (is.null(status)) 0L else as.integer(status),
    output = output
  )
}

# Assert that a script run finished cleanly. On failure the script's own output
# is reported, because that is where the reason will be.
expect_script_ok <- function(run) {
  if (identical(run$status, 0L)) {
    succeed()
  } else {
    fail(paste0(run$script, " exited ", run$status, ":\n",
                paste(run$output, collapse = "\n")))
  }
}

# Exercise both scientific branches on the same constructed library. The output
# follows the original source-row order for direct comparison with the input.
prepare_grouping_fixture <- function(records, library_id, max_n) {
  if (any(!records$locus %in% c("IGH", "IGK", "IGL", "TRA", "TRB", "TRD", "TRG"))) {
    stop("Grouping requires one of the seven receptor loci")
  }
  igh <- prepare_igh_grouping(records[records$locus == "IGH", , drop = FALSE], library_id, max_n)
  others <- lapply(c("IGK", "IGL", "TRA", "TRB", "TRD", "TRG"), function(locus) {
    build_rearrangement_group_table(records[records$locus == locus, , drop = FALSE], library_id, locus)$mapping
  })
  out <- do.call(rbind, c(list(igh), others))
  out[match(records$sequence_id, out$sequence_id), , drop = FALSE]
}

# Run a stage script with named options, as the workflow's shell blocks do.
#
# `options` is a named list. A character or numeric vector joins with commas,
# and a named list or named vector renders as NAME=value pairs, which is how
# the rules pass file lists and locus-keyed settings. An empty vector or list
# renders as an empty value. Returns run_pipeline_script()'s result and, when
# `expect_success`, fails the test on a non-zero exit.
run_stage_with_options <- function(script, options, expect_success = TRUE) {
  render <- function(value) {
    if (is.list(value) || !is.null(names(value))) {
      if (!length(value)) return("")
      return(paste0(names(value), "=", unlist(value), collapse = ","))
    }
    paste(as.character(value), collapse = ",")
  }
  arguments <- vapply(names(options), function(name) {
    paste0("--", name, "=", render(options[[name]]))
  }, character(1))
  run <- run_pipeline_script(script, arguments)
  if (expect_success) expect_script_ok(run)
  invisible(run)
}

# Construct 300 independent junction pairs for successful threshold fitting.
# Each pair has its own V gene. Half differ by 1 to 4 bases and half by 14 to
# 18 bases in a 60-base junction, leaving a known gap between the two modes.
# Every row names a separate TRUST4 assembly, without alternative CDR3 rows.
separated_distance_fixture <- function() {
  n_pairs <- 300L
  pair_index <- rep(seq_len(n_pairs), each = 2L)
  mutation_count <- c(1L + (seq_len(n_pairs / 2L) %% 4L),
                      14L + (seq_len(n_pairs / 2L) %% 5L))
  mutated <- vapply(mutation_count, function(n) {
    paste0(strrep("C", n), strrep("A", 60L - n))
  }, character(1))
  data.frame(
    sequence_id = sprintf("assemble%04d_0", seq_len(2L * n_pairs)),
    junction = as.vector(rbind(rep(strrep("A", 60L), n_pairs), mutated)),
    v_call = sprintf("IGHVTEST%04d*01", pair_index),
    j_call = "IGHJ4*01", locus = "IGH", consensus_count = 1L,
    stringsAsFactors = FALSE
  )
}

# Construct aligned IGH records and votes for unit and script measurements.
shm_feature_fixture <- function() {
  MIN_C_POSITIONS <- 40L
  MIN_POSITION_VOTES <- 2L
  # The released column set has one home, workflow/lib/igh_features.R, and the stage projects
  # its output onto that same definition. Asserting against it here is what makes
  # the stage's header and the module's contract one statement rather than two.
  EXPECTED_COLUMNS <- igh_feature_columns()
  # The IMGT V sub-region widths, which the germline of a full-length synthetic
  # record occupies completely.
  REGION_WIDTH <- c(fwr1 = 78L, cdr1 = 36L, fwr2 = 51L, cdr2 = 30L, fwr3 = 117L)
  N_POSITIONS <- as.integer(shazam::IMGT_V@seqLength)
  REGION_OF_POSITION <- imgt_v_region_of_position()
  # The stage writes a text table, and R serialises a double at 15 significant
  # digits, so a rate read back is within about 1e-15 relative of the one computed
  # in memory. The comparisons below run against SHazaM's in-memory value, so they
  # allow that and nothing wider: this is a serialisation difference, not a method
  # difference, and a real method difference would be orders of magnitude larger.
  SERIALISATION_TOLERANCE <- 1e-12

  # --------------------------------------------------------------- fixtures ---

  # A germline of 104 alanine codons. Substituting the first base of any codon
  # gives ACT (threonine), so every mutation below is a plain replacement and none
  # of them lands on a stop codon, which SHazaM would decline to score.
  GERMLINE <- paste(rep("GCT", 104), collapse = "")
  stopifnot(identical(nchar(GERMLINE), 312L))

  #' Copy the germline and substitute one base per requested position.
  mutate_positions <- function(sequence, positions, base = "A") {
    for (position in positions) substr(sequence, position, position) <- base
    sequence
  }

  #' Replace a closed range of positions with the IMGT spacer character.
  blank_positions <- function(sequence, from, to) {
    substr(sequence, from, to) <- paste(rep(".", to - from + 1L), collapse = "")
    sequence
  }

  # One distinct vote total per sub-region, so a per-region median is the region's
  # own value and a mix-up between regions cannot pass. Non-consensus is one vote
  # everywhere, so the majority is the total less one.
  REGION_TOTAL_VOTES <- c(fwr1 = 4L, cdr1 = 6L, fwr2 = 8L, cdr2 = 12L, fwr3 = 20L)
  REGION_MAJORITY_VOTES <- REGION_TOTAL_VOTES - 1L

  # Two positions carry no usable vote, for the two distinct reasons the released
  # vectors express. CDR2 166-175 have no mapped contig coordinate at all, so the
  # total itself is absent. FWR2 120 has a total but no base to agree with, which is
  # an N at a mate-pair scaffold join, so only the non-consensus count is absent.
  UNMAPPED_CDR2 <- 166:175
  SCAFFOLD_N_FWR2 <- 120L

  vote_cells <- function() {
    total <- REGION_TOTAL_VOTES[REGION_OF_POSITION]
    non_consensus <- rep(1L, N_POSITIONS)
    total[UNMAPPED_CDR2] <- NA_integer_
    non_consensus[UNMAPPED_CDR2] <- NA_integer_
    non_consensus[SCAFFOLD_N_FWR2] <- NA_integer_
    list(total = format_vote_vector(unname(total), N_POSITIONS),
         non_consensus = format_vote_vector(non_consensus, N_POSITIONS))
  }
  VOTES <- vote_cells()
  # A record read in the opposite direction to its contig is emptied by the vote
  # stage, so every position arrives absent.
  EMPTY_VOTES <- format_vote_vector(rep(NA_integer_, N_POSITIONS), N_POSITIONS)

  # The vote-count positions each sub-region contributes once the two absent kinds
  # are removed. FWR2 loses the scaffold N, CDR2 loses the ten unmapped positions.
  # Every majority above is at least two, so in this fixture a position that carries
  # a vote at all is also read-supported and this is both counts. THIN_VOTES below is
  # the fixture that separates them.
  VOTED_WIDTH <- REGION_WIDTH
  VOTED_WIDTH[["fwr2"]] <- REGION_WIDTH[["fwr2"]] - 1L
  VOTED_WIDTH[["cdr2"]] <- REGION_WIDTH[["cdr2"]] - length(UNMAPPED_CDR2)

  # A contig whose FWR1 rests on a single agreeing read, with every other sub-region
  # unchanged. TRUST4 writes a literal 1 on the consensus base wherever a position's
  # tally came out zero, so one vote is where no-evidence and one-read evidence become
  # indistinguishable; the configured vote floor draws the line above it. A total of two
  # against one dissenter leaves a majority of exactly one.
  THIN_VOTES <- local({
    total <- REGION_TOTAL_VOTES[REGION_OF_POSITION]
    total[REGION_OF_POSITION == "fwr1"] <- 2L
    non_consensus <- rep(1L, N_POSITIONS)
    total[UNMAPPED_CDR2] <- NA_integer_
    non_consensus[UNMAPPED_CDR2] <- NA_integer_
    non_consensus[SCAFFOLD_N_FWR2] <- NA_integer_
    list(total = format_vote_vector(unname(total), N_POSITIONS),
         non_consensus = format_vote_vector(non_consensus, N_POSITIONS))
  })

  airr_row <- function(sequence_id, sequence_alignment, germline_alignment_d_mask,
                       library_id = "GUSS_TEST", v_call = "IGHV1-8*01") {
    data.frame(
      sequence_id = sequence_id, library_id = library_id, locus = "IGH",
      sequence = gsub(".", "", sequence_alignment, fixed = TRUE), rev_comp = "F",
      v_call = v_call, d_call = "IGHD3-10*01", j_call = "IGHJ5*02",
      c_call = "IGHM*01", qualified_c_call = "IGHM*01", c_informative_positions = 300L,
      junction = "", junction_aa = "", cdr3 = "", cdr3_aa = "",
      sequence_alignment = sequence_alignment,
      germline_alignment = germline_alignment_d_mask,
      germline_alignment_d_mask = germline_alignment_d_mask,
      v_cigar = "", d_cigar = "", j_cigar = "",
      v_identity = 0.98, j_identity = 0.99, productive = "T", vj_in_frame = "T",
      stop_codon = "F", junction_length = 0L, consensus_count = 12,
      stringsAsFactors = FALSE
    )
  }

  vote_rows <- function(sequence_ids, total = VOTES$total,
                        non_consensus = VOTES$non_consensus) {
    data.frame(
      library_id = "GUSS_TEST", sequence_id = sequence_ids,
      total_votes_by_imgt_position = total,
      non_consensus_votes_by_imgt_position = non_consensus,
      stringsAsFactors = FALSE
    )
  }

  #' The value this stage used to release, computed the way it used to be.
  shazam_whole_v_mu_freq <- function(sequence_alignment, germline_alignment) {
    frame <- data.frame(
      sequence_id = "probe", sequence_alignment = sequence_alignment,
      germline_alignment_d_mask = germline_alignment, stringsAsFactors = FALSE
    )
    suppressWarnings(shazam::observedMutations(
      frame, sequenceColumn = "sequence_alignment",
      germlineColumn = "germline_alignment_d_mask",
      regionDefinition = shazam::IMGT_V, mutationDefinition = NULL,
      ambiguousMode = "eitherOr", frequency = TRUE, combine = TRUE, nproc = 1
    )$mu_freq)
  }

  region_values <- function(table, prefix) {
    as.integer(table[paste0(prefix, V_REGION_LEVELS)][1, ])
  }

  list(MIN_C_POSITIONS = MIN_C_POSITIONS,
       MIN_POSITION_VOTES = MIN_POSITION_VOTES,
       EXPECTED_COLUMNS = EXPECTED_COLUMNS,
       REGION_WIDTH = REGION_WIDTH,
       N_POSITIONS = N_POSITIONS,
       REGION_OF_POSITION = REGION_OF_POSITION,
       SERIALISATION_TOLERANCE = SERIALISATION_TOLERANCE,
       GERMLINE = GERMLINE,
       mutate_positions = mutate_positions,
       blank_positions = blank_positions,
       REGION_TOTAL_VOTES = REGION_TOTAL_VOTES,
       REGION_MAJORITY_VOTES = REGION_MAJORITY_VOTES,
       UNMAPPED_CDR2 = UNMAPPED_CDR2,
       SCAFFOLD_N_FWR2 = SCAFFOLD_N_FWR2,
       vote_cells = vote_cells,
       VOTES = VOTES,
       EMPTY_VOTES = EMPTY_VOTES,
       VOTED_WIDTH = VOTED_WIDTH,
       THIN_VOTES = THIN_VOTES,
       airr_row = airr_row,
       vote_rows = vote_rows,
       shazam_whole_v_mu_freq = shazam_whole_v_mu_freq,
       region_values = region_values)
}
