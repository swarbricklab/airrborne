# Check TRUST4 votes mapped into IMGT V positions, including gaps and orientation.
# Run: python test/run_tests.py --component shm --level unit

load_test_modules(c(
  "table_io", "gene_calls", "assembly_backbone", "trust4_consensus_votes"))

work_dir <- tempfile("trust4-votes-")
dir.create(work_dir)

WINDOW <- 12L

# --------------------------------------------------------------- fixtures ---

# Build a 4-by-position vote matrix in which `consensus_votes` sit on the base the
# consensus carries and `other_votes` sit on a different base, so a position's
# total and its non-consensus share can be set independently.
vote_matrix_for <- function(consensus, consensus_votes, other_votes = NULL) {
  bases <- strsplit(consensus, "", fixed = TRUE)[[1]]
  n <- length(bases)
  m <- matrix(0L, nrow = 4L, ncol = n, dimnames = list(VOTE_BASES, NULL))
  row_of <- match(bases, VOTE_BASES)
  for (j in seq_len(n)) {
    if (!is.na(row_of[j])) m[row_of[j], j] <- as.integer(consensus_votes[j])
    if (!is.null(other_votes) && other_votes[j] > 0L) {
      alternative <- if (is.na(row_of[j])) 1L else (row_of[j] %% 4L) + 1L
      m[alternative, j] <- as.integer(other_votes[j])
    }
  }
  m
}

# Write count lines exactly as TRUST4 does, with a trailing space before the newline.
count_lines_from_matrix <- function(m) {
  vapply(seq_len(nrow(m)),
         function(k) paste0(paste(m[k, ], collapse = " "), " "),
         character(1))
}

# A whole `_final.out` file. Each record is the six lines TRUST4 writes.
write_final_out <- function(path, records) {
  lines <- unlist(lapply(records, function(record) {
    c(paste0(">", record$id, " ", record$name %||% "IGHV1-8*01"),
      record$consensus,
      count_lines_from_matrix(record$matrix))
  }), use.names = FALSE)
  writeLines(lines, path)
}
`%||%` <- function(x, y) if (is.null(x)) y else x

# One QC-pass IGH record, carrying every column build_consensus_vote_rows needs
# plus the rest of the reannotated schema so the file reads as AIRR.
airr_row <- function(sequence_id, sequence, sequence_alignment, v_cigar,
                     v_sequence_start, v_sequence_end, junction = "",
                     rev_comp = "F", library_id = "GUSS_TEST") {
  data.frame(
    sequence_id = sequence_id, library_id = library_id, locus = "IGH",
    sequence = sequence, rev_comp = rev_comp,
    v_call = "IGHV1-8*01", d_call = "IGHD3-10*01", j_call = "IGHJ5*02",
    c_call = "IGHM", c_informative_positions = 300L,
    junction = junction, junction_aa = "", cdr3 = "", cdr3_aa = "",
    sequence_alignment = sequence_alignment, germline_alignment = "",
    germline_alignment_d_mask = "",
    v_cigar = v_cigar, d_cigar = "", j_cigar = "",
    v_sequence_start = v_sequence_start, v_sequence_end = v_sequence_end,
    v_germline_start = 1L, v_germline_end = 320L,
    v_identity = 0.98, j_identity = 0.99, productive = "T", vj_in_frame = "T",
    stop_codon = "F", junction_length = nchar(junction), consensus_count = 12,
    stringsAsFactors = FALSE
  )
}

# Read votes for one handcrafted record and map its IMGT V window.
map_one <- function(record, row, window = WINDOW) {
  totals <- position_vote_totals(record$matrix)
  consensus_votes <- position_consensus_votes(record$consensus, record$matrix)
  map_imgt_v_votes(
    row$sequence_alignment, row$sequence, row$v_cigar, row$v_sequence_start,
    row$v_sequence_end, record$consensus, totals, consensus_votes, window
  )
}

# --- Exact match, no gaps ----------------------------------------------------
plain <- list(id = "assemble0", consensus = "ACGTACGTACGT",
              matrix = vote_matrix_for("ACGTACGTACGT", 10:21))
plain_row <- airr_row("assemble0_0", "ACGTACGTACGT", "ACGTACGTACGT", "12=", 1L, 12L)
plain_map <- map_one(plain, plain_row)
test_that("an exact ungapped alignment maps every position", {
  expect_identical(plain_map$n_mapped, 12L)
  expect_identical(plain_map$n_mismatch, 0L)
  expect_identical(plain_map$total, 10:21)
  expect_identical(plain_map$non_consensus, rep(0L, 12L))
})

# --- Non-consensus votes are the dissenting share, not the winner ------------
dissent <- list(id = "assemble0", consensus = "ACGTACGTACGT",
                matrix = vote_matrix_for("ACGTACGTACGT", rep(10L, 12L),
                                         c(0L, 3L, rep(0L, 10L))))
dissent_map <- map_one(dissent, plain_row)
test_that("non-consensus votes are the dissenting share, not the winner", {
  expect_identical(dissent_map$total, c(10L, 13L, rep(10L, 10L)))
  expect_identical(dissent_map$non_consensus, c(0L, 3L, rep(0L, 10L)))
})

# --- The substituted CDR3 rewrote the last codon: null those, keep the rest ---
spliced_row <- airr_row("assemble0_3", "ACGTACGTAAAA", "ACGTACGTAAAA", "12=",
                        1L, 12L)
spliced_map <- map_one(plain, spliced_row)
test_that("a substituted CDR3 nulls the positions it rewrote", {
  expect_identical(spliced_map$n_mapped, 9L)
  expect_identical(spliced_map$n_mismatch, 3L)
  expect_identical(spliced_map$total, c(10:18, NA_integer_, NA_integer_, NA_integer_))
  expect_identical(spliced_map$non_consensus,
                   c(rep(0L, 9L), NA_integer_, NA_integer_, NA_integer_))
})

# --- A mis-anchored CDR3 window reaches further 5', and is recorded not refused --
# A substituted CDR3 window can extend into FWR3. Disagreeing positions lose
# their support mapping. The first disagreement identifies where that loss
# begins, including when it precedes the conserved Cys104 position.
early_row <- airr_row("assemble0_1", "ACGAACGTACGT", "ACGAACGTACGT", "12=", 1L, 12L)
early_map <- map_one(plain, early_row)
test_that("a mis-anchored CDR3 window is recorded rather than refused", {
  expect_identical(early_map$n_mismatch, 1L)
  expect_identical(early_map$first_mismatch, 4L)
  expect_identical(early_map$n_mapped, 11L)
  expect_true(is.na(early_map$total[4]))
  expect_identical(early_map$total[c(1:3, 5:12)], c(10L, 11L, 12L, 14:21))
})
test_that("the ordinary case reports its 5'-most mismatch inside the last codon", {
  expect_identical(spliced_map$first_mismatch, 10L)
  expect_true(is.na(plain_map$first_mismatch))
})

# --- I skips inserted query bases from the alignment -------------------------
# 4=3I5= consumes 12 query bases but only maps 9 of them; the alignment carries the
# 9 it kept, with the three dropped positions showing as gap columns.
dropped <- list(id = "assemble1", consensus = "ACGTTTTACGTA",
                matrix = vote_matrix_for("ACGTTTTACGTA", 21:32))
dropped_row <- airr_row("assemble1_0", "ACGTTTTACGTA", "ACGT...ACGTA", "4=3I5=",
                        1L, 12L)
dropped_map <- map_one(dropped, dropped_row)
test_that("I skips inserted query bases from the alignment", {
  expect_identical(dropped_map$n_mapped, 9L)
  expect_identical(dropped_map$total,
                   c(21L, 22L, 23L, 24L, NA_integer_, NA_integer_, NA_integer_,
                     28L, 29L, 30L, 31L, 32L))
})

# --- D is germline-only, so it consumes no query base ------------------------
inserted <- list(id = "assemble2", consensus = "ACGTACGTACGT",
                 matrix = vote_matrix_for("ACGTACGTACGT", 41:52))
inserted_row <- airr_row("assemble2_0", "ACGTACGTACGT", "ACGT--ACGTAC",
                         "4=2D6=", 1L, 10L)
inserted_map <- map_one(inserted, inserted_row)
test_that("D is germline-only, so it consumes no query base", {
  expect_identical(inserted_map$n_mapped, 10L)
  expect_identical(inserted_map$total,
                   c(41L, 42L, 43L, 44L, NA_integer_, NA_integer_,
                     45L, 46L, 47L, 48L, 49L, 50L))
  # Treating D as a query-consuming operation shifts every base after it, which the
  # per-base identity check is there to catch.
  expect_error(
    map_cigar_to_sequence_index("4=2D6=", 1L, 12L),
    "consumes 10 query base(s) but v_sequence_end is 12", fixed = TRUE
  )
})

# --- A leading soft clip starts the alignment inside the sequence ------------
clipped <- list(id = "assemble3", consensus = "TTTACGTACGTA",
                matrix = vote_matrix_for("TTTACGTACGTA", 61:72))
clipped_row <- airr_row("assemble3_0", "TTTACGTACGTA", "ACGTACGTA", "3S9=",
                        4L, 12L)
clipped_map <- map_one(clipped, clipped_row)
test_that("a leading soft clip starts the alignment inside the sequence", {
  expect_identical(clipped_map$n_mapped, 9L)
  expect_identical(clipped_map$total, c(64:72, NA_integer_, NA_integer_, NA_integer_))
  expect_error(
    map_cigar_to_sequence_index("3S9=", 1L, 12L),
    "starts the V alignment at sequence position 4 but v_sequence_start is 1",
    fixed = TRUE
  )
})

# --- A mate-pair scaffold gap has no called base, so dissent is undefined -----
# TRUST4 writes N where it joined two mate pairs across a gap it could not read.
# Reads are still tallied at some of those positions -- position 4 below carries
# five votes -- but with no consensus base to agree with, the split into agreeing
# and dissenting votes does not exist. The total is reported and the non-consensus
# share is null. Reporting zero would make every one of those votes read as a vote
# against the consensus.
gapped <- list(id = "assemble4", consensus = "ACGNNNNNNNGT",
               matrix = vote_matrix_for("ACGNNNNNNNGT",
                                        c(80L, 81L, 82L, rep(0L, 7L), 83L, 84L),
                                        other_votes = c(0L, 0L, 0L, 5L,
                                                        rep(0L, 6L), 0L, 0L)))
gapped_row <- airr_row("assemble4_0", "ACGNNNNNNNGT", "ACGNNNNNNNGT", "12=",
                       1L, 12L)
gapped_map <- map_one(gapped, gapped_row)
test_that("a scaffold gap reports its total but no dissenting share", {
  expect_identical(gapped_map$n_mapped, 12L)
  expect_identical(gapped_map$total, c(80L, 81L, 82L, 5L, rep(0L, 6L), 83L, 84L))
  expect_identical(gapped_map$non_consensus,
                   c(0L, 0L, 0L, rep(NA_integer_, 7L), 0L, 0L))
})

# --- A pseudo-count position is reported as one vote, not filtered -----------
pseudo <- list(id = "assemble5", consensus = "ACGTACGTACGT",
               matrix = vote_matrix_for("ACGTACGTACGT",
                                        c(500L, 1L, rep(500L, 10L))))
pseudo_map <- map_one(pseudo, airr_row("assemble5_0", "ACGTACGTACGT",
                                       "ACGTACGTACGT", "12=", 1L, 12L))
test_that("a pseudo-count position is reported as one vote, not filtered", {
  expect_identical(pseudo_map$total[2], 1L)
})

# --- The supporting count comes from the emitted base, never from which.max ---
# Real files never disagree, but the contract is the count of the base TRUST4
# emitted, so a tie or a future change cannot swap in another base's count.
skewed <- matrix(0L, nrow = 4L, ncol = 1L, dimnames = list(VOTE_BASES, NULL))
skewed["A", 1] <- 9L
skewed["T", 1] <- 2L
test_that("the supporting count comes from the emitted base", {
  expect_identical(position_consensus_votes("T", skewed), 2L)
  expect_identical(position_vote_totals(skewed), 11L)
})

# --- A V alignment that ends inside the window leaves the tail unmapped ------
short <- list(id = "assemble6", consensus = "ACGTACGTACGT",
              matrix = vote_matrix_for("ACGTACGTACGT", 91:102))
short_row <- airr_row("assemble6_0", "ACGTACGTACGT", "ACGTACGTACGT", "6=", 1L, 6L)
short_map <- map_one(short, short_row)
test_that("a V alignment ending inside the window leaves the tail unmapped", {
  expect_identical(short_map$n_mapped, 6L)
  expect_identical(short_map$total, c(91:96, rep(NA_integer_, 6L)))
})

# --- Per-base identity is the check that makes the map trustworthy -----------
test_that("per-base identity is the check that makes the map trustworthy", {
  expect_error(
    map_alignment_columns_to_sequence("ACGTACGTACGT", "ACGTACGTACGA", 0:11, WINDOW),
    "disagree at 1 mapped base(s); first at alignment column 12", fixed = TRUE
  )
  expect_error(
    map_alignment_columns_to_sequence("ACGT", "ACGTACGTACGT", 0:11, WINDOW),
    "holds 4 column(s) with an observed base but the V alignment maps 12",
    fixed = TRUE
  )
})

# --- Cigar contract ---------------------------------------------------------
test_that("the cigar contract is enforced", {
  expect_error(map_cigar_to_sequence_index("10H2=", 1L, 12L),
               "unsupported operation(s) H", fixed = TRUE)
  expect_error(map_cigar_to_sequence_index("not-a-cigar", 1L, 12L),
               "does not parse as (length, operation) pairs", fixed = TRUE)
  expect_error(map_cigar_to_sequence_index("", 1L, 12L),
               "must be one non-empty string", fixed = TRUE)
  expect_identical(map_cigar_to_sequence_index("2=1X3=", 1L, 6L), 0:5)
})

# --- final.out format contract ----------------------------------------------
test_that("the final.out format contract is enforced", {
  expect_identical(names(parse_final_out_records(
    c(">assemble7 IGHV1-8*01", "AC", "1 0 ", "0 2 ", "0 0 ", "0 0 "),
    "assemble7")), "assemble7")
  # Records that were not asked for are dropped without being parsed.
  expect_identical(length(parse_final_out_records(
    c(">assemble7 IGHV1-8*01", "AC", "1 0 ", "0 2 ", "0 0 ", "0 0 "),
    "assemble9")), 0L)
  expect_error(
    parse_final_out_records(c(">assemble7 IGHV1-8*01", "AC", "1 0 "), "assemble7"),
    "not a whole number of 6-line records", fixed = TRUE
  )
  expect_error(
    parse_final_out_records(
      c(">ACGTACGT_7 IGHV1-8*01", "AC", "1 0 ", "0 2 ", "0 0 ", "0 0 "),
      "ACGTACGT_7"),
    "record header is not '>assemble<N>'", fixed = TRUE
  )
  expect_error(parse_vote_matrix(c("1 0 ", "0 2 ", "0 0 "), 2L),
               "has 3 count line(s); expected 4", fixed = TRUE)
  expect_error(parse_vote_matrix(c("1 ", "0 2 ", "0 0 ", "0 0 "), 2L),
               "holds 7 vote(s) across 4 count lines; expected 8", fixed = TRUE)
})

# --- Junction votes are located in the contig, not in the record ------------
junction_totals <- position_vote_totals(plain$matrix)
junction_consensus <- position_consensus_votes(plain$consensus, plain$matrix)
unique_junction <- map_junction_votes("CGTACGTACG", plain$consensus,
                                      junction_totals, junction_consensus)
test_that("a junction is located in the contig, not in the record", {
  expect_identical(unique_junction$status, "mapped")
  expect_identical(unique_junction$start, 2L)
  expect_identical(unique_junction$total, "11,12,13,14,15,16,17,18,19,20")
})
test_that("a junction that cannot be located reports why", {
  expect_identical(map_junction_votes("AC", plain$consensus, junction_totals,
                                      junction_consensus)$status,
                   "ambiguous_in_consensus")
  expect_identical(map_junction_votes("TTTTTTTT", plain$consensus, junction_totals,
                                      junction_consensus)$status,
                   "absent_from_consensus")
  expect_identical(map_junction_votes("", plain$consensus, junction_totals,
                                      junction_consensus)$status,
                   "junction_missing")
  expect_identical(map_junction_votes(NA, plain$consensus, junction_totals,
                                      junction_consensus)$status,
                   "junction_missing")
})

# --- Vote vectors serialise to a fixed token count, nulls included ----------
# `parse_vote_vector` is the reader paired with the writer: it restores the nulls
# and, unlike `strsplit`, keeps trailing tokens, so a vector that ends in nulls
# still has 312 values. A null is not a zero -- it means no coordinate mapped
# there, while zero is real evidence that no read voted -- so the round trip has
# to preserve both distinctly.
round_trip <- format_vote_vector(c(1L, NA_integer_, 3L), 3L)
test_that("vote vectors serialise to a fixed token count, nulls included", {
  expect_identical(round_trip, "1,NA,3")
  expect_identical(parse_vote_vector(round_trip, 3L), c(1L, NA_integer_, 3L))
  expect_identical(parse_vote_vector("0,NA,0", 3L), c(0L, NA_integer_, 0L))
  expect_identical(length(parse_vote_vector(
    format_vote_vector(rep(NA_integer_, IMGT_V_POSITIONS), IMGT_V_POSITIONS),
    IMGT_V_POSITIONS)),
    IMGT_V_POSITIONS)
  # An absent cell is a record with no mapped position, not a parse failure.
  expect_identical(parse_vote_vector(NA_character_, 3L), rep(NA_integer_, 3L))
  expect_identical(parse_vote_vector("", 3L), rep(NA_integer_, 3L))
  expect_error(format_vote_vector(c(1L, 2L), 3L),
               "has 2 value(s); expected 3", fixed = TRUE)
  expect_error(parse_vote_vector("1,2", 3L),
               "has 2 token(s); expected 3", fixed = TRUE)
})

# --- Majority votes: the reads that agreed with the assembled base ----------
#
# This is the quantity a support threshold has to be applied to, so the two ways
# it can be wrong are both asserted: silently turning an unmapped position into a
# zero, and accepting two vectors that cannot describe the same record.
test_that("majority votes are the total less the dissenting share", {
  expect_identical(position_majority_votes(c(10L, 4L, 0L), c(0L, 1L, 0L)),
                   c(10L, 3L, 0L))
})

test_that("an unmapped position stays unmapped rather than becoming zero", {
  expect_identical(
    position_majority_votes(c(10L, NA_integer_, 5L), c(1L, NA_integer_, NA_integer_)),
    c(9L, NA_integer_, NA_integer_)
  )
})

test_that("majority votes cannot be negative or come from a different record", {
  expect_error(position_majority_votes(c(2L, 3L), c(1L, 9L)),
               "non-consensus votes exceed total votes at 1 position(s)",
               fixed = TRUE)
  expect_error(position_majority_votes(c(1L, 2L, 3L), c(0L, 0L)),
               "total votes have 3 position(s) and non-consensus votes 2",
               fixed = TRUE)
})

test_that("an empty window gives an empty result", {
  expect_identical(position_majority_votes(integer(0), integer(0)), integer(0))
})

# --- The backbone id is the final.out record id -----------------------------
test_that("the backbone id is the final.out record id", {
  expect_identical(trust4_assembly_backbone_id(c("assemble0_0", "assemble12_7")),
                   c("assemble0", "assemble12"))
  expect_identical(trust4_assembly_backbone_id("GUSS_A~assemble3_11"),
                   "GUSS_A~assemble3")
})

# --- Reading a file: votes must come from the same assembly run -------------
final_out_file <- file.path(work_dir, "votes_final.out")
write_final_out(final_out_file, list(plain, dropped))
read_votes <- read_backbone_votes(final_out_file, c("assemble0", "assemble1"))
test_that("votes must come from the same assembly run", {
  expect_identical(sort(names(read_votes)), c("assemble0", "assemble1"))
  expect_identical(read_votes$assemble0$total, 10:21)
  expect_identical(read_votes$assemble0$consensus, "ACGTACGTACGT")
  expect_error(read_backbone_votes(final_out_file, "assemble999"),
               "missing 1 assembly backbone(s) referenced by the AIRR table",
               fixed = TRUE)
  # A chunk size below the record count must give the same answer.
  expect_identical(
    read_backbone_votes(final_out_file, "assemble1", chunk_records = 1L),
    read_backbone_votes(final_out_file, "assemble1")
  )
})

# --- A record on the opposite strand to its contig is emptied, not mapped -----
# IgBLAST sometimes finds the receptor on the other strand and Change-O reports a
# reverse complemented `sequence`. The alignment is reverse complemented with it, so
# it still agrees with `sequence` base for base and no alignment check can see the
# problem -- but the record shares no position with its contig, so nothing can be
# attributed to it. It is emptied and the reason is recorded.
reversed_contig <- list(id = "assemble7", consensus = "AACCGGTTACGC",
                        matrix = vote_matrix_for("AACCGGTTACGC", 30:41))
reversed_final_out <- file.path(work_dir, "reversed_final.out")
write_final_out(reversed_final_out, list(reversed_contig))
reversed_votes <- read_backbone_votes(reversed_final_out, "assemble7")
flipped <- .reverse_complement("AACCGGTTACGC")
reversed_rows <- build_consensus_vote_rows(
  airr_row("assemble7_0", flipped, flipped, "12=", 1L, 12L,
           junction = substr(flipped, 1L, 6L), rev_comp = "T"),
  reversed_votes, WINDOW)
test_that("a record on the opposite strand to its contig is emptied", {
  expect_identical(nrow(reversed_rows), 1L)
  expect_identical(reversed_rows$imgt_v_positions_mapped, 0L)
  expect_identical(reversed_rows$junction_vote_status, "opposite_contig_orientation")
  # Both vectors keep their token count, so the reader contract is unchanged.
  expect_identical(reversed_rows$total_votes_by_imgt_position,
                   format_vote_vector(rep(NA_integer_, WINDOW), WINDOW))
  expect_identical(reversed_rows$non_consensus_votes_by_imgt_position,
                   format_vote_vector(rep(NA_integer_, WINDOW), WINDOW))
  # No attribution was attempted, so there is no per-position mismatch count.
  expect_true(is.na(reversed_rows$imgt_v_positions_consensus_mismatch))
  expect_true(is.na(reversed_rows$junction_start_in_consensus))
})
test_that("orientation is decided by comparison, not by a cutoff", {
  # The criterion is a comparison, not a cutoff: whichever orientation matches the
  # contig better wins. A record differing from its contig only over a substituted
  # CDR3 stays forward even though it is far from a perfect match.
  expect_false(record_is_forward_of_contig(flipped, "AACCGGTTACGC"))
  expect_true(record_is_forward_of_contig("AACCGGTTACGC", "AACCGGTTACGC"))
  expect_true(record_is_forward_of_contig("AACCGGTTTTTT", "AACCGGTTACGC"))

  # Why it cannot be a cutoff. Some contigs are close to their own reverse
  # complement, and a record on the opposite strand of one of those still matches
  # it forwards very well. This contig is 10 of 12 self-reverse-complementary, so
  # the reversed record below matches forwards at 0.83 while matching its own
  # reverse complement perfectly. Any forward cutoff loose enough to admit an
  # ordinary record would call this one forward and hand it another strand's votes.
  self_similar_contig <- "AACCGGCCGGTA"
  reversed_record <- .reverse_complement(self_similar_contig)
  expect_identical(.fraction_matching(reversed_record, self_similar_contig), 10 / 12)
  expect_identical(
    .fraction_matching(.reverse_complement(reversed_record), self_similar_contig), 1)
  expect_false(record_is_forward_of_contig(reversed_record, self_similar_contig))

  # A perfect palindrome matches equally either way round. The comparison keeps
  # such a record forward, which costs nothing: both orientations name the same
  # bases, so the votes attributed to it are the right ones.
  palindrome <- "AACCGGCCGGTT"
  expect_identical(.reverse_complement(palindrome), palindrome)
  expect_true(record_is_forward_of_contig(palindrome, palindrome))
  # An alignment in the opposite frame to its own `sequence` is a different failure,
  # and that one is fatal because it means the coordinate map is wrong.
  expect_error(
    map_alignment_columns_to_sequence("ACGTACGTACGT", "ACGTACGTACGT",
                                      rev(0:11), WINDOW),
    "disagree at", fixed = TRUE
  )
})

# --- Sibling records of one contig share its V votes ------------------------
siblings <- rbind(plain_row, spliced_row)
sibling_rows <- build_consensus_vote_rows(siblings, read_votes, WINDOW)
test_that("sibling records of one contig share its V votes", {
  expect_identical(sibling_rows$sequence_id, c("assemble0_0", "assemble0_3"))
  expect_identical(sibling_rows$trust4_backbone_id, c("assemble0", "assemble0"))
  expect_identical(sibling_rows$consensus_length, c(12L, 12L))
  expect_identical(sibling_rows$imgt_v_positions_mapped, c(12L, 9L))
  expect_identical(sibling_rows$imgt_v_positions_consensus_mismatch, c(0L, 3L))
  # Both carry the same votes over the positions the substitution cannot reach.
  expect_identical(
    strsplit(sibling_rows$total_votes_by_imgt_position[1], ",",
             fixed = TRUE)[[1]][1:9],
    strsplit(sibling_rows$total_votes_by_imgt_position[2], ",",
             fixed = TRUE)[[1]][1:9]
  )
})
