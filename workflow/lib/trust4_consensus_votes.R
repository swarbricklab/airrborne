# Read TRUST4 A/C/G/T vote tallies and map contig positions into IMGT V coordinates.
# Unmapped positions remain missing. TRUST4 replaces a zero tally with one vote on
# the consensus base, so one vote has ambiguous support (ExtendSeqFromReads in
# SeqSet.hpp at v1.1.10:
# https://github.com/liulab-dfci/TRUST4/blob/v1.1.10/SeqSet.hpp#L10781-L10784).
# These helpers use base R.

FINAL_OUT_RECORD_LINES <- 6L

# Row order of the four count lines.
VOTE_BASES <- c("A", "C", "G", "T")

# Both characters IMGT-aligned AIRR strings use for a column with no observed base:
# `.` for numbering spacers and non-coverage, `-` for a deletion.
ALIGNMENT_GAP_CHARACTERS <- c(".", "-")

# Length of the IMGT V region SHM is measured over, which is `shazam::IMGT_V@seqLength`.
# Using the same window keeps these vectors index-for-index with `shm_frequency`.
IMGT_V_POSITIONS <- 312L

# IMGT unique numbering places the conserved Cys104 codon at V positions 310-312.
# TRUST4 normally starts overwriting the consensus there when it substitutes a
# CDR3. A mismatch before this position signals that the substitution extends
# into the V region. The distinction is recorded, not enforced; see
# `map_imgt_v_votes()`.
IMGT_CYS104_FIRST_POSITION <- 310L

# --------------------------------------------------------- reading final.out ---

#' Split `_final.out` text into the records that were asked for.
#'
#' `lines` must begin at a record boundary and hold a whole number of six-line
#' records; a chunked reader that takes a multiple of six lines at a time
#' satisfies both. `wanted_backbone_ids` are contig ids of the form `assemble<N>`.
#' Returns a named list, one element per requested id present in `lines`, each
#' holding the consensus sequence and its four raw count lines. Records that were
#' not asked for are discarded without being parsed, which is what keeps a
#' cohort-wide pass affordable.
parse_final_out_records <- function(lines, wanted_backbone_ids) {
  n <- length(lines)
  if (n %% FINAL_OUT_RECORD_LINES != 0L) {
    stop("TRUST4 final.out is not a whole number of ", FINAL_OUT_RECORD_LINES,
         "-line records: read ", n, " line(s)")
  }
  if (!n) return(list())

  header_at <- seq.int(1L, n, by = FINAL_OUT_RECORD_LINES)
  headers <- lines[header_at]
  # A barcoded run names records `><barcode>_<i>`, which would silently fail to
  # match the ids in the AIRR table, so reject anything but the bulk spelling.
  malformed <- !grepl("^>assemble[0-9]+( |$)", headers)
  if (any(malformed)) {
    stop("TRUST4 final.out record header is not '>assemble<N>': ",
         headers[which(malformed)[1]])
  }

  ids <- sub(" .*$", "", substring(headers, 2L))
  keep <- which(ids %in% as.character(wanted_backbone_ids))
  records <- lapply(keep, function(k) {
    first <- header_at[k]
    list(consensus = lines[first + 1L],
         count_lines = lines[(first + 2L):(first + 5L)])
  })
  names(records) <- ids[keep]
  records
}

#' Turn one record's four count lines into a base-by-position vote matrix.
#'
#' Returns a 4 by `consensus_length` integer matrix whose rows are A, C, G, T in
#' that order and whose column j holds the votes at sequence character j. Errors
#' unless there are exactly four count lines carrying exactly `consensus_length`
#' integers each, because a short line would shift every position after it.
parse_vote_matrix <- function(count_lines, consensus_length) {
  if (length(count_lines) != length(VOTE_BASES)) {
    stop("TRUST4 final.out record has ", length(count_lines),
         " count line(s); expected ", length(VOTE_BASES))
  }
  values <- scan(text = paste(count_lines, collapse = " "), what = integer(),
                 quiet = TRUE)
  expected <- length(VOTE_BASES) * consensus_length
  if (length(values) != expected) {
    stop("TRUST4 final.out record holds ", length(values), " vote(s) across ",
         length(VOTE_BASES), " count lines; expected ", expected, " for a ",
         consensus_length, "-base consensus")
  }
  matrix(values, nrow = length(VOTE_BASES), byrow = TRUE,
         dimnames = list(VOTE_BASES, NULL))
}

#' Total votes cast at each position of a consensus.
position_vote_totals <- function(vote_matrix) {
  as.integer(colSums(vote_matrix))
}

#' Votes for the base TRUST4 actually wrote at each position of a consensus.
#'
#' Deliberately not `which.max`: ties keep the incumbent base, so the written
#' base is a maximum but need not be the arithmetic argmax.
#'
#' A position whose written base is not A, C, G or T has no called base: TRUST4
#' writes N where it scaffolded two mate pairs together across a gap it could not
#' read. Reads are still sometimes tallied at such a position, but with no base to
#' agree with, the split of those votes into agreeing and dissenting does not
#' exist, so this returns NA there rather than zero. Zero would make every vote
#' cast at a scaffold gap read as a vote against the consensus.
position_consensus_votes <- function(consensus, vote_matrix) {
  bases <- strsplit(as.character(consensus), "", fixed = TRUE)[[1]]
  row_of <- match(bases, VOTE_BASES)
  votes <- rep(NA_integer_, length(bases))
  known <- !is.na(row_of)
  if (any(known)) {
    votes[known] <- vote_matrix[cbind(row_of[known], which(known))]
  }
  votes
}

#' Read the consensus and per-position votes of the requested contigs.
#'
#' Streams `final_out_path` in `chunk_records` records at a time so the whole file
#' never has to be held in memory, keeping only the requested contigs. Returns a
#' named list, one element per requested id, each with the consensus sequence and
#' its per-position total and consensus-base vote counts. Errors if any requested
#' id is absent, which is what stops votes from one assembly run being attached to
#' an AIRR table produced from a different one.
#'
#' Side effects: opens and closes the file, and reports the record count read.
read_backbone_votes <- function(final_out_path, wanted_backbone_ids,
                                chunk_records = 5000L) {
  wanted <- unique(as.character(wanted_backbone_ids))
  if (!length(wanted)) return(list())

  connection <- file(final_out_path, "r")
  on.exit(close(connection), add = TRUE)

  votes <- list()
  records_read <- 0L
  repeat {
    lines <- readLines(connection,
                       n = FINAL_OUT_RECORD_LINES * chunk_records, warn = FALSE)
    if (!length(lines)) break
    records_read <- records_read + length(lines) %/% FINAL_OUT_RECORD_LINES
    records <- parse_final_out_records(lines, wanted)
    for (id in names(records)) {
      record <- records[[id]]
      matrix_of_votes <- parse_vote_matrix(record$count_lines,
                                           nchar(record$consensus))
      votes[[id]] <- list(
        consensus = record$consensus,
        total = position_vote_totals(matrix_of_votes),
        consensus_votes = position_consensus_votes(record$consensus,
                                                   matrix_of_votes)
      )
    }
  }

  absent <- setdiff(wanted, names(votes))
  if (length(absent)) {
    stop("TRUST4 final.out is missing ", length(absent),
         " assembly backbone(s) referenced by the AIRR table, starting with ",
         absent[1], "; the votes and the reannotated table must come from the ",
         "same assembly run")
  }
  message("Read ", records_read, " final.out record(s) from ", final_out_path,
          ", keeping ", length(votes), " referenced backbone(s)")
  votes
}

# ------------------------------------------------------------ coordinate map ---

#' Sequence positions the V alignment consumes, in alignment order.
#'
#' Returns the 0-based indices into `sequence` of the V alignment's bases, in
#' order, so element k is where the k-th ungapped alignment base came from.
#'
#' CIGARs use AIRR semantics: I advances the query and D does not. Reannotation
#' normalises the pinned Change-O convention before writing this table. Inserted
#' query bases are absent from the IMGT substitution alignment. Coordinate checks
#' stop inconsistent inputs before an indel can shift the mapped votes.
map_cigar_to_sequence_index <- function(v_cigar, v_sequence_start, v_sequence_end) {
  cigar <- trimws(as.character(v_cigar))
  if (length(cigar) != 1L || is.na(cigar) || !nzchar(cigar)) {
    stop("v_cigar must be one non-empty string")
  }
  pairs <- regmatches(cigar, gregexpr("[0-9]+[A-Z=]", cigar))[[1]]
  if (!identical(paste(pairs, collapse = ""), cigar)) {
    stop("v_cigar does not parse as (length, operation) pairs: ", cigar)
  }
  run_lengths <- as.integer(sub("[A-Z=]$", "", pairs))
  operations <- sub("^[0-9]+", "", pairs)
  unsupported <- setdiff(operations, c("S", "N", "M", "=", "X", "I", "D"))
  if (length(unsupported)) {
    stop("v_cigar contains unsupported operation(s) ",
         paste(unsupported, collapse = ", "), ": ", cigar)
  }

  consumed <- vector("list", length(operations))
  query <- 0L
  for (i in seq_along(operations)) {
    run <- run_lengths[i]
    if (operations[i] %in% c("S", "I")) {
      # Query bases the alignment does not represent: the soft-clipped prefix, and
      # the bases Change-O drops to keep IMGT numbering.
      query <- query + run
    } else if (operations[i] %in% c("N", "D")) {
      # Germline positions with no query base; they carry no sequence coordinate.
    } else {
      consumed[[i]] <- query + seq_len(run) - 1L
      query <- query + run
    }
  }
  mapped <- unlist(consumed, use.names = FALSE)
  if (is.null(mapped)) mapped <- integer(0)

  start <- as.integer(v_sequence_start)
  end <- as.integer(v_sequence_end)
  if (!length(mapped)) stop("v_cigar aligns no query bases: ", cigar)
  if (mapped[1L] != start - 1L) {
    stop("v_cigar starts the V alignment at sequence position ", mapped[1L] + 1L,
         " but v_sequence_start is ", start)
  }
  if (query != end) {
    stop("v_cigar consumes ", query, " query base(s) but v_sequence_end is ", end)
  }
  mapped
}

#' Where each alignment column of the window sits in `sequence`, verified base by base.
#'
#' Returns one 0-based `sequence` index per alignment column 1..`n_positions`, or
#' NA where the column is a gap or falls beyond the V alignment.
#'
#' Every base the V alignment maps is checked against `sequence` before any index
#' is returned, over the whole V alignment rather than just the window. That check
#' is what makes the map trustworthy: the alternative formulation, which advances an
#' alignment-column cursor on the germline-only operations, returns output of the
#' right type, length and plausible magnitude while being wrong on every record,
#' and only per-base identity catches it.
#'
#' It fixes the frame between the alignment and `sequence` and nothing more. It
#' cannot see a record whose `sequence` is the reverse complement of its contig,
#' because the alignment is reverse complemented with it and the two still agree;
#' `record_is_forward_of_contig()` below is what catches that.
map_alignment_columns_to_sequence <- function(sequence_alignment, sequence,
                                              query_indices, n_positions) {
  characters <- strsplit(as.character(sequence_alignment), "", fixed = TRUE)[[1]]
  observed <- !(characters %in% ALIGNMENT_GAP_CHARACTERS)
  n_query <- length(query_indices)
  if (sum(observed) < n_query) {
    stop("sequence_alignment holds ", sum(observed),
         " column(s) with an observed base but the V alignment maps ", n_query,
         " query base(s)")
  }

  observed_columns <- which(observed)[seq_len(n_query)]
  disagree <- characters[observed_columns] !=
    substring(sequence, query_indices + 1L, query_indices + 1L)
  if (any(disagree)) {
    first <- which(disagree)[1L]
    stop("sequence_alignment and sequence disagree at ", sum(disagree),
         " mapped base(s); first at alignment column ", observed_columns[first],
         " (sequence position ", query_indices[first] + 1L, ")")
  }

  # Rank each window column among the columns holding an observed base, then read
  # off the sequence position the V alignment put there. A column past the end of
  # the alignment, or holding a gap, has no rank and so stays unmapped.
  ranks <- rep(NA_integer_, n_positions)
  take <- min(length(characters), n_positions)
  if (take > 0L) {
    in_window <- observed[seq_len(take)]
    ranks[seq_len(take)][in_window] <- cumsum(in_window)[in_window]
  }
  mapped <- rep(NA_integer_, n_positions)
  within_v <- !is.na(ranks) & ranks <= n_query
  mapped[within_v] <- query_indices[ranks[within_v]]
  mapped
}

#' Reverse complement of a nucleotide string, leaving other characters as they are.
.reverse_complement <- function(sequence) {
  bases <- rev(strsplit(as.character(sequence), "", fixed = TRUE)[[1]])
  paste(chartr("ACGT", "TGCA", bases), collapse = "")
}

#' Fraction of positions at which two sequences carry the same character.
.fraction_matching <- function(left, right) {
  left_bases <- strsplit(as.character(left), "", fixed = TRUE)[[1]]
  right_bases <- strsplit(as.character(right), "", fixed = TRUE)[[1]]
  shared <- min(length(left_bases), length(right_bases))
  if (!shared) return(0)
  sum(left_bases[seq_len(shared)] == right_bases[seq_len(shared)]) / shared
}

#' Test whether a record's sequence agrees with its contig in the forward direction.
#'
#' Every vote is indexed by contig position, so attributing one to a record assumes
#' the record's `sequence` is that contig read in the same direction. Usually it is
#' the contig verbatim, give or take the substituted CDR3. But IgBLAST sometimes
#' finds the receptor on the opposite strand and Change-O then reports a reverse
#' complemented `sequence`. Such a record agrees with its own `sequence_alignment`
#' base for base, so no alignment check can see it, yet it shares no position with
#' its contig, and the positions that happen to carry the same base by chance would
#' be given votes describing unrelated bases.
#'
#' Compare forward and reverse-complement agreement directly. A fixed forward
#' identity cutoff cannot reliably distinguish near-palindromic contigs.
record_is_forward_of_contig <- function(sequence, consensus) {
  .fraction_matching(sequence, consensus) >=
    .fraction_matching(.reverse_complement(sequence), consensus)
}

# ------------------------------------------------------------- vote vectors ---

#' Serialise a vote vector as one comma-separated cell.
#'
#' Missing values are written as the literal token `NA` rather than left empty,
#' because `strsplit` drops empty trailing tokens and the reader must be able to
#' rely on the token count. Read a cell back with `parse_vote_vector()`.
format_vote_vector <- function(values, n_positions) {
  if (length(values) != n_positions) {
    stop("vote vector has ", length(values), " value(s); expected ", n_positions)
  }
  paste(ifelse(is.na(values), "NA", as.character(as.integer(values))),
        collapse = ",")
}

#' Read one serialised vote vector back into an integer vector.
#'
#' The inverse of `format_vote_vector()`. `cell` is one comma-separated cell from
#' a released vote column and `n_positions` is the token count it must carry.
#' Returns an integer vector of that length; an empty or missing cell gives all
#' `NA`, which is how a record with no mapped position is stored.
#'
#' `NA` and `0` mean different things here and the distinction is the point of
#' having one reader: `NA` is "no number belongs at this position", while `0` is
#' real evidence that no read voted. A total is `NA` when no coordinate mapped
#' there: an IMGT numbering spacer, a truncated end of the V region, a base the
#' substituted CDR3 rewrote, or a record read in the opposite direction to its
#' contig. A non-consensus count is also `NA` where the contig carries no
#' called base (an `N` at a mate-pair scaffold join), because votes cannot be split
#' into agreeing and dissenting without a base to agree with. `scan` restores
#' the literal `NA` token without the coercion warning `as.integer` raises, and
#' keeps trailing tokens that `strsplit` would drop.
parse_vote_vector <- function(cell, n_positions) {
  value <- as.character(cell)
  if (length(value) != 1L) {
    stop("parse_vote_vector expects one cell; got ", length(value))
  }
  if (is.na(value) || !nzchar(value)) return(rep(NA_integer_, n_positions))
  values <- scan(text = value, sep = ",", what = integer(), quiet = TRUE)
  if (length(values) != n_positions) {
    stop("vote vector has ", length(values), " token(s); expected ", n_positions)
  }
  values
}

#' Reads supporting the base a position was actually called from.
#'
#' `total` and `non_consensus` are two parsed vote vectors of equal length, as
#' released in `total_votes_by_imgt_position` and
#' `non_consensus_votes_by_imgt_position`. Returns their difference: the votes
#' that agreed with the contig's base. `NA` in either input gives `NA`, so a
#' position with no mapped coordinate stays unmeasured rather than becoming a
#' zero.
#'
#' This is the quantity to threshold on when asking how well a position is
#' supported, because a position can carry many votes and still be split. The
#' pipeline's per-sub-region vote summary takes it from here, so "supported" is
#' defined in one place. A negative difference means the two vectors do not
#' describe the same record and is an error rather than a value.
position_majority_votes <- function(total, non_consensus) {
  if (length(total) != length(non_consensus)) {
    stop("total votes have ", length(total), " position(s) and non-consensus ",
         "votes ", length(non_consensus))
  }
  majority <- as.integer(total) - as.integer(non_consensus)
  if (any(!is.na(majority) & majority < 0L)) {
    stop("non-consensus votes exceed total votes at ",
         sum(!is.na(majority) & majority < 0L), " position(s)")
  }
  majority
}

#' Votes behind one record's IMGT V positions.
#'
#' Returns the per-position total and non-consensus vote counts over alignment
#' columns 1..`n_positions`, how many positions carry a count, how many were
#' dropped because the record's base disagreed with its contig, and the 5'-most
#' position that disagreed. A position is NA when the alignment does not cover it:
#' an IMGT numbering spacer, a 5'-truncated V region, a V alignment ending inside
#' the window, or a base the substituted CDR3 rewrote (see
#' `IMGT_CYS104_FIRST_POSITION`).
#'
#' The 5'-most disagreement is returned because it is what separates the ordinary
#' case from a damaged one. Disagreement at or after `IMGT_CYS104_FIRST_POSITION`
#' is expected, being where the substituted CDR3 starts. A record that starts
#' disagreeing well before it has part of its V region replaced by a read-level
#' CDR3 variant, and its surviving positions carry that caveat. Either way this is
#' a property of the assembly, not of the coordinate map, which
#' `map_alignment_columns_to_sequence()` has already verified base by base.
map_imgt_v_votes <- function(sequence_alignment, sequence, v_cigar,
                             v_sequence_start, v_sequence_end,
                             consensus, total_votes, consensus_votes,
                             n_positions) {
  query_indices <- map_cigar_to_sequence_index(v_cigar, v_sequence_start,
                                               v_sequence_end)
  mapped <- map_alignment_columns_to_sequence(sequence_alignment, sequence,
                                              query_indices, n_positions)

  consensus_bases <- strsplit(as.character(consensus), "", fixed = TRUE)[[1]]
  beyond <- !is.na(mapped) & mapped >= length(consensus_bases)
  if (any(beyond)) {
    stop("V position ", which(beyond)[1L], " maps to sequence position ",
         mapped[which(beyond)[1L]] + 1L, ", past the end of a ",
         length(consensus_bases), "-base consensus")
  }

  # The votes describe the contig's base. Where TRUST4's substituted CDR3 rewrote
  # it, the record's base is a different one and the votes do not belong to it.
  mismatch <- !is.na(mapped) &
    substring(sequence, mapped + 1L, mapped + 1L) != consensus_bases[mapped + 1L]
  mapped[mismatch] <- NA_integer_

  total <- rep(NA_integer_, n_positions)
  non_consensus <- rep(NA_integer_, n_positions)
  present <- !is.na(mapped)
  if (any(present)) {
    index <- mapped[present] + 1L
    total[present] <- total_votes[index]
    non_consensus[present] <- total_votes[index] - consensus_votes[index]
  }
  list(total = total, non_consensus = non_consensus,
       n_mapped = sum(present), n_mismatch = sum(mismatch),
       first_mismatch = if (any(mismatch)) which(mismatch)[1L] else NA_integer_)
}

#' Votes behind one record's junction, located in its contig consensus.
#'
#' The junction is searched for in the consensus rather than in the record's own
#' sequence, so a record whose CDR3 is not the contig's is reported as having no
#' junction here instead of being given votes that describe other bases. A
#' junction found twice is reported as ambiguous rather than resolved by taking
#' the first, because nothing here says which copy the record means.
#'
#' Returns the 1-based start of the junction in the consensus, the two vote
#' vectors over the junction window, and the status explaining an empty result.
map_junction_votes <- function(junction, consensus, total_votes, consensus_votes) {
  empty <- list(start = NA_integer_, total = NA_character_,
                non_consensus = NA_character_)
  if (is_airr_null(junction)) {
    return(c(empty, list(status = "junction_missing")))
  }
  # gregexpr reports every non-overlapping match 1-based, or a single -1 for none.
  hits <- gregexpr(as.character(junction), as.character(consensus),
                   fixed = TRUE)[[1]]
  if (hits[1L] == -1L) {
    return(c(empty, list(status = "absent_from_consensus")))
  }
  if (length(hits) > 1L) {
    return(c(empty, list(status = "ambiguous_in_consensus")))
  }
  start <- as.integer(hits[1L])
  width <- nchar(as.character(junction))
  index <- start - 1L + seq_len(width)
  list(start = start,
       total = format_vote_vector(total_votes[index], width),
       non_consensus = format_vote_vector(total_votes[index] -
                                            consensus_votes[index], width),
       status = "mapped")
}

#' Build one library's IGH consensus-vote table.
#'
#' `igh_airr` is the QC-pass IGH subset of the reannotated table and
#' `backbone_votes` the output of `read_backbone_votes()`. Returns one row per
#' record, ordered by `sequence_id`. Records sharing a contig share its V votes,
#' because they read the same contig bases. Where TRUST4's substituted CDR3
#' rewrote a base, `map_imgt_v_votes()` drops that position instead of giving it
#' the contig's votes, and `imgt_v_first_consensus_mismatch_position` reports how
#' far five-prime the substitution reached.
build_consensus_vote_rows <- function(igh_airr, backbone_votes, n_positions) {
  required <- c("library_id", "sequence_id", "sequence", "sequence_alignment",
                "v_cigar", "v_sequence_start", "v_sequence_end", "junction")
  missing <- setdiff(required, names(igh_airr))
  if (length(missing)) {
    stop("IGH AIRR lacks: ", paste(missing, collapse = ", "))
  }

  n <- nrow(igh_airr)
  rows <- data.frame(
    library_id = as.character(igh_airr$library_id),
    sequence_id = as.character(igh_airr$sequence_id),
    trust4_backbone_id = trust4_assembly_backbone_id(igh_airr$sequence_id),
    consensus_length = rep(NA_integer_, n),
    imgt_v_positions_mapped = rep(NA_integer_, n),
    imgt_v_positions_consensus_mismatch = rep(NA_integer_, n),
    imgt_v_first_consensus_mismatch_position = rep(NA_integer_, n),
    total_votes_by_imgt_position = rep(NA_character_, n),
    non_consensus_votes_by_imgt_position = rep(NA_character_, n),
    junction_start_in_consensus = rep(NA_integer_, n),
    total_votes_by_junction_position = rep(NA_character_, n),
    non_consensus_votes_by_junction_position = rep(NA_character_, n),
    junction_vote_status = rep(NA_character_, n),
    stringsAsFactors = FALSE
  )
  if (!n) return(rows)

  for (i in seq_len(n)) {
    backbone <- rows$trust4_backbone_id[i]
    votes <- backbone_votes[[backbone]]
    if (is.null(votes)) stop("no final.out votes were read for ", backbone)
    rows$consensus_length[i] <- nchar(votes$consensus)

    # A record whose sequence agrees better with the reverse complement of its
    # contig shares no reliable position with it, so no vote can be attributed to
    # it. Empty the record and record why, rather than stopping the library:
    # IgBLAST does sometimes report the receptor on the opposite strand, and that
    # is a property of the data, not a misconfiguration a rerun could fix.
    if (!record_is_forward_of_contig(igh_airr$sequence[i], votes$consensus)) {
      empty_vector <- format_vote_vector(rep(NA_integer_, n_positions), n_positions)
      rows$imgt_v_positions_mapped[i] <- 0L
      rows$total_votes_by_imgt_position[i] <- empty_vector
      rows$non_consensus_votes_by_imgt_position[i] <- empty_vector
      rows$junction_vote_status[i] <- "opposite_contig_orientation"
      next
    }

    v <- map_imgt_v_votes(
      igh_airr$sequence_alignment[i], igh_airr$sequence[i],
      igh_airr$v_cigar[i], igh_airr$v_sequence_start[i],
      igh_airr$v_sequence_end[i], votes$consensus, votes$total,
      votes$consensus_votes, n_positions
    )
    rows$imgt_v_positions_mapped[i] <- v$n_mapped
    rows$imgt_v_positions_consensus_mismatch[i] <- v$n_mismatch
    rows$imgt_v_first_consensus_mismatch_position[i] <- v$first_mismatch
    rows$total_votes_by_imgt_position[i] <- format_vote_vector(v$total, n_positions)
    rows$non_consensus_votes_by_imgt_position[i] <-
      format_vote_vector(v$non_consensus, n_positions)

    j <- map_junction_votes(igh_airr$junction[i], votes$consensus, votes$total,
                            votes$consensus_votes)
    rows$junction_start_in_consensus[i] <- j$start
    rows$total_votes_by_junction_position[i] <- j$total
    rows$non_consensus_votes_by_junction_position[i] <- j$non_consensus
    rows$junction_vote_status[i] <- j$status
  }

  rows[order(rows$sequence_id), , drop = FALSE]
}
