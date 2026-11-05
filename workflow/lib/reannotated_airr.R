# Join receptor annotations and evidence before common record QC.
# IgBLAST compares each assembled sequence with reference germline genes;
# MakeDb converts those results to AIRR rows. These helpers retain the raw
# allele calls and attach germline and constant-region alignment evidence by
# sequence ID. Missing evidence stays missing for later eligibility checks.
# build_reannotated_airr_frame() also attaches V/J assignments from gene_calls.R
# using reference metadata made by germline_reference.py. It does not filter
# records for common QC or choose assembly representatives.
#
# Inputs come from IgBLAST 1.22.0 (Ye et al. 2013, doi:10.1093/nar/gkt382) and
# Change-O 1.3.4 MakeDb and CreateGermlines (Gupta et al. 2015,
# doi:10.1093/bioinformatics/btv359). AIRR readers need the pinned image;
# parsing and candidate-completeness helpers use base R.

REANNOTATED_COLUMNS <- c(
  "sequence_id", "library_id", "locus", "sequence", "rev_comp",
  "v_call", "d_call", "j_call",
  "v_gene_status", "j_gene_status", "v_gene", "j_gene",
  "v_gene_label", "j_gene_label",
  "v_call_candidates_complete", "d_call_candidates_complete",
  "j_call_candidates_complete", "c_call_candidates_complete",
  "v_reference_functionality", "d_reference_functionality",
  "j_reference_functionality", "c_reference_functionality",
  "v_reference_orphon", "d_reference_orphon", "j_reference_orphon", "c_reference_orphon",
  "c_call", "c_informative_positions",
  "junction", "junction_aa", "cdr3", "cdr3_aa",
  "sequence_alignment", "germline_alignment", "germline_alignment_d_mask",
  "v_cigar", "d_cigar", "j_cigar",
  "v_sequence_start", "v_sequence_end", "v_germline_start", "v_germline_end",
  "d_sequence_start", "d_sequence_end", "j_sequence_start", "j_sequence_end",
  "v_identity", "j_identity", "productive", "vj_in_frame", "stop_codon",
  "junction_length", "consensus_count"
)

#' Convert pinned MakeDb CIGARs to AIRR query/reference semantics.
#'
#' Change-O 1.3.4 decodeBTOP reverses I and D. This adapter applies only to its
#' output, before it enters the canonical pipeline table. Empty segment evidence
#' stays missing. Coordinates describe the oriented query, with 1-based closed
#' intervals; IMGT-gapped germline coordinates are not query spans.
#' Raw MakeDb tables are retained separately by the annotation wrapper.
normalise_makedb_cigar <- function(cigar, start, end, sequence_length, sequence_id, segment) {
  cigar <- chartr("ID", "DI", as.character(cigar))
  for (i in seq_along(cigar)) {
    if (is_airr_null(cigar[i])) next
    context <- paste0("MakeDb ", sequence_id[i], " ", segment, "_cigar: ")
    # S/N occur only at the ends; they record offsets outside the aligned span.
    if (!grepl("^([1-9][0-9]*S)?([1-9][0-9]*N)?([1-9][0-9]*[=XMID])+([1-9][0-9]*S)?([1-9][0-9]*N)?$",
               cigar[i])) stop(context, "invalid AIRR operations: ", cigar[i])
    tokens <- regmatches(cigar[i], gregexpr("[0-9]+[=XMIDNS]", cigar[i]))[[1]]
    lengths <- as.numeric(sub("[=XMIDNS]$", "", tokens))
    operations <- sub("^[0-9]+", "", tokens)
    query_span <- sum(lengths[operations %in% c("=", "X", "M", "I")])
    query_start <- if (operations[1] == "S") lengths[1] + 1 else 1
    first <- suppressWarnings(as.numeric(start[i]))
    last <- suppressWarnings(as.numeric(end[i]))
    if (!is.finite(first) || !is.finite(last) || first != query_start ||
        last != first + query_span - 1 || last > sequence_length[i]) {
      stop(context, "query coordinates contradict ", cigar[i],
           " (", start[i], "..", end[i], ")")
    }
  }
  cigar
}

#' Resolve a TRA/TRD V call from concordant J and optional D/C evidence.
#'
#' IgBLAST can label a delta rearrangement TRA from its shared V hit. Every V
#' candidate must support both TRA and TRD before the J call can decide between
#' them. Present D and C calls must agree with that decision; missing D or C
#' calls are allowed. Other annotations and segment conflicts stay unchanged
#' for common QC. Apply this before the locus-agreement audit as well as export.
resolve_shared_trav_dv_locus <- function(records) {
  require_columns(records, c("locus", "v_call", "j_call", "d_call", "c_call"),
                  "MakeDb-pass records")
  shared_v <- gene_call_locus_compatible(records$v_call, "TRA", "V") &
    gene_call_locus_compatible(records$v_call, "TRD", "V")
  for (locus in c("TRA", "TRD")) {
    compatible <- shared_v & gene_call_locus_compatible(records$j_call, locus, "J") &
      (is_airr_null(records$d_call) |
         gene_call_locus_compatible(records$d_call, locus, "D")) &
      (is_airr_null(records$c_call) |
         gene_call_locus_compatible(records$c_call, locus, "C"))
    records$locus[compatible] <- locus
  }
  records
}

#' Read an AIRR Rearrangement table that the step before it may not have written.
#'
#' `path` is the file, or the empty string when the step that would have produced
#' it had nothing to run on; `label` names the table in the error message so the
#' reader knows which of the stage's several inputs is missing. Returns the table,
#' or a frame with no columns when there is nothing to read.
read_optional_airr_table <- function(path, label) {
  if (!nzchar(path)) return(data.frame())
  if (!file.exists(path)) stop(label, " file does not exist: ", path)
  airr::read_rearrangement(path)
}

#' Read the C-alignment evidence out of an IgBLAST fmt7 hit table that may be absent.
#'
#' `path` is the fmt7 file, or the empty string when IgBLAST was not run for that
#' receptor family; `label` names it in the error message. Returns whatever
#' parse_igblast_c_evidence() makes of the lines, which for no file is its
#' zero-row frame.
read_optional_c_evidence <- function(path, label) {
  if (!nzchar(path)) return(parse_igblast_c_evidence(character(0)))
  if (!file.exists(path)) stop(label, " file does not exist: ", path)
  parse_igblast_c_evidence(readLines(path, warn = FALSE))
}

#' Translate the CDR3s whose nucleotides make whole codons.
#'
#' `nt` are CDR3 nucleotide sequences. A CDR3 whose length is not a multiple of
#' three cannot be read as codons, and neither can an AIRR null or an empty
#' string, so those come back NA rather than as a truncated or shifted peptide.
#' The rest go through Alakazam 1.4.3 translateDNA(trim = FALSE). Returns a
#' character vector as long as `nt`.
translate_in_frame_cdr3 <- function(nt) {
  nt <- as.character(nt)
  out <- rep(NA_character_, length(nt))
  usable <- !is_airr_null(nt) & nchar(nt) > 0L & nchar(nt) %% 3L == 0L
  if (any(usable)) {
    out[usable] <- suppressWarnings(alakazam::translateDNA(nt[usable], trim = FALSE))
  }
  out
}

#' Attach the D-masked germline alignment to the records it was reconstructed for.
#'
#' `records` are the MakeDb-pass rows and `germlines` the CreateGermlines
#' germline-pass table, joined on `sequence_id`. Germlines are reconstructed for
#' the IG loci only, so a TR record, and an IG record CreateGermlines could not
#' reconstruct, keeps an NA. Downstream IGH SHM and BASELINe read this
#' column, which is why it is carried through the shared schema rather than being
#' re-derived per analysis.
#'
#' Returns `records` with `germline_alignment_d_mask` set, and `n_attached`, the
#' number of records a germline was found for.
attach_d_masked_germline <- function(records, germlines) {
  germ_index <- match(as.character(records$sequence_id), as.character(germlines$sequence_id))
  dmask <- rep(NA_character_, nrow(records))
  has_germ <- !is.na(germ_index)
  if (any(has_germ)) {
    dmask[has_germ] <- as.character(
      column_or_default(germlines, "germline_alignment_d_mask")[germ_index[has_germ]])
  }
  records$germline_alignment_d_mask <- dmask
  list(records = records, n_attached = sum(has_germ))
}

#' Attach the informative-position count of each record's best C alignment.
#'
#' `records` are the MakeDb-pass rows and `evidence` the frame
#' parse_igblast_c_evidence() returned, joined on `sequence_id`. That function's
#' documentation below is where the fmt7 hit ordering, and why
#' the first C row is the right one to read, are explained. The count is carried
#' forward so the isotype-call qualification stage can mask poorly supported C
#' calls without discarding otherwise usable receptor records.
#'
#' A summary C call can lack a fmt7 alignment. Preserve its missing evidence
#' count so isotype_calls.R masks qualification while retaining the record.
#'
#' Returns `records` with `c_informative_positions` set, `n_attached`, the number
#' of records evidence was found for, and `n_call_without_evidence`, the number
#' of C calls left with no alignment behind them for the caller to report.
attach_constant_region_evidence <- function(records, evidence) {
  c_call <- as.character(column_or_default(records, "c_call"))
  c_index <- match(as.character(records$sequence_id), as.character(evidence$sequence_id))
  has_c_call <- !is_airr_null(c_call)
  c_call_without_evidence <- sum(has_c_call & is.na(c_index))
  c_informative_positions <- rep(NA_integer_, nrow(records))
  has_c_evidence <- !is.na(c_index)
  if (any(has_c_evidence)) {
    c_informative_positions[has_c_evidence] <-
      evidence$c_informative_positions[c_index[has_c_evidence]]
  }
  records$c_informative_positions <- c_informative_positions
  list(records = records, n_attached = sum(has_c_evidence),
       n_call_without_evidence = c_call_without_evidence)
}

#' Project the joined records onto the released all-locus schema.
#'
#' `records` are the MakeDb-pass rows with the germline and C-alignment evidence
#' already attached, `library_id` the library every row is stamped with,
#' `consensus_count` the TRUST4 read support of each row, in the same order, and
#' `identity_classes` the lookup from original allele names to groups of identical
#' complete reference sequences, read from the reference metadata.
#'
#' The raw allele lists remain in `v_call` and `j_call`. gene_calls.R supplies
#' each call's status, grouping token and display label. Several named genes
#' share a token only when their named alleles have the same complete,
#' unambiguous reference sequence. Partial reference entries or a potentially
#' truncated candidate list cannot establish that assignment. Downstream
#' grouping reads these supplied assignments rather than choosing a gene from
#' the first entry in a tied list.
#'
#' Every field is taken from `records` by name and coerced here rather than
#' wherever it came from, because the pipeline reads its own tables entirely as
#' character: the coordinates and the junction length are made integer, the
#' two identities numeric, and the four AIRR booleans are passed through untouched
#' for the AIRR writer to spell. Present CIGARs must agree with their query
#' coordinates. Unparseable identities and junction length become NA; common QC
#' decides whether the record's remaining evidence is usable.
#'
#' The MakeDb fields are required rather than defaulted. MakeDb writes both
#' receptor families through one AIRR writer, so the IG and TR tables carry the
#' same columns and a locus that has no D segment gets an empty `d_call` and
#' `d_cigar` rather than a missing column. An absent column therefore means the
#' annotation contract changed, which should stop the library instead of quietly
#' releasing a column of nulls.
#'
#' Returns a frame carrying exactly REANNOTATED_COLUMNS, ordered by locus and
#' then sequence_id as text under the C collation set_pipeline_collation()
#' fixes, so two runs produce the same row order.
build_reannotated_airr_frame <- function(records, library_id, consensus_count,
                                         identity_classes) {
  require_columns(records, c(
    "sequence_id", "locus", "sequence", "rev_comp",
    "v_call", "d_call", "j_call", "c_call", "c_informative_positions",
    "junction", "junction_aa", "cdr3",
    "sequence_alignment", "germline_alignment", "germline_alignment_d_mask",
    "v_cigar", "d_cigar", "j_cigar",
    "v_sequence_start", "v_sequence_end", "v_germline_start", "v_germline_end",
    "d_sequence_start", "d_sequence_end", "j_sequence_start", "j_sequence_end",
    "v_identity", "j_identity", "productive", "vj_in_frame", "stop_codon",
    "junction_length"
  ), "MakeDb-pass records joined with the germline and C-alignment evidence")

  for (segment in c("v", "d", "j")) {
    column <- paste0(segment, "_cigar")
    records[[column]] <- normalise_makedb_cigar(
      records[[column]], records[[paste0(segment, "_sequence_start")]],
      records[[paste0(segment, "_sequence_end")]], nchar(records$sequence),
      records$sequence_id, segment)
  }

  cdr3 <- as.character(records$cdr3)
  v_call <- as.character(records$v_call)
  j_call <- as.character(records$j_call)
  complete <- lapply(c("v", "d", "j", "c"), function(segment) {
    field <- paste0(segment, "_call_candidates_complete")
    if (field %in% names(records)) records[[field]] else rep(TRUE, nrow(records))
  })
  names(complete) <- c("v", "d", "j", "c")
  frame <- data.frame(
    sequence_id = as.character(records$sequence_id),
    library_id = rep(library_id, nrow(records)),
    locus = as.character(records$locus),
    sequence = as.character(records$sequence),
    rev_comp = records$rev_comp,
    v_call = v_call,
    d_call = as.character(records$d_call),
    j_call = j_call,
    v_gene_status = vj_tie_kind(v_call, identity_classes, complete$v),
    j_gene_status = vj_tie_kind(j_call, identity_classes, complete$j),
    v_gene = vj_canonical_token(v_call, identity_classes, complete$v),
    j_gene = vj_canonical_token(j_call, identity_classes, complete$j),
    c_call = as.character(records$c_call),
    c_informative_positions = records$c_informative_positions,
    junction = as.character(records$junction),
    junction_aa = as.character(records$junction_aa),
    cdr3 = cdr3,
    cdr3_aa = translate_in_frame_cdr3(cdr3),
    sequence_alignment = as.character(records$sequence_alignment),
    germline_alignment = as.character(records$germline_alignment),
    germline_alignment_d_mask = as.character(records$germline_alignment_d_mask),
    v_cigar = as.character(records$v_cigar),
    d_cigar = as.character(records$d_cigar),
    j_cigar = as.character(records$j_cigar),
    v_sequence_start = suppressWarnings(as.integer(records$v_sequence_start)),
    v_sequence_end = suppressWarnings(as.integer(records$v_sequence_end)),
    v_germline_start = suppressWarnings(as.integer(records$v_germline_start)),
    v_germline_end = suppressWarnings(as.integer(records$v_germline_end)),
    d_sequence_start = suppressWarnings(as.integer(records$d_sequence_start)),
    d_sequence_end = suppressWarnings(as.integer(records$d_sequence_end)),
    j_sequence_start = suppressWarnings(as.integer(records$j_sequence_start)),
    j_sequence_end = suppressWarnings(as.integer(records$j_sequence_end)),
    v_identity = suppressWarnings(as.numeric(records$v_identity)),
    j_identity = suppressWarnings(as.numeric(records$j_identity)),
    productive = records$productive,
    vj_in_frame = records$vj_in_frame,
    stop_codon = records$stop_codon,
    junction_length = suppressWarnings(as.integer(records$junction_length)),
    consensus_count = consensus_count,
    stringsAsFactors = FALSE
  )
  frame$v_gene_label <- vj_entity_label(v_call, identity_classes, complete$v)
  frame$j_gene_label <- vj_entity_label(j_call, identity_classes, complete$j)
  for (segment in c("v", "d", "j", "c")) {
    call <- records[[paste0(segment, "_call")]]
    frame[[paste0(segment, "_call_candidates_complete")]] <- complete[[segment]]
    frame[[paste0(segment, "_reference_functionality")]] <- reference_call_functionality(call, identity_classes)
    frame[[paste0(segment, "_reference_orphon")]] <- reference_call_orphon(call)
  }
  if (nrow(frame)) frame <- frame[order(frame$locus, frame$sequence_id), , drop = FALSE]
  frame[, REANNOTATED_COLUMNS, drop = FALSE]
}

#
# Constant-region alignment evidence, read out of IgBLAST fmt7 output.
#
# MakeDb records which C gene it called but not how much of the alignment
# supported the call, so the evidence has to be read from the hit table IgBLAST
# writes alongside it.

#' Extract C-hit evidence from IgBLAST fmt7 lines.
#'
#' The fmt7 hit table orders alignments by descending bit score. The first C
#' row for each query supplies the best-ranked C alignment. This function
#' counts its observed A/C/G/T query positions; it does not choose a class or
#' subclass from tied C genes. isotype_calls.R qualifies those calls later.
#'
#' A C hit row has 17 tab-separated fields, because run_igblast_all_loci.sh asks
#' IgBLAST for `7 std qseq sseq btop` and IgBLAST prepends the chain-type letter
#' and inserts its own `gaps` column:
#'   1 chain type, 2 query id, 3 subject id, 4 % identity, 5 alignment length,
#'   6 mismatches, 7 gap opens, 8 gaps, 9-10 query start/end,
#'   11-12 subject start/end, 13 evalue, 14 bit score, 15 query seq,
#'   16 subject seq, 17 BTOP.
#' Field 15 is therefore the aligned QUERY; field 16 would be the germline and
#' would report a short observation as fully callable. Only informative A/C/G/T
#' positions are counted, so gaps and ambiguous bases do not make a short
#' C-region observation look callable.
parse_igblast_c_evidence <- function(lines) {
  c_lines <- lines[startsWith(lines, "C\t")]
  if (!length(c_lines)) {
    return(data.frame(
      sequence_id = character(0),
      c_informative_positions = integer(0),
      stringsAsFactors = FALSE
    ))
  }
  fields <- strsplit(c_lines, "\t", fixed = TRUE)
  # Exact, not a lower bound: a different field count means the output format
  # changed and field 15 is no longer the aligned query.
  if (any(lengths(fields) != 17L)) {
    stop("IgBLAST fmt7 C hit row(s) do not have the expected 17 fields; ",
         "the IgBLAST output format has changed")
  }
  # IgBLAST marks minus-strand hits with this prefix, while MakeDb preserves
  # the original TRUST4 query ID. Remove it before selecting and joining hits.
  sequence_id <- sub("^reversed\\|", "", vapply(fields, `[[`, character(1), 2L))
  if (any(!nzchar(sequence_id))) stop("IgBLAST fmt7 contains a C hit without query ID")
  top <- !duplicated(sequence_id)
  query_alignment <- vapply(fields[top], `[[`, character(1), 15L)
  data.frame(
    sequence_id = sequence_id[top],
    c_informative_positions = as.integer(nchar(gsub(
      "[^ACGTacgt]", "", query_alignment
    ))),
    stringsAsFactors = FALSE
  )
}

#' Read original allele metadata and its complete-reference sequence classes.
#'
#' The lookup names alleles, before any call is reduced to gene names. Missing
#' identity tokens record partial or ambiguous reference sequence evidence.
read_germline_identity_classes <- function(path) {
  metadata <- read_pipeline_table(path, required_columns = c(
    "allele", "gene", "region", "functionality", "completeness",
    "identity_id", "identity_label", "orphon"), context = "Reference allele metadata")
  if (anyDuplicated(metadata$allele)) stop("Reference metadata names an allele more than once")
  classes <- stats::setNames(as.character(metadata$identity_id), metadata$allele)
  usable <- !is_airr_null(metadata$identity_id)
  attr(classes, "labels") <- stats::setNames(metadata$identity_label[usable], metadata$identity_id[usable])
  attr(classes, "metadata") <- metadata
  classes
}

#' Report functionality for each named reference allele, preserving mixed calls.
reference_call_functionality <- function(call, identity_classes) {
  metadata <- attr(identity_classes, "metadata")
  vapply(call, function(cell) {
    alleles <- sort(unique(.split_calls(cell)))
    if (!length(alleles)) return(NA_character_)
    status <- metadata$functionality[match(alleles, metadata$allele)]
    if (!length(status)) status <- rep(NA_character_, length(alleles))
    status[is_airr_null(status)] <- "unknown"
    paste(paste(alleles, status, sep = ":"), collapse = ";")
  }, character(1), USE.NAMES = FALSE)
}

#' Identify calls naming an IMGT orphon location without inferring receptor function.
reference_call_orphon <- function(call) {
  vapply(call, function(cell) {
    alleles <- .split_calls(cell)
    if (!length(alleles)) return(NA)
    any(grepl("/OR", alleles, fixed = TRUE))
  }, logical(1), USE.NAMES = FALSE)
}

#' Count equally best-scoring hits before deciding whether a candidate list is complete.
#'
#' Each fmt7 V/D/J/C hit supplies an allele and bit score for one query.
#' Return the number of distinct alleles at the highest score per query and
#' segment. Reaching the configured segment cap can hide additional tied
#' candidates. Lower-scoring hits at the cap do not make the best-score list
#' incomplete. Empty input returns a header-only frame.
parse_igblast_candidate_counts <- function(lines) {
  rows <- lines[grepl("^[VDJC]\\t", lines)]
  if (!length(rows)) return(data.frame(sequence_id = character(), segment = character(), n_candidates = integer()))
  fields <- strsplit(rows, "\t", fixed = TRUE)
  if (any(lengths(fields) != 17L)) stop("IgBLAST hit rows must contain 17 fields")
  hits <- data.frame(sequence_id = sub("^reversed\\|", "", vapply(fields, `[[`, "", 2L)),
    segment = vapply(fields, `[[`, "", 1L), allele = vapply(fields, `[[`, "", 3L),
    score = as.numeric(vapply(fields, `[[`, "", 14L)), stringsAsFactors = FALSE)
  if (any(!is.finite(hits$score))) stop("IgBLAST hit scores must be finite")
  grouped <- split(seq_len(nrow(hits)), paste(hits$sequence_id, hits$segment, sep = "\x1f"))
  do.call(rbind, lapply(grouped, function(indices) {
    group <- hits[indices, , drop = FALSE]
    data.frame(sequence_id = group$sequence_id[1], segment = group$segment[1],
      n_candidates = length(unique(group$allele[group$score == max(group$score)])))
  }))
}

#' Mark whether each segment's best-score candidates fit below its alignment cap.
#'
#' `records` holds MakeDb calls, `counts` comes from parse_igblast_candidate_counts()
#' and `caps` holds the four params.igblast.num_alignments_* settings. Use the
#' larger of the hit count and the number of named call alleles, so absent hit
#' evidence cannot make a capped raw call appear complete. Preserve the calls
#' and add one logical completeness column per segment. V/J classification
#' and constant-call qualification read these flags; this function drops no rows.
attach_candidate_completeness <- function(records, counts, caps) {
  for (segment in c("v", "d", "j", "c")) {
    local <- counts[counts$segment == toupper(segment), , drop = FALSE]
    count <- local$n_candidates[match(records$sequence_id, local$sequence_id)]
    reported <- vapply(records[[paste0(segment, "_call")]], function(call) length(unique(.split_calls(call))), integer(1))
    count <- pmax(count, reported, na.rm = TRUE)
    records[[paste0(segment, "_call_candidates_complete")]] <- count < caps[[segment]]
  }
  records
}
