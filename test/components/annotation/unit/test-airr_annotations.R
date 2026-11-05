# Check annotation projection, germline joins and AIRR coordinates.
# Run: python test/run_tests.py --component annotation --level unit

load_test_modules(c("table_io", "airr_io", "gene_calls", "reannotated_airr"))

# MakeDb-pass rows as the pipeline sees them: coordinates and identities arrive
# as text, booleans as logicals, and a locus's inapplicable fields as NA.
MAKEDB <- data.frame(
  sequence_id = c("s1", "s2", "s3"),
  locus = c("IGH", "IGH", "TRB"),
  sequence = c("AAACCC", "GGGTTT", "ACGTAC"),
  rev_comp = c(FALSE, FALSE, TRUE),
  v_call = c("IGHV1-8*01", "IGHV3-23*01", "TRBV5-1*01"),
  d_call = c("IGHD3-3*01", NA, NA),
  j_call = c("IGHJ5*02", "IGHJ4*02", "TRBJ2-1*01"),
  c_call = c("IGHM", "IGHG1", NA),
  junction = c("TGTGCCTGG", "TGTGCC", "TGTGCCAGG"),
  junction_aa = c("CAW", "CA", "CAR"),
  cdr3 = c("TGTGCCTGG", "TGTGC", "TGTGCCAGG"),
  sequence_alignment = c("AAACCC", "GGGTTT", "ACGTAC"),
  germline_alignment = c("AAACCC", "GGGTTT", "ACGTAC"),
  v_cigar = "6M", d_cigar = c("2M", NA, NA), j_cigar = "6M",
  d_sequence_start = c("1", NA, NA), d_sequence_end = c("2", NA, NA),
  j_sequence_start = "1", j_sequence_end = "6",
  v_sequence_start = "1", v_sequence_end = "6",
  v_germline_start = "1", v_germline_end = "6",
  v_identity = c("0.98", "bogus", "1.0"), j_identity = "1.0",
  productive = c(TRUE, FALSE, TRUE),
  vj_in_frame = c(TRUE, TRUE, FALSE),
  stop_codon = c(FALSE, FALSE, TRUE),
  junction_length = c("9", "6", "9"),
  stringsAsFactors = FALSE
)

# --- The released schema ------------------------------------------------------

# A constructed lookup that treats two IMGT gene names as one germline entity,
# as the pinned human reference does for six V/J gene pairs. Passed to the
# projection the way the reannotation stage passes the lookup it read from the
# reference.
IDENTITY_CLASSES <- c("IGHV1-8*01" = "VSEQ_TEST", "IGHV1-69*01" = "VSEQ_TEST")

test_that("the released column set and its order are fixed", {
  expect_identical(REANNOTATED_COLUMNS, c(
    "sequence_id", "library_id", "locus", "sequence", "rev_comp",
    "v_call", "d_call", "j_call",
    "v_gene_status", "j_gene_status", "v_gene", "j_gene",
    "v_gene_label", "j_gene_label",
    "v_call_candidates_complete", "d_call_candidates_complete", "j_call_candidates_complete", "c_call_candidates_complete",
    "v_reference_functionality", "d_reference_functionality", "j_reference_functionality", "c_reference_functionality",
    "v_reference_orphon", "d_reference_orphon", "j_reference_orphon", "c_reference_orphon",
    "c_call", "c_informative_positions",
    "junction", "junction_aa", "cdr3", "cdr3_aa",
    "sequence_alignment", "germline_alignment", "germline_alignment_d_mask",
    "v_cigar", "d_cigar", "j_cigar",
    "v_sequence_start", "v_sequence_end", "v_germline_start", "v_germline_end",
    "d_sequence_start", "d_sequence_end", "j_sequence_start", "j_sequence_end",
    "v_identity", "j_identity", "productive", "vj_in_frame", "stop_codon",
    "junction_length", "consensus_count"
  ))
})

# --- translate_in_frame_cdr3() ------------------------------------------------

describe("translating CDR3 nucleotides", {
  it("translates only the lengths that make whole codons", {
    expect_identical(
      translate_in_frame_cdr3(c("TGTGCCTGG", "TGTGCC", "TGTGC", "", NA,
                                "TGTNNNTGG")),
      c("CAW", "CA", NA, NA, NA, "CXW")
    )
  })
  it("returns nothing for no sequences", {
    expect_identical(translate_in_frame_cdr3(character(0)), character(0))
  })
})

# --- read_optional_airr_table() and read_optional_c_evidence() ----------------

describe("an input the step before it never wrote", {
  it("reads an absent AIRR table as no columns and no rows", {
    empty <- read_optional_airr_table("", "IG MakeDb-pass")
    expect_identical(nrow(empty), 0L)
    expect_identical(ncol(empty), 0L)
  })
  it("reads absent C evidence as the parser's own empty frame", {
    expect_identical(read_optional_c_evidence("", "IG IgBLAST fmt7"),
                     parse_igblast_c_evidence(character(0)))
  })
  it("names which input is missing when a path was given but does not exist", {
    expect_error(read_optional_airr_table("/nope/nope.tsv", "IG MakeDb-pass"),
                 "IG MakeDb-pass file does not exist", fixed = TRUE)
    expect_error(read_optional_c_evidence("/nope/x.fmt7", "IG IgBLAST fmt7"),
                 "IG IgBLAST fmt7 file does not exist", fixed = TRUE)
  })
})

# --- attach_d_masked_germline() -----------------------------------------------

describe("joining the D-masked germline", {
  germlines <- data.frame(
    sequence_id = c("s2", "sX"),
    germline_alignment_d_mask = c("GGGNNN", "ZZZ"),
    stringsAsFactors = FALSE
  )
  attached <- attach_d_masked_germline(MAKEDB, germlines)

  it("carries the germline onto the record it belongs to", {
    expect_identical(attached$records$germline_alignment_d_mask,
                     c(NA, "GGGNNN", NA))
  })
  it("counts what it attached, for the caller to log", {
    expect_identical(attached$n_attached, 1L)
  })
  it("ignores a germline for a record that is not here", {
    expect_identical(nrow(attached$records), 3L)
  })
})

test_that("no germlines at all leaves every record with an NA", {
  attached <- attach_d_masked_germline(MAKEDB, data.frame())
  expect_identical(attached$records$germline_alignment_d_mask,
                   rep(NA_character_, 3L))
  expect_identical(attached$n_attached, 0L)
})

# --- attach_constant_region_evidence() ----------------------------------------

describe("joining the first best-ranked C alignment", {
  evidence <- data.frame(
    sequence_id = c("s1", "s3"),
    c_informative_positions = c(120L, 45L),
    stringsAsFactors = FALSE
  )
  attached <- attach_constant_region_evidence(MAKEDB, evidence)

  it("carries the informative-position count onto its record", {
    expect_identical(attached$records$c_informative_positions,
                     c(120L, NA, 45L))
  })
  it("counts what it attached, for the caller to log", {
    expect_identical(attached$n_attached, 2L)
  })
  # s2 has a c_call but no C hit row; s3 has evidence but no c_call, so only s2
  # is a call the isotype-call qualification stage will have to mask.
  it("counts the C calls left with no alignment behind them", {
    expect_identical(attached$n_call_without_evidence, 1L)
  })
})

test_that("no C evidence at all leaves every count NA", {
  attached <- attach_constant_region_evidence(
    MAKEDB, parse_igblast_c_evidence(character(0)))
  expect_identical(attached$records$c_informative_positions,
                   rep(NA_integer_, 3L))
  expect_identical(attached$n_attached, 0L)
  expect_identical(attached$n_call_without_evidence, 2L)
})

# --- build_reannotated_airr_frame() -------------------------------------------

describe("projecting three joined records onto the released schema", {
  records <- attach_constant_region_evidence(
    attach_d_masked_germline(
      MAKEDB,
      data.frame(sequence_id = "s2", germline_alignment_d_mask = "GGGNNN",
                 stringsAsFactors = FALSE)
    )$records,
    data.frame(sequence_id = c("s1", "s3"),
               c_informative_positions = c(120L, 45L), stringsAsFactors = FALSE)
  )$records
  frame <- build_reannotated_airr_frame(records, "GUSS_TEST", c(11, 7, 3),
                                        IDENTITY_CLASSES)

  it("emits exactly the released columns, in order", {
    expect_identical(names(frame), REANNOTATED_COLUMNS)
  })
  it("orders rows by locus and then sequence_id", {
    expect_identical(frame$sequence_id, c("s1", "s2", "s3"))
    expect_identical(frame$locus, c("IGH", "IGH", "TRB"))
  })
  it("stamps the library on every row", {
    expect_identical(frame$library_id, rep("GUSS_TEST", 3L))
  })
  it("coerces the coordinates and the junction length to integer", {
    expect_identical(frame$v_sequence_start, rep(1L, 3L))
    expect_identical(frame$v_sequence_end, rep(6L, 3L))
    expect_identical(frame$v_germline_start, rep(1L, 3L))
    expect_identical(frame$v_germline_end, rep(6L, 3L))
    expect_identical(frame$junction_length, c(9L, 6L, 9L))
  })
  it("coerces the identities to numeric, leaving an unparseable cell NA", {
    expect_identical(frame$v_identity, c(0.98, NA, 1.00))
    expect_identical(frame$j_identity, rep(1.00, 3L))
  })
  it("leaves the AIRR booleans as logicals for the writer to spell", {
    expect_identical(frame$rev_comp, c(FALSE, FALSE, TRUE))
    expect_identical(frame$productive, c(TRUE, FALSE, TRUE))
    expect_identical(frame$vj_in_frame, c(TRUE, TRUE, FALSE))
    expect_identical(frame$stop_codon, c(FALSE, FALSE, TRUE))
  })
  it("translates the CDR3 it can and leaves the out-of-frame one NA", {
    expect_identical(frame$cdr3_aa, c("CAW", NA, "CAR"))
  })
  it("carries the joined germline and C evidence through", {
    expect_identical(frame$germline_alignment_d_mask, c(NA, "GGGNNN", NA))
    expect_identical(frame$c_informative_positions, c(120L, NA, 45L))
  })
  it("carries TRUST4's read support through", {
    expect_identical(frame$consensus_count, c(11, 7, 3))
  })
  # Each call names one gene, so each resolves; the V calls are mapped through
  # their identity class, which for s1 and s2 is the two-gene one above.
  it("says how far each call reduces to one germline entity", {
    expect_identical(frame$v_gene_status, rep("resolved", 3L))
    expect_identical(frame$j_gene_status, rep("resolved", 3L))
  })
  it("names the germline entity each call resolves to", {
    expect_identical(frame$v_gene,
                     c("IGHV1-8", "IGHV3-23", "TRBV5-1"))
    expect_identical(frame$j_gene, c("IGHJ5", "IGHJ4", "TRBJ2-1"))
  })
})

# A tie the reference cannot separate is still one entity, and a tie it can is
# none: the projection is where both are decided, once, for every stage after it.
test_that("a multi-gene call is classified against the reference", {
  ambiguous <- MAKEDB
  ambiguous$v_call <- c("IGHV1-8*01,IGHV1-69*01", "IGHV3-23*01,IGHV3-30*01",
                        "TRBV5-1*01")
  frame <- build_reannotated_airr_frame(
    attach_constant_region_evidence(
      attach_d_masked_germline(ambiguous, data.frame())$records,
      parse_igblast_c_evidence(character(0)))$records,
    "GUSS_TEST", c(11, 7, 3), IDENTITY_CLASSES)

  expect_identical(frame$v_gene_status,
                   c("assigned_group", "ambiguous", "resolved"))
  expect_identical(frame$v_gene,
                   c("VSEQ_TEST", "", "TRBV5-1"))
})

test_that("one record still yields the whole released schema", {
  single <- attach_constant_region_evidence(
    attach_d_masked_germline(MAKEDB[1, , drop = FALSE], data.frame())$records,
    parse_igblast_c_evidence(character(0))
  )$records
  frame <- build_reannotated_airr_frame(single, "GUSS_TEST", 11, IDENTITY_CLASSES)
  expect_identical(names(frame), REANNOTATED_COLUMNS)
  expect_identical(nrow(frame), 1L)
  expect_identical(frame$sequence_id, "s1")
})

# A library IgBLAST annotated nothing in still owes its consumers the header, and
# the empty frame the script writes must carry exactly the same one.
test_that("no records yield the released header and nothing else", {
  path <- tempfile(fileext = ".tsv")
  write_airr_table(make_empty_frame(REANNOTATED_COLUMNS), path)
  expect_identical(readLines(path),
                   paste(REANNOTATED_COLUMNS, collapse = "\t"))
})

# ------------------------------------------------------- raw IgBLAST C evidence ---

# Two queries' C hits in IgBLAST fmt7 order, best hit first. seq1 has a
# second, lower-scoring C hit that must be ignored.
igblast_c_hit_lines <- function() {
  c(
    "# IgBLAST tabular output",
    paste(c("V", "seq1", "IGHV1-2*01", "100", "4", "0", "0", "0",
            "1", "4", "1", "4", "1e-5", "50", "ACGT", "ACGT", "4"),
          collapse = "\t"),
    paste(c("C", "seq1", "IGHM*01", "100", "9", "0", "0", "0",
            "1", "9", "1", "9", "1e-9", "90", "ACGTNN--A", "ACGTAAAAT", "4NN2-A"),
          collapse = "\t"),
    paste(c("C", "seq1", "IGHD*01", "100", "8", "0", "0", "0",
            "1", "8", "1", "8", "1e-8", "80", "ACGTACGT", "ACGTACGT", "8"),
          collapse = "\t"),
    paste(c("C", "seq2", "TRBC1*01", "100", "9", "0", "0", "0",
            "1", "9", "1", "9", "1e-7", "70", "aCgT.NN--", "ACGTAAAAT", "4NN2-A"),
          collapse = "\t")
  )
}

describe("parse_igblast_c_evidence", {
  evidence <- parse_igblast_c_evidence(igblast_c_hit_lines())

  it("keeps the first, best-scoring C row per sequence and nothing else",
     expect_identical(evidence$sequence_id, c("seq1", "seq2")))
  # Field 15 is the aligned QUERY, so gaps and ambiguous bases keep a short
  # observation from looking fully callable.
  it("counts only the A/C/G/T query positions of the best hit",
     expect_identical(evidence$c_informative_positions[[1]], 5L))
  it("counts lowercase query bases too",
     expect_identical(evidence$c_informative_positions[[2]], 4L))
})

test_that("minus-strand C hits join by original ID and retain their best hit", {
  lines <- igblast_c_hit_lines()
  reversed <- gsub("\tseq1\t", "\treversed|seq1\t", lines, fixed = TRUE)
  expect_identical(parse_igblast_c_evidence(reversed), parse_igblast_c_evidence(lines))
  # Only the leading tool prefix is removed.
  embedded <- gsub("\tseq1\t", "\tseq1_reversed|part\t", lines, fixed = TRUE)
  expect_identical(parse_igblast_c_evidence(embedded)$sequence_id,
                   c("seq1_reversed|part", "seq2"))
})


test_that("MakeDb CIGARs become AIRR CIGARs for every segment", {
  for (segment in c("v", "d", "j")) {
    raw <- c("60=3D233=", "60=3I230=", "10S4N30=3D20=3I40=", "293=", "", NA)
    expected <- c("60=3I233=", "60=3D230=", "10S4N30=3I20=3D40=", "293=", "", NA)
    result <- normalise_makedb_cigar(raw, c(1, 1, 11, 1, NA, NA),
      c(296, 290, 103, 293, NA, NA), rep(400, 6), paste0("seq", 1:6), segment)
    expect_identical(result, expected)
  }
})

test_that("a producer convention change or corrupt alignment fails with context", {
  for (cigar in c("60=3I233=", "3S60=3D233=", "60=3Q233=", "60=3D233=garbage")) {
    expect_error(normalise_makedb_cigar(cigar, 1, 296, 400, "seq1", "v"),
                 "MakeDb seq1 v_cigar", fixed = TRUE)
  }
  expect_error(normalise_makedb_cigar("20=", NA, 20, 40, "seq2", "j"),
               "query coordinates contradict", fixed = TRUE)
})
