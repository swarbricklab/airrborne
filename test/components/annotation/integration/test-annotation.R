# Run annotation tools on all seven loci and constructed alignment boundaries.
# Run: python test/run_tests.py --component annotation --level integration

load_test_modules(c("table_io", "gene_calls", "reannotated_airr", "record_qc"))

# Public germline-derived inputs exercise real tool calls and annotation losses.
# Expected loci, junctions and V/J genes follow the declared sequence joins.

annotation_fixture <- function() {
  directory <- file.path(MODULE_ROOT, "test", "data", "annotation")
  lines <- readLines(file.path(directory, "segments.fasta"))
  reference <- setNames(lines[seq(2, length(lines), 2)],
                         substring(lines[seq(1, length(lines), 2)], 2))
  joins <- read.delim(file.path(directory, "joins.tsv"), stringsAsFactors = FALSE)
  sequences <- vapply(seq_len(nrow(joins)), function(i) {
    row <- joins[i, ]
    v <- gsub(".", "", substring(reference[[row$v_call]], 1, row$v_end), fixed = TRUE)
    paste0(v, row$insert, substring(reference[[row$j_call]], row$j_start + 1))
  }, character(1))
  list(records = data.frame(sequence_id = paste0(joins$locus, "_0"),
    sequence = sequences, locus = joins$locus, consensus_count = c(2, 3, 5, 7, 11, 13, 17)),
    joins = joins, alleles = reference, reference = file.path(directory, "segments.fasta"))
}

# Build the pipeline's reference and IgBLAST auxiliary data from the image, as
# the obtain_imgt_reference rule does. Returns the two written paths and the
# image's own auxiliary data file, which the written one extends.
pipeline_reference_files <- function(directory) {
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  script <- file.path(MODULE_ROOT, "workflow", "scripts", "obtain_imgt_reference.py")
  paths <- list(reference = file.path(directory, "human_IMGT+C.fa"),
                auxiliary = file.path(directory, "human_gl.aux"),
                metadata = file.path(directory, "reference.metadata.tsv"),
                shipped = "/usr/local/share/igblast/optional_file/human_gl.aux")
  output <- suppressWarnings(system2("python3", shQuote(c(script,
    "--germline-dir", "/usr/local/share/germlines/imgt/human/vdj",
    "--constant-dir", "/usr/local/share/germlines/imgt/human/constant",
    "--auxiliary-data", paths$shipped, "--output", paths$reference,
    "--auxiliary-output", paths$auxiliary, "--metadata-output", paths$metadata)), stdout = TRUE, stderr = TRUE))
  status <- attr(output, "status")
  expect_script_ok(list(script = script, output = output,
    status = if (is.null(status)) 0L else as.integer(status)))
  paths
}

# Entries of an IgBLAST auxiliary data file, named by allele: the frame, chain
# type, CDR3 end and extra bases as the file lists them.
auxiliary_entries <- function(path) {
  rows <- trimws(readLines(path))
  rows <- rows[nzchar(rows) & !startsWith(rows, "#")]
  fields <- lapply(strsplit(rows, "[[:space:]]+"), function(x) x[-1])
  setNames(fields, vapply(strsplit(rows, "[[:space:]]+"), `[`, "", 1))
}

# Write one native contig with ten agreeing votes per called base and none at N.
# The caller supplies TRUST4 assembly IDs so the vote mapper can find each row.
write_native_votes <- function(native, path) {
  votes <- unlist(lapply(seq_len(nrow(native)), function(i) {
    bases <- strsplit(native$sequence[i], "", fixed = TRUE)[[1]]
    c(paste0(">", sub("_0$", "", native$sequence_id[i])), native$sequence[i],
      vapply(c("A", "C", "G", "T"), function(base) {
        paste(ifelse(bases == base, 10, 0), collapse = " ")
      }, character(1)))
  }), use.names = FALSE)
  writeLines(votes, path)
}

# `auxiliary` defaults to the file the pipeline writes, built into `directory`.
run_annotation <- function(records, directory, auxiliary = NULL) {
  dir.create(directory)
  reference_files <- pipeline_reference_files(file.path(directory, "reference"))
  if (is.null(auxiliary)) auxiliary <- reference_files$auxiliary
  scratch <- file.path(directory, "scratch")
  dir.create(scratch)
  native <- file.path(directory, "native.tsv")
  cdr3 <- file.path(directory, "cdr3.out")
  write_pipeline_table(records, native)
  scores <- data.frame(consensus = sub("_0$", "", records$sequence_id),
    index = rep(0, nrow(records)), v = rep(".", nrow(records)),
    d = rep(".", nrow(records)), j = rep(".", nrow(records)), c = rep(".", nrow(records)),
    cdr1 = rep(".", nrow(records)), cdr2 = rep(".", nrow(records)),
    cdr3 = rep("TGTGCC", nrow(records)), score = rep(1, nrow(records)),
    reads = records$consensus_count, similarity = rep(100, nrow(records)),
    complete = rep(1, nrow(records)))
  write.table(scores, cdr3, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  paths <- file.path(directory, c("annotated.tsv", "audit.tsv", "failed.tsv", "makedb.log.gz"))
  script <- file.path(MODULE_ROOT, "workflow", "scripts", "run_igblast_all_loci.sh")
  output <- suppressWarnings(system2("bash", shQuote(c(script, native, paths[1], scratch,
    "0040", "1", "/usr/local/share/igblast", "/usr/local/share/germlines/imgt/human/vdj",
    cdr3, paths[2:4], "50", "20", "3", "125",
    file.path(directory, "makedb_pass.tsv.gz"), auxiliary, reference_files$metadata, file.path(directory, "igblast_hits.fmt7.gz"))), stdout = TRUE, stderr = TRUE))
  status <- attr(output, "status")
  expect_script_ok(list(script = script, output = output,
    status = if (is.null(status)) 0L else as.integer(status)))
  list(records = read_pipeline_table(paths[1]), audit = read_pipeline_table(paths[2]),
       failed = read_pipeline_table(paths[3]), scratch = scratch, output = output)
}

test_that("all loci are annotated and losses retain their native denominators", {
  fixture <- annotation_fixture()
  records <- fixture$records
  moved <- records[records$locus == "TRD", ]
  moved$sequence_id <- "reassigned_0"
  moved$locus <- "TRA"
  moved$sequence <- sub("GCCGCCGCCGCCGCC", "GCCGCCGCCGCCGCCGCC", moved$sequence, fixed = TRUE)
  moved$consensus_count <- 29
  records <- rbind(records, moved,
    data.frame(sequence_id = c("unannotated_0", "unrouted_0"),
      sequence = c(strrep("A", 360), strrep("C", 360)), locus = c("IGH", "OTHER"),
      consensus_count = c(19, 23)))
  directory <- tempfile("annotation-all-")
  on.exit(unlink(directory, recursive = TRUE), add = TRUE)
  result <- run_annotation(records, directory)
  expected_ids <- c(fixture$records$sequence_id, "reassigned_0")
  expect_setequal(result$records$sequence_id, expected_ids)
  matched <- result$records[match(expected_ids, result$records$sequence_id), ]
  expect_identical(matched$locus, c(fixture$joins$locus, "TRD"))
  expect_equal(as.numeric(matched$consensus_count), c(2, 3, 5, 7, 11, 13, 17, 29))
  for (i in seq_len(nrow(fixture$joins))) {
    genes <- strsplit(gsub("\\*[^,]*", "", matched$v_call[i]), ",", fixed = TRUE)[[1]]
    expect_true(sub("\\*.*", "", fixture$joins$v_call[i]) %in% genes)
    expect_true(grepl(sub("\\*.*", "", fixture$joins$j_call[i]), matched$j_call[i], fixed = TRUE))
    expect_true(matched$productive[i] %in% c("T", "TRUE", "true"))
    expect_identical(matched$junction_aa[i], if (fixture$joins$locus[i] == "IGH") "CAAAAAW" else "CAAAAAF")
  }
  expect_identical(result$failed$sequence_id, "unannotated_0")
  expect_identical(result$failed$receptor_family, "ig")
  for (family in c("ig", "tr")) {
    fasta <- readLines(file.path(result$scratch, paste0(family, ".fasta")))
    ids <- substring(fasta[startsWith(fasta, ">")], 2)
    expected <- records$sequence_id[if (family == "ig") records$locus %in% c("IGH", "IGK", "IGL") else records$locus %in% c("TRA", "TRB", "TRD", "TRG")]
    expect_setequal(ids, expected)
  }
  counts <- function(stage, locus = "all") {
    row <- result$audit[result$audit$stage == stage & result$audit$locus == locus, ]
    expect_equal(nrow(row), 1)
    as.numeric(row[c("n_units_in", "n_units_out", "consensus_count_in", "consensus_count_out")])
  }
  expect_equal(counts("igblast_locus_routing"), c(10, 9, 129, 106))
  expect_equal(counts("igblast_makedb"), c(9, 8, 106, 87))
  expect_equal(counts("igblast_makedb", "IGH"), c(2, 1, 21, 2))
  expect_equal(counts("igblast_locus_agreement", "TRD"), c(2, 1, 42, 13))
  agreement <- result$audit[result$audit$stage == "igblast_locus_agreement", ]
  expect_true(all(agreement$kind == "diagnostic"))
  expect_true(all(agreement$locus_basis == "annotated"))
  # A native/annotated disagreement is diagnostic. Common QC must use the
  # returned genes and locus, including the reassigned TRD record.
  qc <- file.path(directory, "qc_pass.tsv")
  qc_audit <- file.path(directory, "qc_audit.tsv")
  expect_script_ok(run_pipeline_script("filter_airr.R", c(
    paste0("--input_file=", file.path(directory, "annotated.tsv")),
    paste0("--output_file=", qc), "--library_id=0040", "--max_n_fraction=0.1",
    paste0("--qc_audit_file=", qc_audit))))
  expect_setequal(read_pipeline_table(qc)$sequence_id, expected_ids)
  decisions <- read_pipeline_table(qc_audit)
  expect_true(all(decisions$n_units_in == decisions$n_units_out))
})

test_that("J alleles without an NCBI auxiliary entry receive a copied entry and keep their junctions", {
  directory <- tempfile("annotation-auxiliary-")
  on.exit(unlink(directory, recursive = TRUE), add = TRUE)
  files <- pipeline_reference_files(file.path(directory, "reference"))
  lines <- readLines(files$reference)
  reference <- setNames(lines[seq(2, length(lines), 2)], substring(lines[seq(1, length(lines), 2)], 2))
  shipped <- auxiliary_entries(files$shipped)
  written <- auxiliary_entries(files$auxiliary)
  # Each allele the image's file lacks, the allele whose shipped entry the
  # pipeline copies for it, and a V partner from the constructed panel.
  cases <- data.frame(
    allele = c("IGHJ5*04", "IGKJ4*03", "IGLJ2A*01", "TRAJ41*02", "TRAJ57*02", "TRDJ4*02", "TRGJP2*02"),
    source = c("IGHJ5*02", "IGKJ4*01", "IGLJ2*01", "TRAJ41*01", "TRAJ57*01", "TRDJ4*01", "TRGJP2*01"),
    v_call = c("IGHV2-5*01", "IGKV1-12*01", "IGLV1-40*01", "TRAV1-1*01", "TRAV1-1*01", "TRDV1*01", "TRGV1*01"),
    locus = c("IGH", "IGK", "IGL", "TRA", "TRA", "TRD", "TRG"), stringsAsFactors = FALSE)
  for (i in seq_len(nrow(cases))) {
    expect_false(cases$allele[i] %in% names(shipped))
    expect_identical(written[[cases$allele[i]]], shipped[[cases$source[i]]])
  }
  # A pseudogene of a gene with no other allele keeps no entry, as in NCBI's file.
  expect_false("TRAJ55*01" %in% names(written))
  expect_false("TRAJ55*01" %in% names(shipped))

  # V through the conserved Cys104 (gapped position 312), five alanine codons
  # padded so the J's coding frame from the shipped entry is in phase, then the
  # whole J. The expected junction runs from Cys104 to the conserved codon that
  # follows the shipped CDR3 end, so it depends on the copied values.
  construct <- function(v_call, j_call, entry) {
    v <- gsub(".", "", substring(reference[[v_call]], 1, 312), fixed = TRUE)
    j <- gsub(".", "", reference[[j_call]], fixed = TRUE)
    frame <- as.integer(entry[1])
    stop <- as.integer(entry[3])
    insert <- paste0("GCCGCCGCCGCCGCC", strrep("G", (3 - frame) %% 3))
    list(sequence = paste0(v, insert, j),
         junction = paste0(substring(v, nchar(v) - 2), insert, substring(j, 1, stop + 4)))
  }
  id <- function(allele) paste0(gsub("[*/]", "_", allele), "_0")
  built <- list()
  for (i in seq_len(nrow(cases))) {
    entry <- shipped[[cases$source[i]]]
    built[[id(cases$allele[i])]] <- c(construct(cases$v_call[i], cases$allele[i], entry), locus = cases$locus[i])
    built[[id(cases$source[i])]] <- c(construct(cases$v_call[i], cases$source[i], entry), locus = cases$locus[i])
  }
  pseudogene <- construct("TRAV1-1*01", "TRAJ55*01", c("0", "JA", "0", "0"))
  built[[id("TRAJ55*01")]] <- c(pseudogene, locus = "TRA")
  records <- data.frame(sequence_id = names(built),
    sequence = vapply(built, `[[`, "", "sequence"), locus = vapply(built, `[[`, "", "locus"),
    consensus_count = 10, stringsAsFactors = FALSE)

  result <- run_annotation(records, file.path(directory, "supplemented"), files$auxiliary)
  kept <- names(built) != id("TRAJ55*01")
  expect_identical(result$failed$sequence_id, id("TRAJ55*01"))
  expect_setequal(result$records$sequence_id, names(built)[kept])
  annotated <- result$records
  for (name in names(built)[kept]) {
    row <- annotated[which(annotated$sequence_id == name), ]
    expect_identical(nrow(row), 1L, info = name)
    expect_identical(row$junction, built[[name]]$junction, info = name)
    expect_identical(row$productive, "T", info = name)
    expect_identical(nchar(row$junction) %% 3L, 0L, info = name)
    expect_identical(substring(row$junction_aa, nchar(row$junction_aa)),
                     if (built[[name]]$locus == "IGH") "W" else "F", info = name)
  }
  for (i in seq_len(nrow(cases))) {
    row <- annotated[which(annotated$sequence_id == id(cases$allele[i])), ]
    expect_true(grepl(cases$allele[i], row$j_call, fixed = TRUE), info = cases$allele[i])
  }

  # With the image's own file the same queries lose their junctions, except
  # IGLJ2A*01, whose identical sibling IGLJ2*01 is listed first by IgBLAST.
  before <- run_annotation(records, file.path(directory, "shipped"), files$shipped)
  lost <- setdiff(cases$allele, "IGLJ2A*01")
  expect_setequal(before$failed$sequence_id, c(id(lost), id("TRAJ55*01")))
  expect_true(all(before$failed$junction == "" | is.na(before$failed$junction)))
  survivor <- before$records[which(before$records$sequence_id == id("IGLJ2A*01")), ]
  expect_identical(survivor$junction, built[[id("IGLJ2A*01")]]$junction)
  expect_true(startsWith(survivor$j_call, "IGLJ2*01"))
})

test_that("productivity and SHM use the aligned sequence when a terminal stop is clipped", {
  fixture <- annotation_fixture()
  sequence <- fixture$records$sequence[fixture$records$locus == "IGH"]
  terminal <- internal <- sequence
  substr(terminal, 1, 1) <- "T"
  substr(internal, 37, 37) <- "T"
  expect_identical(substr(sequence, 1, 3), "CAG")
  expect_identical(substr(terminal, 1, 3), "TAG")
  expect_identical(substr(internal, 37, 39), "TAA")
  ids <- c("assemble1_0", "assemble2_0", "assemble3_0")
  native <- data.frame(sequence_id = ids, sequence = c(sequence, terminal, internal),
                       locus = "IGH", consensus_count = 10)
  directory <- tempfile("annotation-stops-")
  on.exit(unlink(directory, recursive = TRUE), add = TRUE)
  result <- run_annotation(native, directory)
  expect_setequal(result$records$sequence_id, ids)
  annotated <- result$records[match(ids, result$records$sequence_id), ]
  expect_identical(annotated$productive, c("T", "T", "F"))
  expect_identical(annotated$stop_codon, c("F", "F", "T"))
  expect_identical(annotated$v_sequence_start, c("1", "2", "1"))
  expect_identical(annotated$v_cigar, c("293=", "1S1N292=", "36=1X256="))
  expect_identical(substr(annotated$sequence_alignment[2], 1, 3), ".AG")
  raw <- read_pipeline_table(file.path(directory, "makedb_pass.tsv.gz"))
  raw <- raw[match(ids, raw$sequence_id), ]
  for (field in c("productive", "stop_codon", "v_sequence_start", "v_cigar")) {
    expect_identical(annotated[[field]], raw[[field]])
  }

  path <- function(name) file.path(directory, paste0(name, ".tsv"))
  expect_script_ok(run_pipeline_script("filter_airr.R", c(
    paste0("--input_file=", path("annotated")), paste0("--output_file=", path("qc")),
    "--library_id=0040", "--max_n_fraction=0.1",
    paste0("--qc_audit_file=", path("qc_audit")))))
  retained <- read_pipeline_table(path("qc"))
  expect_setequal(retained$sequence_id, ids[1:2])

  # Construct ten agreeing native votes at each query base. Measurement must
  # still exclude the clipped base from its aligned mutation opportunities.
  write_native_votes(native, file.path(directory, "final.out"))
  expect_script_ok(run_pipeline_script("build_igh_consensus_votes.R", c(
    paste0("--input_file=", path("qc")),
    paste0("--final_out_file=", file.path(directory, "final.out")),
    paste0("--output_file=", path("votes")), "--library_id=0040")))
  write_pipeline_table(data.frame(sequence_id = ids[1:2], library_id = "0040",
                                  qualified_c_call = ""), path("calls"))
  expect_script_ok(run_pipeline_script("calculate_igh_features.R", c(
    paste0("--input_file=", path("qc")), paste0("--votes_file=", path("votes")),
    paste0("--isotype_calls_file=", path("calls")),
    paste0("--output_file=", path("features")), "--library_id=0040",
    "--min_position_votes=2", paste0("--qc_audit_file=", path("shm_audit")))))
  features <- read_pipeline_table(path("features"))
  expect_setequal(features$sequence_id, ids[1:2])
  features <- features[match(ids[1:2], features$sequence_id), ]
  expect_identical(features$informative_v, c("291", "290"))
  expect_equal(as.numeric(features$shm_frequency), c(0, 0))
})

test_that("an empty family and an empty library write complete annotation outputs", {
  fixture <- annotation_fixture()
  for (population in list(fixture$records[1, ], fixture$records[FALSE, ])) {
    directory <- tempfile("annotation-empty-")
    result <- run_annotation(population, directory)
    expect_setequal(result$records$sequence_id, population$sequence_id)
    expect_length(readLines(file.path(result$scratch, "tr.fasta")), 0)
    expect_false(file.exists(file.path(result$scratch, "tr_igblast.fmt7")))
    expect_equal(nrow(result$failed), 0)
    expect_setequal(result$audit$locus, c(fixture$records$locus, "all"))
    if (!nrow(population)) {
      expect_true(all(as.numeric(result$audit$n_units_out) == 0))
      expect_true(all(result$audit$consensus_count_out %in% c("", "0")))
    }
    unlink(directory, recursive = TRUE)
  }
})

test_that("a nonempty family with no MakeDb returns retains failure evidence", {
  records <- data.frame(sequence_id = "unannotated_0", sequence = strrep("A", 360),
                        locus = "IGH", consensus_count = 19)
  directory <- tempfile("annotation-no-returns-")
  on.exit(unlink(directory, recursive = TRUE), add = TRUE)
  result <- run_annotation(records, directory)
  expect_equal(nrow(result$records), 0)
  expect_identical(result$failed$sequence_id, "unannotated_0")
  expect_identical(result$failed$receptor_family, "ig")
  row <- result$audit[result$audit$stage == "igblast_makedb" & result$audit$locus == "IGH", ]
  expect_equal(as.numeric(row[c("n_units_in", "n_units_out",
                                "consensus_count_in", "consensus_count_out")]), c(1, 0, 19, 0))
  expect_false(file.exists(file.path(result$scratch, "ig_germ-pass.tsv")))
  log <- readLines(gzfile(file.path(directory, "makedb.log.gz")))
  expect_true(any(grepl("unannotated_0", log, fixed = TRUE)))
})

test_that("shared TRAV/DV rearrangements reach the correct audit, QC and groups", {
  fixture <- annotation_fixture()
  v <- gsub(".", "", substr(fixture$alleles[["TRAV14/DV4*01"]], 1, 312), fixed = TRUE)
  delta <- paste0(v, strrep("GCC", 5), substring(fixture$alleles[["TRDJ1*01"]], 18))
  alpha <- paste0(v, strrep("GCC", 5), substring(fixture$alleles[["TRAJ1*01"]], 29))
  records <- data.frame(sequence_id = c("shared_delta_0", "shared_delta_no_c_0", "shared_alpha_0"),
    sequence = c(paste0(delta, fixture$alleles[["TRDC*01"]]), delta, alpha),
    locus = c("TRD", "TRD", "TRA"), consensus_count = c(3, 5, 7))
  directory <- tempfile("annotation-shared-v-")
  on.exit(unlink(directory, recursive = TRUE), add = TRUE)
  result <- run_annotation(records, directory)
  annotated <- result$records[match(records$sequence_id, result$records$sequence_id), ]
  expect_identical(annotated$locus, records$locus)
  expect_true(all(grepl("TRAV14/DV4", annotated$v_call, fixed = TRUE)))
  expect_true(all(annotated$junction_aa == "CAAAAAF"))
  expect_identical(annotated$c_call[2:3], c("", ""))
  expect_true(grepl("TRDC", annotated$c_call[1], fixed = TRUE))
  for (locus in c("TRA", "TRD")) {
    row <- result$audit[result$audit$stage == "igblast_locus_agreement" &
                        result$audit$locus == locus, ]
    expect_equal(as.numeric(row$n_units_in), sum(records$locus == locus))
    expect_identical(row$n_units_out, row$n_units_in)
    expect_equal(as.numeric(row$consensus_count_out), sum(records$consensus_count[records$locus == locus]))
  }
  qc <- file.path(directory, "qc.tsv")
  audit <- file.path(directory, "qc_audit.tsv")
  expect_script_ok(run_pipeline_script("filter_airr.R", c(
    paste0("--input_file=", file.path(directory, "annotated.tsv")),
    paste0("--output_file=", qc), "--library_id=0040", "--max_n_fraction=0.1",
    paste0("--qc_audit_file=", audit))))
  expect_setequal(read_pipeline_table(qc)$sequence_id, records$sequence_id)
  decisions <- read_pipeline_table(audit)
  expect_true(all(decisions$n_units_in == decisions$n_units_out))
  loci <- c("igk", "igl", "tra", "trb", "trd", "trg")
  outputs <- setNames(as.list(file.path(directory, paste0(loci, ".tsv"))), paste0("groups_", loci))
  outputs$grouping <- file.path(directory, "grouping.tsv")
  outputs$qc_audit <- file.path(directory, "grouping_audit.tsv")
  run_stage_with_options("build_clonotypes.R", c(
    list(membership_file = file.path(directory, "clonotype_members.tsv"), input_file = qc, library_id = "0040", grouping_file = outputs$grouping,
         qc_audit_file = outputs$qc_audit),
    setNames(outputs[paste0("groups_", loci)], paste0(loci, "_groups_file"))))
  expect_identical(read_pipeline_table(outputs$groups_trd)$consensus_count, "8")
  expect_identical(read_pipeline_table(outputs$groups_tra)$consensus_count, "7")
  expect_true(all(read_pipeline_table(outputs$grouping)$grouping_eligible == "TRUE"))
})

test_that("forward and reverse IgBLAST hits retain the same callable C evidence", {
  fixture <- annotation_fixture()
  sequence <- paste0(fixture$records$sequence[fixture$records$locus == "IGH"],
                      substring(fixture$alleles[["IGHG1*01"]], 1, 270))
  reverse <- chartr("ACGT", "TGCA", paste(rev(strsplit(sequence, "")[[1]]), collapse = ""))
  records <- data.frame(sequence_id = c("forward_0", "reverse_0"),
                        sequence = c(sequence, reverse), locus = "IGH", consensus_count = c(3, 5))
  directory <- tempfile("annotation-reverse-c-")
  on.exit(unlink(directory, recursive = TRUE), add = TRUE)
  result <- run_annotation(records, directory)
  annotated <- result$records[match(records$sequence_id, result$records$sequence_id), ]
  expect_identical(annotated$rev_comp, c("F", "T"))
  expect_equal(as.integer(annotated$c_informative_positions), c(270L, 270L))
  fmt7 <- readLines(file.path(result$scratch, "ig_igblast.fmt7"))
  expect_true(any(startsWith(fmt7, "C\treversed|reverse_0\t")))
  qc <- file.path(directory, "qc.tsv")
  expect_script_ok(run_pipeline_script("filter_airr.R", c(
    paste0("--input_file=", file.path(directory, "annotated.tsv")),
    paste0("--output_file=", qc), "--library_id=0040", "--max_n_fraction=0.1",
    paste0("--qc_audit_file=", file.path(directory, "qc_audit.tsv")))))
  calls <- file.path(directory, "calls.tsv")
  run_stage_with_options("qualify_isotype_calls.R", list(
    input_file = qc, library_id = "0040", min_c_positions = 40,
    qc_audit_file = file.path(directory, "isotype_audit.tsv"), output_file = calls))
  qualified <- read_pipeline_table(calls)
  expect_setequal(qualified$sequence_id, records$sequence_id)
  expect_true(all(qualified$c_call_callable == "TRUE"))
  expect_true(all(grepl("IGHG1", qualified$qualified_c_call, fixed = TRUE)))
})

# Reconstruct the operation path from the raw alignment strings, independently
# of Change-O's BTOP decoder and the pipeline's CIGAR adapter/vote mapper.
alignment_cigar <- function(fields) {
  query <- strsplit(fields[15], "", fixed = TRUE)[[1]]
  reference <- strsplit(fields[16], "", fixed = TRUE)[[1]]
  stopifnot(length(query) == length(reference))
  operations <- ifelse(query == "-", "D", ifelse(reference == "-", "I",
    ifelse(query == reference, "=", "X")))
  runs <- rle(operations)
  offsets <- c(as.integer(fields[9]), as.integer(fields[11])) - 1L
  prefix <- paste0(offsets[offsets > 0], c("S", "N")[offsets > 0], collapse = "")
  paste0(prefix, paste0(runs$lengths, runs$values, collapse = ""))
}

test_that("real V and J indels retain AIRR operations, coordinates and raw evidence", {
  fixture <- annotation_fixture()
  sequence <- fixture$records$sequence[fixture$records$locus == "IGH"]
  insert <- paste0(substr(sequence, 1, 60), "GCC", substring(sequence, 61))
  delete <- paste0(substr(sequence, 1, 60), substring(sequence, 64))
  balanced <- paste0(substr(sequence, 1, 60), "GCC", substr(sequence, 61, 120),
                      substring(sequence, 124))
  reverse <- chartr("ACGT", "TGCA", paste(rev(strsplit(insert, "")[[1]]), collapse = ""))
  j_offset <- 306 + 18
  sequences <- c(sequence, insert, delete, balanced, substring(insert, 31), reverse,
    paste0(substr(sequence, 1, j_offset), "GCC", substring(sequence, j_offset + 1)),
    paste0(substr(sequence, 1, j_offset), substring(sequence, j_offset + 4)))
  ids <- paste0(c("control", "insert", "delete", "balanced", "partial", "reverse",
                  "j_insert", "j_delete"), "_0")
  directory <- tempfile("annotation-indels-")
  on.exit(unlink(directory, recursive = TRUE), add = TRUE)
  result <- run_annotation(data.frame(sequence_id = ids, sequence = sequences,
    locus = "IGH", consensus_count = 10), directory)
  expect_setequal(result$records$sequence_id, ids)
  annotated <- result$records[match(ids, result$records$sequence_id), ]
  expect_identical(annotated$v_cigar[1:3], c("293=", "60=3I233=", "60=3D230="))
  expect_true(grepl("I", annotated$v_cigar[4]) && grepl("D", annotated$v_cigar[4]))
  expect_true(grepl("I", annotated$j_cigar[7]))
  expect_true(grepl("D", annotated$j_cigar[8]))
  expect_identical(annotated$rev_comp[6], "T")
  raw <- read_pipeline_table(file.path(directory, "makedb_pass.tsv.gz"))
  raw <- raw[match(ids, raw$sequence_id), ]
  expect_identical(raw$v_cigar[1:3], c("293=", "60=3D233=", "60=3I230="))
  hits <- strsplit(readLines(file.path(result$scratch, "ig_igblast.fmt7")), "\t", fixed = TRUE)
  seen <- character()
  for (fields in hits) {
    if (!fields[1] %in% c("V", "J")) next
    id <- sub("^reversed\\|", "", fields[2])
    key <- paste(fields[1], id)
    if (key %in% seen) next
    seen <- c(seen, key)
    column <- paste0(tolower(fields[1]), "_cigar")
    expect_identical(annotated[[column]][match(id, ids)], alignment_cigar(fields))
  }
  expect_length(seen, 2L * length(ids))
  roundtrip <- file.path(MODULE_ROOT, "test/components/annotation/integration/check_changeo_roundtrip.py")
  output <- suppressWarnings(system2("python3", shQuote(c(roundtrip, directory)), stdout = TRUE, stderr = TRUE))
  status <- attr(output, "status")
  expect_script_ok(list(script = roundtrip, output = output,
    status = if (is.null(status)) 0L else as.integer(status)))
  # Removing scratch must leave the original producer evidence available.
  unlink(result$scratch, recursive = TRUE)
  expect_equal(nrow(read_pipeline_table(file.path(directory, "makedb_pass.tsv.gz"))), length(ids))
})

# Real IgBLAST/MakeDb shapes differ from the constructed unit cases: MakeDb drops
# inserted bases from sequence_alignment, writes a deletion as gap characters
# in the observed alignment and stores a reversed query in V orientation. Each
# contig below uses its native sequence, with ten agreeing votes at each called
# base and no votes at N. A correct map retains forward positions and withholds
# votes from the opposite-oriented record. SHM still measures aligned bases.
test_that("real indels, an N run and a reversed query map votes and receive SHM", {
  fixture <- annotation_fixture()
  sequence <- fixture$records$sequence[fixture$records$locus == "IGH"]
  j_offset <- 306 + 18
  unread <- sequence
  substr(unread, 61, 63) <- "NNN"
  sequences <- c(sequence,
    paste0(substr(sequence, 1, 60), "GCC", substring(sequence, 61)),
    paste0(substr(sequence, 1, 60), substring(sequence, 64)),
    paste0(substr(sequence, 1, j_offset), "GCC", substring(sequence, j_offset + 1)),
    paste0(substr(sequence, 1, j_offset), substring(sequence, j_offset + 4)),
    unread,
    chartr("ACGT", "TGCA", paste(rev(strsplit(sequence, "")[[1]]), collapse = "")))
  # Control, V insertion, V deletion, J insertion, J deletion, an N run and the
  # reversed control. Vote parsing requires TRUST4 assembly identifiers.
  ids <- paste0("assemble", seq_along(sequences), "_0")
  native <- data.frame(sequence_id = ids, sequence = sequences, locus = "IGH",
                       consensus_count = 10)
  directory <- tempfile("annotation-indel-votes-")
  on.exit(unlink(directory, recursive = TRUE), add = TRUE)
  run_annotation(native, directory)
  path <- function(name) file.path(directory, paste0(name, ".tsv"))
  expect_script_ok(run_pipeline_script("filter_airr.R", c(
    paste0("--input_file=", path("annotated")), paste0("--output_file=", path("qc")),
    "--library_id=0040", "--max_n_fraction=0.1",
    paste0("--qc_audit_file=", path("qc_audit")))))
  retained <- read_pipeline_table(path("qc"))
  expect_setequal(retained$sequence_id, ids)
  retained <- retained[match(ids, retained$sequence_id), ]
  expect_identical(retained$v_cigar[2:3], c("60=3I233=", "60=3D230="))
  expect_identical(retained$rev_comp[7], "T")
  expect_identical(retained$sequence[7], sequence)
  expect_identical(retained$v_cigar[7], retained$v_cigar[1])

  write_native_votes(native, file.path(directory, "final.out"))
  expect_script_ok(run_pipeline_script("build_igh_consensus_votes.R", c(
    paste0("--input_file=", path("qc")),
    paste0("--final_out_file=", file.path(directory, "final.out")),
    paste0("--output_file=", path("votes")), "--library_id=0040")))
  write_pipeline_table(data.frame(sequence_id = ids, library_id = "0040",
                                  qualified_c_call = ""), path("calls"))
  expect_script_ok(run_pipeline_script("calculate_igh_features.R", c(
    paste0("--input_file=", path("qc")), paste0("--votes_file=", path("votes")),
    paste0("--isotype_calls_file=", path("calls")),
    paste0("--output_file=", path("features")), "--library_id=0040",
    "--min_position_votes=2", paste0("--qc_audit_file=", path("shm_audit")))))

  # Positions 1 to 312 observed in both the record and its D-masked germline,
  # counted independently of the pipeline. An N is not an observed base.
  window <- function(alignment) strsplit(substr(alignment, 1, 312), "", fixed = TRUE)[[1]]
  informative <- vapply(seq_along(ids), function(i) {
    sum(window(retained$sequence_alignment[i]) %in% c("A", "C", "G", "T") &
        window(retained$germline_alignment_d_mask[i]) %in% c("A", "C", "G", "T"))
  }, integer(1))
  expect_identical(informative[c(2, 4, 5, 7)], rep(informative[1], 4L))
  expect_identical(informative[c(3, 6)], rep(informative[1] - 3L, 2L))
  features <- read_pipeline_table(path("features"))
  features <- features[match(ids, features$sequence_id), ]
  expect_identical(as.integer(features$informative_v), informative)
  expect_equal(as.numeric(features$shm_frequency), rep(0, length(ids)))

  # N positions retain their mapped coordinates despite zero votes. The reversed
  # record lies on the opposite strand to its contig and keeps no mapped votes.
  mapped <- read_pipeline_table(path("votes"))
  mapped <- mapped[match(ids, mapped$sequence_id), ]
  expect_identical(as.integer(mapped$imgt_v_positions_mapped),
                   c(informative[1:5], informative[6] + 3L, 0L))
  expect_identical(mapped$junction_vote_status,
                   c(rep("mapped", 6), "opposite_contig_orientation"))
  expect_identical(as.integer(mapped$imgt_v_positions_consensus_mismatch[1:6]), rep(0L, 6L))
  expect_identical(as.integer(features$informative_supported_fwr1[1:6]),
                   as.integer(features$informative_fwr1[1:6]))
})

test_that("pinned allele identity classes preserve completeness and exact memberships", {
  directory <- tempfile("reference-identity-")
  on.exit(unlink(directory, recursive = TRUE), add = TRUE)
  files <- pipeline_reference_files(directory)
  classes <- read_germline_identity_classes(files$metadata)
  metadata <- attr(classes, "metadata")
  expect_identical(vj_tie_kind(c("IGLJ2*01,IGLJ3*01", "IGHV3-30*02,IGHV3-30-5*02",
    "IGHV4-59*09,IGHV4-61*07"), classes), c("assigned_group", "assigned_group", "ambiguous"))
  partial <- metadata[metadata$allele %in% c("IGHV4-59*09", "IGHV4-61*07"), ]
  expect_identical(partial$completeness, rep("partial in 5'", 2))
  expect_true(all(is_airr_null(partial$identity_id)))
  code <- paste(
    "import sys, json, hashlib",
    sprintf("sys.path.insert(0, %s)", encodeString(MODULE_ROOT, quote = '"')),
    "from pathlib import Path",
    "from workflow.lib.germline_reference import build_reference_metadata, reference_identity_classes",
    "files = sorted(Path('/usr/local/share/germlines/imgt/human/vdj').glob('*.fasta')) + sorted(Path('/usr/local/share/germlines/imgt/human/constant').glob('*.fasta'))",
    "classes = reference_identity_classes(build_reference_metadata(files))",
    "print(hashlib.sha256(json.dumps(classes, sort_keys=True, separators=(',', ':')).encode()).hexdigest())",
    sep = "\n")
  digest <- system2("python3", c("-c", shQuote(code)), stdout = TRUE)
  # This digest records every class member in the pinned snapshot, including
  # singleton classes. A reference change requires review of the receipt diff.
  expect_identical(digest, "ec543358edb0ce97eeafcc03f615dbe597eb266778662c82a3f26c354f5a0b43")
})

test_that("real annotation resolves TRDV1 alpha and ordinary TRAV delta joins", {
  fixture <- annotation_fixture()
  joins <- data.frame(sequence_id = c("delta_v_alpha_j_0", "alpha_v_delta_j_0"),
    v = c("TRDV1*01", "TRAV1-1*01"), j = c("TRAJ1*01", "TRDJ1*01"), j_start = c(28, 17))
  sequences <- vapply(seq_len(nrow(joins)), function(i) {
    v <- gsub(".", "", substring(fixture$alleles[[joins$v[i]]], 1, 312), fixed = TRUE)
    paste0(v, "GCCGCCGCCGCCGCC", substring(fixture$alleles[[joins$j[i]]], joins$j_start[i] + 1))
  }, character(1))
  records <- data.frame(sequence_id = joins$sequence_id, sequence = sequences,
                        locus = c("TRD", "TRA"), consensus_count = c(3, 5))
  directory <- tempfile("annotation-shared-pool-")
  on.exit(unlink(directory, recursive = TRUE), add = TRUE)
  result <- run_annotation(records, directory)
  expect_setequal(result$records$sequence_id, records$sequence_id)
  matched <- result$records[match(records$sequence_id, result$records$sequence_id), ]
  expect_identical(matched$locus, c("TRA", "TRD"))
  expect_true(all(airr_qc_pass(airr_qc_flags(matched))))
})

test_that("the pinned V/D/J hit parser convention is converted at the handoff", {
  script <- tempfile(fileext = ".py")
  on.exit(unlink(script), add = TRUE)
  writeLines(c(
    "import changeo",
    "from changeo.IO import IgBLASTReader",
    "assert changeo.__version__ == '1.3.4'",
    "reader = IgBLASTReader.__new__(IgBLASTReader)",
    "for segment in ('v', 'd', 'j'):",
    "    for btop in ('30A-C-G-20', '30-A-C-G20', '30A-C-G-20-A-C-G40'):",
    "        hit = {'segment': segment.upper(), 'bit score': '50', '% identity': '90',",
    "               'evalue': '0.001', 'BTOP': btop, 'q. start': '11', 's. start': '5'}",
    "        print(reader._parseHitScores([hit], segment)[segment + '_cigar'])"
  ), script)
  raw <- system2("python", shQuote(script), stdout = TRUE)
  expect_null(attr(raw, "status"))
  expect_identical(raw, rep(c("10S4N30=3D20=", "10S4N30=3I20=", "10S4N30=3D20=3I40="), 3))
  normalized <- unlist(lapply(seq_len(3), function(i) {
    indices <- (i - 1L) * 3L + seq_len(3)
    normalise_makedb_cigar(raw[indices], rep(11, 3), c(63, 60, 103),
      rep(400, 3), paste0("constructed_hit", indices), c("v", "d", "j")[i])
  }), use.names = FALSE)
  expect_identical(normalized,
    rep(c("10S4N30=3I20=", "10S4N30=3D20=", "10S4N30=3I20=3D40="), 3))
})
