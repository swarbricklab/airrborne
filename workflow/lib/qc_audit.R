# Format population counts supplied by scientific producers using base R.
# Each row states its unit, locus basis, analysis and decision status, and
# QC_AUDIT_KEY identifies it. Validate keys, counts and declared conservation
# without applying biological eligibility.

QC_AUDIT_LOCI <- c("IGH", "IGK", "IGL", "TRA", "TRB", "TRD", "TRG")
QC_AUDIT_COLUMNS <- c(
  "scope_id", "grain", "donor_id", "tissue_id", "library_id", "locus",
  "locus_basis", "analysis", "stage", "stage_index", "isotype_group", "v_region",
  "kind", "decision_status", "input_unit", "output_unit", "n_units_in",
  "n_units_out", "consensus_count_in", "consensus_count_out", "support_conserved"
)
QC_AUDIT_NUMERIC_COLUMNS <- c(
  "stage_index", "n_units_in", "n_units_out", "consensus_count_in", "consensus_count_out"
)
QC_AUDIT_KEY <- c("scope_id", "donor_id", "library_id", "locus", "locus_basis",
                  "analysis", "stage", "isotype_group", "v_region")

qc_audit_empty <- function() {
  out <- make_empty_frame(QC_AUDIT_COLUMNS)
  for (column in QC_AUDIT_NUMERIC_COLUMNS) out[[column]] <- numeric()
  out$support_conserved <- logical()
  out
}

# One row per explicitly supplied population. Checkpoints describe one unchanged
# population; eligibility compares a population with its selected subset.
# Counts are in `input_unit` and `output_unit`, and support is summed
# consensus_count. `decision_status` is applied, disabled, not_applicable,
# unavailable or unknown_history. A status other than applied is not a zero:
# the row carries whatever counts its producer supplied, or NA.
qc_audit_rows <- function(library_id, stage, n_in, n_out = n_in,
                           support_in = NA_real_, support_out = support_in,
                           locus = "IGH", analysis, stage_index, kind = "checkpoint",
                           input_unit = "rearrangement", output_unit = input_unit,
                           locus_basis = "annotated", decision_status = "applied",
                           support_conserved = !is.na(support_in),
                           scope_id = NA_character_, donor_id = NA_character_,
                           tissue_id = NA_character_, grain = "library",
                           isotype_group = NA_character_, v_region = NA_character_) {
  out <- data.frame(
    scope_id, grain, donor_id, tissue_id, library_id, locus, locus_basis, analysis,
    stage, stage_index, isotype_group, v_region, kind, decision_status,
    input_unit, output_unit, n_units_in = as.numeric(n_in),
    n_units_out = as.numeric(n_out), consensus_count_in = as.numeric(support_in),
    consensus_count_out = as.numeric(support_out), support_conserved,
    stringsAsFactors = FALSE
  )
  validate_qc_audit(out)
  out[, QC_AUDIT_COLUMNS, drop = FALSE]
}

# Sum x within each requested locus. A value outside `loci` contributes to no row,
# which is what keeps these rows a canonical-locus partition; the routing stage
# accounts for the records that fall outside it.
qc_audit_locus_sums <- function(locus, x = NULL, loci = QC_AUDIT_LOCI) {
  if (is.null(x)) x <- rep(1, length(locus))
  if (length(x) != length(locus) || anyNA(locus)) stop("Invalid audit locus vector")
  vapply(loci, function(value) sum(x[locus == value]), numeric(1), USE.NAMES = FALSE)
}

# Totals are meaningful only for a complete canonical-locus partition. The
# collector combines IGH and non-IGH partitions before adding their totals.
qc_audit_locus_rows <- function(library_id, stage, n_in, n_out = n_in,
                                 support_in = NA_real_, support_out = support_in,
                                 loci = QC_AUDIT_LOCI, total = TRUE, ...) {
  rows <- qc_audit_rows(library_id, stage, n_in, n_out, support_in, support_out,
                        locus = loci, ...)
  if (total) {
    if (!setequal(loci, QC_AUDIT_LOCI)) stop("Audit total needs all canonical loci")
    rows <- rbind(rows, qc_audit_rows(
      library_id, stage, sum(rows$n_units_in), sum(rows$n_units_out),
      sum(rows$consensus_count_in), sum(rows$consensus_count_out), locus = "all", ...
    ))
  }
  rows
}

# Rows for one population and the subset `keep` selects from it: one per locus
# in `loci`, plus an `all` total when `total` is TRUE. `support` is each
# record's read support and `locus` its locus on the row's locus basis.
qc_audit_gate_rows <- function(library_id, stage, keep, support, locus,
                                loci = QC_AUDIT_LOCI, total = TRUE, ...) {
  if (!is.logical(keep) || anyNA(keep) || length(keep) != length(support) ||
      length(locus) != length(keep)) stop("Invalid audit population or decision mask")
  qc_audit_locus_rows(
    library_id, stage, qc_audit_locus_sums(locus, loci = loci),
    qc_audit_locus_sums(locus[keep], loci = loci),
    qc_audit_locus_sums(locus, support, loci),
    qc_audit_locus_sums(locus[keep], support[keep], loci),
    loci = loci, total = total, ...
  )
}

# The flags and their conjunction are the exact decisions used by common QC.
# Individual failures overlap because each criterion uses the same input, so
# the criterion rows are not sequential losses and must not be summed.
# A maximum N fraction of 1 removes the fraction ceiling only, so low_n is
# recorded as disabled while it still fails an alignment with no gap-excluded
# positions. See airr_qc_flags() in record_qc.R.
qc_audit_qc_rows <- function(library_id, flags, pass, support, locus, max_n_fraction) {
  rows <- lapply(seq_along(flags), function(index) {
    criterion <- names(flags)[index]
    qc_audit_gate_rows(
      library_id, paste0("qc_", criterion), flags[[criterion]], support, locus,
      analysis = "common", stage_index = 10 + index, kind = "criterion",
      support_conserved = FALSE,
      decision_status = if (criterion == "low_n" && max_n_fraction == 1) "disabled" else "applied"
    )
  })
  rows[[length(rows) + 1L]] <- qc_audit_gate_rows(
    library_id, "qc_pass_total", pass, support, locus, analysis = "common",
    stage_index = 20, kind = "eligibility", support_conserved = FALSE
  )
  do.call(rbind, rows)
}

# Mapping rows come from the grouping operation, with source support unchanged.
# grouping_eligible counts common-QC records against those grouping can use.
# assembly_backbone_collapse, at IGH, IGK and IGL only, counts eligible records
# against backbone representatives. exact_sequence_collapse counts backbone
# representatives against inference representatives, whose transferred support
# equals the eligible support. inference_representatives restates that
# population as a checkpoint.
qc_audit_grouping_rows <- function(library_id, records, mapping, loci) {
  mapping <- validate_grouping_mapping(records, mapping)
  eligible <- mapping$grouping_eligible
  backbone <- mapping$is_backbone_representative
  representative <- mapping$is_inference_representative
  support <- as_read_support(records$consensus_count)
  counts <- qc_audit_locus_sums(records$locus[representative], loci = loci)
  inferred_support <- qc_audit_locus_sums(
    records$locus[representative], mapping$inference_consensus_count[representative], loci
  )
  rows <- qc_audit_gate_rows(
    library_id, "grouping_eligible", eligible, support, records$locus, loci, FALSE,
    analysis = "grouping", stage_index = 30, kind = "eligibility", support_conserved = FALSE
  )
  bcr <- intersect(loci, c("IGH", "IGK", "IGL"))
  if (length(bcr)) {
    rows <- rbind(rows, qc_audit_locus_rows(
      library_id, "assembly_backbone_collapse",
      qc_audit_locus_sums(records$locus[eligible], loci = bcr),
      qc_audit_locus_sums(records$locus[backbone], loci = bcr),
      qc_audit_locus_sums(records$locus[eligible], support[eligible], bcr),
      qc_audit_locus_sums(records$locus[eligible], support[eligible], bcr),
      loci = bcr, total = FALSE, analysis = "grouping", stage_index = 31,
      kind = "grouping", output_unit = "backbone_representative", support_conserved = TRUE
    ))
  }
  rows <- rbind(rows, qc_audit_locus_rows(
    library_id, "exact_sequence_collapse",
    qc_audit_locus_sums(records$locus[backbone], loci = loci), counts,
    qc_audit_locus_sums(records$locus[eligible], support[eligible], loci), inferred_support,
    loci = loci, total = FALSE, analysis = "grouping", stage_index = 32,
    kind = "grouping", input_unit = "backbone_representative",
    output_unit = "inference_representative", support_conserved = TRUE
  ))
  rbind(rows, qc_audit_locus_rows(
    library_id, "inference_representatives", counts, support_in = inferred_support,
    loci = loci, total = FALSE, analysis = "grouping", stage_index = 33,
    input_unit = "inference_representative", support_conserved = TRUE
  ))
}

# Check the vocabulary, key uniqueness, non-negative counts and whole unit
# counts. A checkpoint must keep its unit count, a row declaring conserved
# support must have equal defined support, and a subset in its input's unit
# must not exceed that input. Returns `rows`, invisibly.
validate_qc_audit <- function(rows) {
  require_columns(rows, QC_AUDIT_COLUMNS, "QC audit")
  if (!nrow(rows)) return(invisible(rows))
  if (any(!rows$grain %in% c("library", "donor")) ||
      any(!rows$kind %in% c("checkpoint", "criterion", "eligibility", "grouping", "diagnostic")) ||
      any(!rows$locus_basis %in% c("native", "annotated", "not_applicable")) ||
      any(!rows$decision_status %in% c("applied", "disabled", "not_applicable",
                                     "unknown_history", "unavailable"))) stop("Invalid QC audit vocabulary")
  key <- rows[, QC_AUDIT_KEY, drop = FALSE]
  key[] <- lapply(key, function(x) { x[is.na(x)] <- ""; as.character(x) })
  if (anyDuplicated(key)) stop("Duplicate QC audit row key")
  for (column in QC_AUDIT_NUMERIC_COLUMNS) {
    x <- rows[[column]]
    if (any(!is.na(x) & (!is.finite(x) | x < 0))) stop("Invalid audit count: ", column)
  }
  for (column in c("stage_index", "n_units_in", "n_units_out")) {
    x <- rows[[column]]
    if (any(!is.na(x) & x != floor(x))) stop("Non-integer audit unit count: ", column)
  }
  if (any(!rows$locus %in% c(QC_AUDIT_LOCI, "all"))) stop("Invalid QC audit locus")
  equal <- function(a, b) (is.na(a) & is.na(b)) | (!is.na(a) & !is.na(b) & a == b)
  if (any(rows$kind == "checkpoint" & !equal(rows$n_units_in, rows$n_units_out))) {
    stop("QC audit checkpoint changed unit counts")
  }
  conserved <- !is.na(rows$support_conserved) & rows$support_conserved
  if (any(conserved & (is.na(rows$consensus_count_in) |
                       !equal(rows$consensus_count_in, rows$consensus_count_out)))) {
    stop("QC audit operation did not preserve defined support")
  }
  comparable <- rows$input_unit == rows$output_unit &
    rows$kind %in% c("eligibility", "criterion", "diagnostic", "grouping")
  if (any(comparable & rows$n_units_out > rows$n_units_in, na.rm = TRUE)) {
    stop("QC audit subset exceeds its input population")
  }
  invisible(rows)
}

write_qc_audit <- function(rows, path) {
  validate_qc_audit(rows)
  write_pipeline_table(rows, path, columns = QC_AUDIT_COLUMNS)
  invisible(path)
}

read_qc_audit <- function(path) {
  rows <- read_pipeline_table(path, required_columns = QC_AUDIT_COLUMNS)
  for (column in QC_AUDIT_NUMERIC_COLUMNS) rows[[column]] <- as.numeric(rows[[column]])
  rows$support_conserved <- as.logical(rows$support_conserved)
  validate_qc_audit(rows)
  rows
}
