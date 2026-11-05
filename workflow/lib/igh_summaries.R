# Build one donor's IGH tables within a scope from clone membership and the
# per-library IGH features: clones with their V/J entities and isotypes,
# library-local clone isotype rows, the library SHM summary, clone SHM detail,
# the member export and its SHM rows. Clone support is pooled across the donor;
# library summaries retain specimen context. These helpers use base R.

# The released clones.tsv, in the order the data dictionary declares it.
CLONE_FAMILY_COLUMNS <- c("scope_id", "donor_id", "clone_id", "representative_sequence_id",
                          "representative_library_id", "v_gene_key",
                          "v_gene_names", "j_gene_key", "j_gene_names",
                          "junction_length", "n_members", "n_libraries",
                          "consensus_count", "isotype_classes",
                          "has_unswitched_and_switched_members")

#' Summarise every clonal family in a scope.
#'
#' `membership` is the scope's clone membership with the IGH feature sidecar's
#' `qualified_c_call` already attached, and `scope_id` stamps each row.
#' `consensus_count` sums the members' recovered read support within the donor
#' and scope. It does not count cells or molecules.
#'
#' Two isotype summaries are reported, and they answer different questions.
#' `isotype_classes` is the sorted set of major classes the members resolve to
#' unambiguously, so it describes the recovered constant calls in the clone.
#' `has_unswitched_and_switched_members` says whether the clone holds both an
#' unambiguously unswitched member (IgM or IgD) and an unambiguously switched
#' one. Their co-occurrence does not establish a treatment-associated switch. It is
#' undefined, rather than FALSE, for a clone whose members resolve to no switch
#' state at all: absence of evidence is not evidence that the clone did not
#' switch.
#'
#' `v_gene_key` and `j_gene_key` are the one V and one J analytical entity the
#' members' `v_gene` and `j_gene` share, and `junction_length` their one junction
#' length in nucleotides, because clone inference groups within exact entities
#' and equal junction lengths; a clone whose members disagree on any of the
#' three stops the build. `v_gene_names` and `j_gene_names` name each entity
#' through `entity_labels` from read_vj_entity_labels(), so a composite reads as
#' its member genes. `n_libraries` counts the distinct libraries among the
#' members.
#'
#' Returns one row per clone in CLONE_FAMILY_COLUMNS order, sorted by donor and
#' then clone, or the released header alone when there are no members.
build_clone_family_table <- function(membership, scope_id, entity_labels = character()) {
  if (nrow(membership) == 0L) {
    return(make_empty_frame(CLONE_FAMILY_COLUMNS))
  }
  support <- as_read_support(
    membership$consensus_count, "consensus_count")
  # Extract and group the required columns once, avoiding repeated copies of
  # the wide membership table for each clone.
  donor_id <- as.character(membership$donor_id)
  library_id <- as.character(membership$library_id)
  v_gene <- as.character(membership$v_gene)
  j_gene <- as.character(membership$j_gene)
  junction_length <- as.character(membership$junction_length)
  classes <- resolve_isotype(membership$qualified_c_call, "class")
  states <- resolve_isotype(membership$qualified_c_call, "switch")
  by_lineage <- split(seq_len(nrow(membership)),
                      as.character(membership$clone_id))
  # split() names its groups in sorted order, so this is the clone of each group.
  lineage_ids <- names(by_lineage)
  over_lineages <- function(f, value) {
    vapply(by_lineage, f, value, USE.NAMES = FALSE)
  }
  shared_entity <- function(values, property) {
    over_lineages(function(rows) {
      entity <- unique(values[rows])
      if (length(entity) != 1L) {
        stop(sprintf("Clone %s has members with more than one %s: %s",
                     membership$clone_id[rows[1L]], property, paste(entity, collapse = ", ")))
      }
      entity
    }, character(1))
  }
  clone_v_gene <- shared_entity(v_gene, "V entity")
  clone_j_gene <- shared_entity(j_gene, "J entity")
  clone_junction_length <- as.integer(shared_entity(junction_length, "junction length"))

  best <- observed_representative_rows(membership$clone_id, support,
    membership$junction, membership$sequence, membership$v_gene,
    membership$j_gene, membership$sequence_id)
  best <- best[match(lineage_ids, membership$clone_id[best])]
  families <- data.frame(
    scope_id = scope_id,
    donor_id = over_lineages(function(rows) donor_id[rows[1L]], character(1)),
    clone_id = lineage_ids,
    representative_sequence_id = as.character(membership$sequence_id[best]),
    representative_library_id = as.character(membership$library_id[best]),
    v_gene_key = clone_v_gene,
    v_gene_names = vj_entity_display_label(clone_v_gene, entity_labels),
    j_gene_key = clone_j_gene,
    j_gene_names = vj_entity_display_label(clone_j_gene, entity_labels),
    junction_length = clone_junction_length,
    n_members = lengths(by_lineage, use.names = FALSE),
    n_libraries = over_lineages(function(rows) length(unique(library_id[rows])), integer(1)),
    consensus_count = over_lineages(function(rows) sum(support[rows]), numeric(1)),
    isotype_classes = over_lineages(function(rows) {
      seen <- classes[rows]
      paste(sort(unique(seen[!is.na(seen)])), collapse = ",")
    }, character(1)),
    has_unswitched_and_switched_members = over_lineages(function(rows) {
      seen <- states[rows]
      seen <- unique(seen[!is.na(seen)])
      if (length(seen) == 0L) NA else all(c("unswitched", "switched") %in% seen)
    }, logical(1)),
    stringsAsFactors = FALSE
  )
  families[order(families$donor_id, families$clone_id),
           CLONE_FAMILY_COLUMNS, drop = FALSE]
}

# The released per-library clone isotype table, in the order the data dictionary
# declares it. One row per library, clone, isotype level and isotype.
LINEAGE_ISOTYPE_COLUMNS <- c(
  "scope_id", "library_id", "tissue_id", "donor_id", "clone_id", "isotype_level",
  "isotype", "n_members", "consensus_count", "callable_consensus_count",
  "support_fraction")

#' Summarise each clone's isotype composition as one library saw it.
#'
#' `membership` carries the IGH features joined by
#' attach_igh_features_to_membership(). `manifest` supplies tissue IDs and
#' `scope_id` stamps each row.
#'
#' A clone is pooled across a donor's libraries, so the same `clone_id` can appear
#' under several `library_id`; a row describes only what that library contributed
#' to it. Each level of ISOTYPE_COMPOSITION_LEVELS resolves the qualified C call
#' through resolve_isotype(), and a row exists for each isotype with at least one
#' member resolving to it.
#'
#' `consensus_count` sums the transferred support of those members.
#' `callable_consensus_count` sums, within the library and clone, the support of
#' the members resolving at that level, and `support_fraction` divides the two.
#' Common QC requires positive support, so that denominator is positive on every
#' row. A call resolving to IgG can stay unresolved at subclass level, so the two
#' levels have separate denominators. A clone with no member resolving at a level
#' has no rows there, and group_abundance.tsv keeps its support. These calls describe
#' recovered observations and do not establish the clone's switching history.
#' SHM by isotype within a clone is in clone_shm.tsv, each member's switch state
#' is in clone_members.tsv and its regional SHM is in clone_member_shm.tsv.
#'
#' Returns the rows in LINEAGE_ISOTYPE_COLUMNS order, sorted by library and clone,
#' then by level and isotype in the order of ISOTYPE_COMPOSITION_LEVELS, or the
#' released header alone when no member resolves at either level.
build_lineage_isotype_table <- function(membership, manifest, scope_id) {
  if (nrow(membership) == 0L) {
    return(make_empty_frame(LINEAGE_ISOTYPE_COLUMNS))
  }
  # Extract the required columns once.
  support <- as_read_support(
    membership$consensus_count, "consensus_count")
  library_id <- as.character(membership$library_id)
  donor_id <- as.character(membership$donor_id)
  clone_id <- as.character(membership$clone_id)
  # The library has to enter the key because a clone spans a donor's libraries. The
  # group names are never parsed back apart (the identifiers are read off each group's
  # first row), so the separator only has to be a character that cannot appear in
  # them, which the ASCII unit separator is.
  unit <- paste(library_id, clone_id, sep = "\x1f")
  levels <- names(ISOTYPE_COMPOSITION_LEVELS)
  rows <- lapply(levels, function(level) {
    isotype <- resolve_isotype(membership$qualified_c_call, level)
    callable <- which(!is.na(isotype))
    if (!length(callable)) return(NULL)
    by_isotype <- split(callable, paste(unit[callable], isotype[callable], sep = "\x1f"))
    first <- vapply(by_isotype, function(members) members[1L], integer(1), USE.NAMES = FALSE)
    callable_support <- tapply(support[callable], unit[callable], sum)
    data.frame(
      scope_id = scope_id,
      library_id = library_id[first],
      donor_id = donor_id[first],
      clone_id = clone_id[first],
      isotype_level = level,
      isotype = isotype[first],
      n_members = lengths(by_isotype, use.names = FALSE),
      consensus_count = vapply(by_isotype, function(members) sum(support[members]),
                               numeric(1), USE.NAMES = FALSE),
      callable_consensus_count = as.numeric(callable_support[unit[first]]),
      stringsAsFactors = FALSE
    )
  })
  units <- do.call(rbind, rows)
  if (is.null(units)) return(make_empty_frame(LINEAGE_ISOTYPE_COLUMNS))
  units$support_fraction <- units$consensus_count / units$callable_consensus_count

  # The manifest is read for the tissue alone, and looked up rather than merged so
  # that the row count is a property of the membership and nothing the manifest
  # says can change it.
  units$tissue_id <- as.character(
    manifest$tissue_id[match(units$library_id, as.character(manifest$library_id))])
  level_order <- unlist(lapply(levels, function(level) {
    paste(level, ISOTYPE_COMPOSITION_LEVELS[[level]])
  }))
  units[order(units$library_id, units$clone_id,
              match(paste(units$isotype_level, units$isotype), level_order)),
        LINEAGE_ISOTYPE_COLUMNS, drop = FALSE]
}

#
# Join one donor's IGH evidence once and build its released summaries.
#
# Scope-specific membership remains the source of clonal assignments and
# transferred support. The record shard supplies the same members' sequences
# and alignments. Their shared frame receives library-level features once,
# before SHM, clone and member tables are formed. Memory use follows one
# donor's records rather than all records in the scope.

#' The donor member shard's columns, including the selected SHM regions.
#'
#' Scope collection releases the subset DONOR_METRIC_COLUMNS lists in scope.smk,
#' under its names. The shard also keeps each member's sequence, alignments,
#' qualified C call and mutation measurements, and every region's eligibility and
#' evidence ID, for the release checks and reference calculations. `locus` is IGH on
#' every row. `v_gene_names` and `j_gene_names` name the `v_gene` and `j_gene`
#' entities as vj_entity_display_label() does.
donor_member_export_columns <- function() {
  c(
    "scope_id", "donor_id", "tissue_id", "library_id", "locus", "sequence_id", "clone_id",
    "clonotype_key", "v_call", "j_call", "v_gene", "v_gene_names",
    "j_gene", "j_gene_names", "qualified_c_call", "junction", "junction_aa",
    "consensus_count", "sequence", "junction_length", "sequence_alignment",
    "germline_alignment_d_mask", "isotype_class", "isotype_subclass",
    "isotype_switch_state", member_shm_columns()
  )
}

#' Attach one donor's library features to its authoritative clone records.
#'
#' records and membership are the two projections of the same assignment
#' output. Check their keys and shared fields before using the wider records
#' frame. features may include common-QC records that entered no clone, but
#' must describe libraries belonging to this donor within the scope. Return
#' the enriched records and the donor's manifest, retaining libraries with no
#' assigned sequence so their SHM summaries still receive empty rows.
enrich_donor_igh_records <- function(records, membership, features, manifest,
                                    scope_id, donor_id, region_levels) {
  require_columns(records, CLONE_RECORD_COLUMNS, "Lineage records")
  require_columns(membership, CLONE_MEMBER_COLUMNS, "Lineage membership")
  require_columns(manifest, c("library_id", "donor_id", "tissue_id"), "Scope manifest")
  manifest <- manifest[as.character(manifest$donor_id) == donor_id,
                       c("library_id", "donor_id", "tissue_id"), drop = FALSE]
  manifest <- manifest[order(manifest$library_id), , drop = FALSE]
  if (!nrow(manifest)) stop("Scope manifest has no libraries for donor ", donor_id)
  if (any(records$scope_id != scope_id | records$donor_id != donor_id) ||
      any(membership$scope_id != scope_id | membership$donor_id != donor_id)) {
    stop("Lineage inputs contain another scope or donor")
  }
  if (anyDuplicated(records$sequence_id) || anyDuplicated(membership$sequence_id)) {
    stop("sequence_id must be unique in donor lineage inputs")
  }
  if (!setequal(records$sequence_id, membership$sequence_id)) {
    stop("Lineage records and membership have different sequence_id keys")
  }
  records <- records[match(membership$sequence_id, records$sequence_id), , drop = FALSE]
  for (column in CLONE_MEMBER_COLUMNS) {
    if (!identical(as.character(records[[column]]), as.character(membership[[column]]))) {
      stop("Lineage records and membership disagree in ", column)
    }
  }
  if (any(!as.character(features$library_id) %in% as.character(manifest$library_id))) {
    stop("IGH features contain libraries outside the donor's scope manifest")
  }
  records <- add_scope_manifest_identifiers(records, manifest)
  records <- attach_igh_features_to_membership(records, features, region_levels)
  list(records = records, manifest = manifest)
}

#' Derive six donor tables from one enriched membership frame.
#'
#' Eligibility and numerical summaries use their existing producers. Only
#' after those calculations does the member export project away the raw
#' feature counts. Donor shards follow the key order DONOR_METRIC_KEYS declares
#' in scope.smk, so scope collection can merge sorted streams without loading
#' their combined records into memory. The SHM tables order isotype groups as
#' SHM_ISOTYPE_GROUPS and V regions as shm_v_regions(), and the member SHM rows
#' order V regions the same way; other keys sort as text.
#' `entity_labels` names the V and J entities of the clone and member tables.
build_donor_igh_tables <- function(records, manifest, scope_id,
                                   min_subregion_fraction,
                                   min_germline_occupied_v, region_levels, audit_counts = FALSE,
                                   evidence = NULL, entity_labels = character()) {
  if (is.null(evidence)) {
    evidence <- choose_shm_evidence(records, records, records$sequence_id,
      min_subregion_fraction, min_germline_occupied_v, region_levels)
  }
  v_regions <- shm_v_regions(region_levels)
  shm <- build_scope_shm_table(
    records, manifest, scope_id, v_regions,
    min_subregion_fraction, min_germline_occupied_v, region_levels,
    audit_counts = audit_counts, evidence = evidence)
  audit <- attr(shm, "qc_audit")
  clone_shm <- build_clone_shm_table(records, min_subregion_fraction,
                                      min_germline_occupied_v, region_levels, evidence = evidence)
  clones <- build_clone_family_table(records, scope_id, entity_labels)
  clone_isotype <- build_lineage_isotype_table(records, manifest, scope_id)
  clone_member_shm <- build_clone_member_shm_table(records, min_subregion_fraction,
                                                   min_germline_occupied_v, region_levels,
                                                   evidence = evidence)

  records$locus <- rep("IGH", nrow(records))
  records$v_gene_names <- vj_entity_display_label(records$v_gene, entity_labels)
  records$j_gene_names <- vj_entity_display_label(records$j_gene, entity_labels)
  records$isotype_class <- resolve_isotype(records$qualified_c_call, "class")
  records$isotype_subclass <- resolve_isotype(records$qualified_c_call, "subclass")
  records$isotype_switch_state <- resolve_isotype(records$qualified_c_call, "switch")
  records <- add_member_shm_columns(records, min_subregion_fraction,
                                    min_germline_occupied_v, region_levels, evidence)
  tables <- list(
    clone_shm = clone_shm[order(clone_shm$library_id,
                                match(clone_shm$isotype_group, SHM_ISOTYPE_GROUPS),
                                match(clone_shm$v_region, v_regions),
                                clone_shm$clone_id), , drop = FALSE],
    shm = shm[order(shm$library_id, match(shm$isotype_group, SHM_ISOTYPE_GROUPS),
                    match(shm$v_region, v_regions)), , drop = FALSE],
    clones = clones[order(clones$donor_id, clones$clone_id), , drop = FALSE],
    clone_isotype = clone_isotype[order(clone_isotype$library_id,
                                           clone_isotype$clone_id), , drop = FALSE],
    members = records[order(records$library_id, records$clone_id, records$sequence_id),
                       donor_member_export_columns(), drop = FALSE],
    clone_member_shm = clone_member_shm[order(clone_member_shm$library_id,
                                              clone_member_shm$clone_id,
                                              clone_member_shm$sequence_id,
                                              match(clone_member_shm$v_region, v_regions)),
                                        , drop = FALSE]
  )
  tables$clone_member_shm_coverage <- build_member_shm_coverage(records, evidence, region_levels)
  if (audit_counts) attr(tables, "qc_audit") <- audit
  tables
}

#' Read, enrich and write the six summary shards for one donor in one scope.
#'
#' input_paths names records, membership, pooled, grouping, features and manifest;
#' output_paths names
#' qc_audit, clone_shm, shm, clones, clone_isotype, members and clone_member_shm.
#' feature_columns and region_levels come from the annotation modules at the entry
#' point, keeping this summarisation code independent of the packages that measured
#' them.
#' entity_labels names the clone and member V and J entities.
write_donor_igh_summaries <- function(input_paths, output_paths, scope_id, donor_id,
                                      min_subregion_fraction,
                                      min_germline_occupied_v,
                                      region_levels, feature_columns,
                                      entity_labels = character()) {
  records <- read_pipeline_table(input_paths$records,
                                  required_columns = CLONE_RECORD_COLUMNS)
  membership <- read_pipeline_table(input_paths$membership,
                                     required_columns = CLONE_MEMBER_COLUMNS)
  features <- bind_frames_on_column_union(
    lapply(input_paths$features, function(path) {
      read_pipeline_table(path, required_columns = feature_columns,
                          context = paste("IGH feature shard", path))
    }), columns = feature_columns)
  manifest <- read_pipeline_table(input_paths$manifest)
  enriched <- enrich_donor_igh_records(records, membership, features, manifest,
                                       scope_id, donor_id, region_levels)
  source <- read_pipeline_table(input_paths$pooled)
  mapping <- validate_grouping_mapping(source, read_pipeline_table(input_paths$grouping))
  candidates <- attach_igh_features_to_membership(
    source[mapping$is_backbone_representative, , drop = FALSE], features, region_levels)
  evidence <- choose_shm_evidence(enriched$records, candidates,
    mapping$inference_representative_sequence_id[mapping$is_backbone_representative],
    min_subregion_fraction, min_germline_occupied_v, region_levels)
  log_info("Scope %s donor %s: %d members over %d libraries; joined %d feature rows",
           scope_id, donor_id, nrow(enriched$records), nrow(enriched$manifest),
           nrow(features))
  rm(records, membership, features, manifest, source, mapping, candidates)
  tables <- build_donor_igh_tables(
    enriched$records, enriched$manifest, scope_id, min_subregion_fraction,
    min_germline_occupied_v, region_levels, audit_counts = TRUE, evidence = evidence,
    entity_labels = entity_labels)
  write_qc_audit(attr(tables, "qc_audit"), output_paths$qc_audit)
  for (name in names(tables)) {
    table <- order_donor_detail(tables[[name]], name, enriched$records, region_levels)
    write_pipeline_table(table, output_paths[[name]])
    log_info("Wrote %d %s rows to %s", nrow(tables[[name]]), name, output_paths[[name]])
  }
  for (v_region in shm_v_regions(region_levels)) {
    selected <- tables$shm$isotype_group == "all_IGH" & tables$shm$v_region == v_region
    log_info("Scope %s donor %s: %s summarises %d of %d measured lineage members",
             scope_id, donor_id, v_region, sum(tables$shm$n_members_eligible[selected]),
             sum(tables$shm$n_members_measured[selected]))
  }
  invisible(NULL)
}

# Attach scope manifest context to scientific summary rows.

#' Attach the manifest identifiers a scope table reports its rows under.
#'
#' `records` are the scope's rows, carrying `library_id`; `manifest` is the scope
#' manifest; `columns` names the manifest columns to attach, which is `tissue_id`
#' for a table whose shards already carry the donor and both identifiers for one
#' whose shards carry neither. Returns `records` with those columns. One manifest
#' row per library, so this is a lookup rather than a merge and cannot change the
#' row count. A library the manifest does not list is a scoping error rather than
#' a missing value, so it stops.
add_scope_manifest_identifiers <- function(records, manifest,
                                           columns = "tissue_id") {
  require_columns(manifest, c("library_id", columns), "Scope manifest")
  if (anyDuplicated(manifest$library_id)) stop("scope manifest has a duplicate library_id")
  index <- match(records$library_id, manifest$library_id)
  if (anyNA(index)) {
    stop("Scope records contain librar(ies) absent from the scope manifest: ",
         paste(unique(records$library_id[is.na(index)]), collapse = ", "))
  }
  for (column in columns) {
    records[[column]] <- as.character(manifest[[column]][index])
  }
  records
}

# Add internal sort keys to donor shards. Scope collection removes these keys
# after merging, so presentation never requires a scope abundance dependency.
order_donor_detail <- function(table, name, records, region_levels) {
  detail <- c("clone_shm", "clone_isotype", "members", "clone_member_shm",
              "clone_member_shm_coverage")
  if (name == "clones") {
    table$sort_pooled_support <- -as.numeric(table$consensus_count)
    return(table[order(table$donor_id, table$sort_pooled_support, table$clone_id), , drop = FALSE])
  }
  if (!name %in% detail) return(table)
  groups <- split(seq_len(nrow(records)), paste(records$library_id, records$clone_id, sep = "\x1f"))
  support <- vapply(groups, function(ix) sum(as_read_support(records$consensus_count[ix])), numeric(1))
  junction <- vapply(groups, function(ix) min(as.character(records$junction[ix])), character(1))
  key <- paste(table$library_id, table$clone_id, sep = "\x1f")
  table$sort_group_support <- -unname(support[key])
  table$sort_group_junction <- unname(junction[key])
  keys <- list(table$library_id, table$sort_group_support, table$sort_group_junction, table$clone_id)
  if ("isotype_group" %in% names(table)) keys <- c(keys, list(match(table$isotype_group, SHM_ISOTYPE_GROUPS)))
  if ("isotype_level" %in% names(table)) keys <- c(keys, list(match(table$isotype_level, c("class", "subclass")),
    match(table$isotype, c("IgM", "IgD", "IgG", "IgG1", "IgG2", "IgG3", "IgG4", "IgA", "IgA1", "IgA2", "IgE"))))
  if ("sequence_id" %in% names(table)) keys <- c(keys, list(table$sequence_id))
  if ("v_region" %in% names(table)) keys <- c(keys, list(match(table$v_region, shm_v_regions(region_levels))))
  table[do.call(order, keys), , drop = FALSE]
}
