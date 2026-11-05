# Compare two donor libraries at one locus using supplied clone/clonotype keys.
# Compare detection, redistribution and the direction of differences among
# shared groups. Differences are endpoint 2 against endpoint 1; chronology is
# supplied separately. Retain empty schemas and complete support denominators.
# These helpers use base R.

# The clone-level description of a shared group, taken from its best-supported
# observation, in released order.
SHARED_CLONE_DESCRIPTOR_COLUMNS <- c(
  "representative_junction", "representative_junction_aa", "v_gene_key",
  "v_gene_names", "j_gene_key", "j_gene_names", "representative_sequence_id"
)

# The released columns of the three donor tables, keyed by product name. The
# pair table leaves out every value that follows exactly from these columns or
# from diversity.tsv: the overlap status, the Simpson concentration difference,
# the shared-group fractions and balance and the support only at one endpoint.
OVERLAP_SCHEMAS <- list(
  donor_overlap_pairs = c(
    "donor_id", "locus", "group_unit",
    "library_id_1", "tissue_id_1", "grouped_consensus_count_1", "n_groups_1",
    "library_id_2", "tissue_id_2", "grouped_consensus_count_2", "n_groups_2",
    "pair_status", "min_consensus_count_per_library",
    "n_groups_shared", "shared_consensus_count_1", "shared_consensus_count_2",
    "shared_support_fraction_1", "shared_support_fraction_2",
    "jaccard", "morisita_horn", "bray_curtis_proportion",
    "n_shared_higher_in_1", "n_shared_higher_in_2", "n_shared_equal"
  ),
  donor_overlap_summary = c(
    "donor_id", "locus", "group_unit", "n_libraries_in_donor", "n_libraries_nonempty",
    "n_pairs_total", "n_pairs_available", "donor_mean_status",
    "n_groups_union",
    "n_groups_shared_2plus", "n_groups_shared_all",
    "mean_jaccard",
    "mean_morisita_horn"
  ),
  shared_groups = c(
    "library_id", "tissue_id", "donor_id", "locus", "group_unit", "group_key",
    SHARED_CLONE_DESCRIPTOR_COLUMNS,
    "n_libraries_in_donor", "n_libraries_detected",
    "consensus_count", "support_fraction"
  )
)

# Return the three overlap tables with no rows and their released column types,
# for a donor that contributes no comparable pair.
make_empty_overlap_tables <- function() {
  list(
    donor_overlap_pairs = data.frame(
      donor_id = character(), locus = character(), group_unit = character(),
      library_id_1 = character(), tissue_id_1 = character(),
      grouped_consensus_count_1 = numeric(), n_groups_1 = integer(),
      library_id_2 = character(), tissue_id_2 = character(),
      grouped_consensus_count_2 = numeric(), n_groups_2 = integer(),
      pair_status = character(), min_consensus_count_per_library = numeric(),
      n_groups_shared = integer(),
      shared_consensus_count_1 = numeric(), shared_consensus_count_2 = numeric(),
      shared_support_fraction_1 = numeric(), shared_support_fraction_2 = numeric(),
      jaccard = numeric(), morisita_horn = numeric(), bray_curtis_proportion = numeric(),
      n_shared_higher_in_1 = integer(), n_shared_higher_in_2 = integer(),
      n_shared_equal = integer(),
      stringsAsFactors = FALSE
    ),
    donor_overlap_summary = data.frame(
      donor_id = character(), locus = character(), group_unit = character(),
      n_libraries_in_donor = integer(), n_libraries_nonempty = integer(),
      n_pairs_total = integer(), n_pairs_available = integer(),
      donor_mean_status = character(),
      n_groups_union = integer(), n_groups_shared_2plus = integer(),
      n_groups_shared_all = integer(),
      mean_jaccard = numeric(),
      mean_morisita_horn = numeric(),
      stringsAsFactors = FALSE
    ),
    shared_groups = data.frame(
      library_id = character(), tissue_id = character(),
      donor_id = character(), locus = character(),
      group_unit = character(), group_key = character(),
      representative_junction = character(), representative_junction_aa = character(),
      v_gene_key = character(), v_gene_names = character(),
      j_gene_key = character(), j_gene_names = character(),
      representative_sequence_id = character(),
      n_libraries_in_donor = integer(), n_libraries_detected = integer(),
      consensus_count = numeric(), support_fraction = numeric(),
      stringsAsFactors = FALSE
    )
  )
}

#' Jaccard similarity of two libraries' clone detection.
#'
#' Jaccard 1912 (New Phytol 11:37-50, doi:10.1111/j.1469-8137.1912.tb05611.x):
#' |A intersect B| / |A union B|, on clone presence/absence. The
#' arguments are the two presence masks over one shared, aligned key vector, so
#' the caller decides what a clone is and this only counts. Read support does not
#' enter: a clone detected once and a clone detected a thousand times weigh the
#' same. Returns one numeric, or NA when neither library detected anything and
#' the union is therefore empty.
#'
#' This is the same quantity as immunarch::repOverlap(.method = "jaccard"), which
#' computes |A intersect B| / (|A| + |B| - |A intersect B|) on the clone key
#' column. immunarch is installed neither in the pinned Immcantation image nor in
#' the host R the unit suite uses, so the identity is not asserted in process: it
#' was measured once in a separate immunarch image and the numbers committed to
#' test/data/expected/overlap_reference.tsv, which the unit suite checks against.

jaccard_index <- function(present_first, present_second) {
  n_groups_shared <- sum(present_first & present_second)
  n_union <- sum(present_first | present_second)
  if (n_union) n_groups_shared / n_union else NA_real_
}

#' Morisita-Horn similarity of two libraries' relative read support.
#'
#' Horn's 1966 proportion form, 2*sum(p*q) / (sum(p^2) + sum(q^2)), where p and q
#' are within-sample proportions over one shared, aligned key vector. Unlike
#' Jaccard this weighs each clone by how much read support it carries, so two
#' libraries dominated by the same clone score highly even when their rare clones
#' differ. Returns one numeric, or NA when both proportion vectors are empty or
#' all zero and the denominator is therefore zero.
#' calculate_repertoire_pair_overlap() calls this only for a pair with positive
#' support at both endpoints, so that branch does not arise in the pipeline.
#'   Horn 1966, Am Nat 100:419-424, doi:10.1086/282436
#'
#' This is the same quantity as immunarch::repOverlap(.method = "morisita"), which
#' computes it on raw counts as 2*sum(a*b) / ((sum(a^2)/X^2 + sum(b^2)/Y^2) * X * Y)
#' for X = sum(a), Y = sum(b): the same expression with the normalisation moved. See
#' jaccard_index() above for how that identity is checked without immunarch
#' installed.
morisita_horn_index <- function(proportions_first, proportions_second) {
  denominator <- simpson_concentration(proportions_first) +
    simpson_concentration(proportions_second)
  if (denominator > 0) {
    2 * sum(proportions_first * proportions_second) / denominator
  } else {
    NA_real_
  }
}

#' Reduce one library's records to the units that can be shared.
#'
#' `group_keys` identify inferred IGH clones or exact non-IGH clonotypes across
#' libraries; `consensus_count` is the record's retained read support. At IGH that
#' is the support already transferred to the retained representative, so nothing
#' here may transfer it a second time. A record is outside the shared
#' repertoire when it carries no key or no positive support, which is what makes
#' a library contribute to a locus exactly when it has at least one
#' positive-support clone. Support is then summed per key, because a library may
#' report one clone on several rows.
#'
#' Returns a frame of `group_key` and `consensus_count` sorted by key, with no rows
#' when nothing survived.
build_sharing_units <- function(group_keys, consensus_count) {
  units <- data.frame(
    group_key = as.character(group_keys),
    consensus_count = consensus_count,
    stringsAsFactors = FALSE
  )
  units <- units[!is.na(units$group_key) & units$group_key != "" &
                   units$consensus_count > 0, , drop = FALSE]
  if (!nrow(units)) {
    return(data.frame(group_key = character(0), consensus_count = numeric(0),
                      stringsAsFactors = FALSE))
  }
  aggregate(consensus_count ~ group_key, units, sum)
}

#' Compare one pair of libraries at one locus.
#'
#' `first` and `second` are sharing-unit frames from build_sharing_units(). Both
#' are widened onto the union of their keys, with an absent key scored as zero
#' support, so every measurement uses the original endpoint denominators.
#'
#' The unit is a receptor unit (an IGH clone or a non-IGH clonotype) with
#' positive read support at an endpoint. p and q are each unit's share of the
#' first and second endpoint's total support at this locus. A shared unit has
#' positive support at both endpoints.
#'   shared_support_fraction_k    support in shared units over endpoint k's total
#'   n_shared_higher_in_1, n_shared_higher_in_2, n_shared_equal
#'                                shared units whose share is higher at endpoint
#'                                1, higher at endpoint 2 or equal
#'
#' Bray-Curtis (Bray and Curtis 1957, Ecol Monogr 27:325-349,
#' doi:10.2307/1942268) on proportions is total variation: half the absolute
#' difference over the complete union. Shared-group counts compare unrounded
#' proportions and give every shared group one vote. Exact equality means an
#' equal observed share; no effect-size threshold or pseudocount is used. These
#' are descriptive support measurements of endpoint 2 against endpoint 1, which
#' carry no chronology, rather than cell counts. A pair minimum flags values
#' without masking them. An empty endpoint makes comparisons unavailable, and
#' pair_status names it. The shared-group fractions, the support only at one
#' endpoint (1 - shared_support_fraction_k) and the Simpson concentration
#' difference (from each endpoint's diversity row) follow exactly from the
#' released values, so they are not returned.
calculate_repertoire_pair_overlap <- function(first, second,
                                              minimum_consensus_count = 0) {
  keys <- union(first$group_key, second$group_key)
  a <- setNames(first$consensus_count, first$group_key)[keys]
  b <- setNames(second$consensus_count, second$group_key)[keys]
  a[is.na(a)] <- 0
  b[is.na(b)] <- 0
  present_first <- a > 0
  present_second <- b > 0
  shared <- present_first & present_second
  total_a <- sum(a)
  total_b <- sum(b)
  status <- if (total_a == 0 && total_b == 0) "empty_both" else if (
    total_a == 0) "empty_first" else if (total_b == 0) "empty_second" else "available"
  pair_statuses <- vapply(c(total_a, total_b), repertoire_support_status,
                         character(1), minimum_consensus_count)
  pair_status <- if (status != "available") status else if (
    any(pair_statuses == "below_minimum_support")) "below_minimum_support" else "available"
  p <- if (total_a > 0) a / total_a else rep(NA_real_, length(a))
  q <- if (total_b > 0) b / total_b else rep(NA_real_, length(b))
  difference <- q[shared] - p[shared]
  list(
    pair_status = pair_status,
    min_consensus_count_per_library = minimum_consensus_count,
    n_groups_shared = sum(shared),
    shared_consensus_count_1 = sum(a[shared]),
    shared_consensus_count_2 = sum(b[shared]),
    shared_support_fraction_1 = if (total_a > 0) sum(a[shared]) / total_a else NA_real_,
    shared_support_fraction_2 = if (total_b > 0) sum(b[shared]) / total_b else NA_real_,
    jaccard = if (status == "available") jaccard_index(present_first, present_second) else NA_real_,
    morisita_horn = if (status == "available") morisita_horn_index(p, q) else NA_real_,
    bray_curtis_proportion = if (status == "available") sum(abs(q - p)) / 2 else NA_real_,
    n_shared_higher_in_1 = if (status == "available") sum(difference < 0) else NA_integer_,
    n_shared_higher_in_2 = if (status == "available") sum(difference > 0) else NA_integer_,
    n_shared_equal = if (status == "available") sum(difference == 0) else NA_integer_
  )
}

#' Build the released pairwise overlap rows for one donor and locus.
#'
#' `repertoires` is a named list of sharing-unit frames, keyed by library_id and
#' holding only the libraries that contributed; `group_unit` names what a key
#' identifies; `tissue_by_sample` maps a library_id to its tissue. `libraries`
#' names the libraries to compare, and any of them that
#' is absent from `repertoires` is compared as an empty repertoire. The production
#' caller passes all of the donor's in-scope libraries, so every unordered pair of
#' them gets one row and a donor with n of them gets n*(n-1)/2 rows. A pair with an
#' empty endpoint carries that `pair_status` and no similarity value.
#'
#' Returns a frame with the released donor_overlap_pairs columns, empty when
#' `libraries` names fewer than two libraries and there is therefore no pair to
#' compare.
build_pairwise_overlap_rows <- function(repertoires, donor_id, locus, group_unit,
                                        tissue_by_sample,
                                        libraries = names(repertoires),
                                        minimum_consensus_count = 0) {
  for (library in setdiff(libraries, names(repertoires))) {
    repertoires[[library]] <- build_sharing_units(character(), numeric())
  }
  repertoires <- repertoires[libraries]
  if (length(repertoires) < 2L) return(make_empty_overlap_tables()$donor_overlap_pairs)
  sample_pairs <- combn(names(repertoires), 2L, simplify = FALSE)
  pair_rows <- lapply(sample_pairs, function(pair) {
    s1 <- pair[1]
    s2 <- pair[2]
    values <- calculate_repertoire_pair_overlap(repertoires[[s1]],
                                                repertoires[[s2]],
                                                minimum_consensus_count)
    data.frame(
      donor_id = donor_id, locus = locus, group_unit = group_unit,
      library_id_1 = s1, tissue_id_1 = unname(tissue_by_sample[s1]),
      library_id_2 = s2, tissue_id_2 = unname(tissue_by_sample[s2]),
      n_groups_1 = nrow(repertoires[[s1]]),
      n_groups_2 = nrow(repertoires[[s2]]),
      grouped_consensus_count_1 = sum(repertoires[[s1]]$consensus_count),
      grouped_consensus_count_2 = sum(repertoires[[s2]]$consensus_count),
      values,
      stringsAsFactors = FALSE
    )[, OVERLAP_SCHEMAS$donor_overlap_pairs, drop = FALSE]
  })
  do.call(rbind, pair_rows)
}

#' Summarise how a donor's clones are spread across their libraries at one locus.
#'
#' `repertoires` is the same named list build_pairwise_overlap_rows() consumed,
#' `pair_rows` its output, `group_unit` what a key identifies, and
#' `n_libraries_in_donor` how many libraries the donor has in scope. A library is
#' nonempty at a locus when it has at least one positive-support unit, which is
#' when `repertoires` holds it. Occupancy is how many nonempty libraries each
#' clone key appears in, which gives the union, the clones seen in two or more,
#' and the clones seen in all of the donor's libraries. n_pairs_total counts every
#' unordered pair of the donor's libraries, and n_pairs_available those whose
#' pair_status is available.
#'
#' Returns one row with the released donor_overlap_summary columns.
summarise_donor_locus_occupancy <- function(repertoires, pair_rows, donor_id,
                                            locus, group_unit, n_libraries_in_donor) {
  all_keys <- unlist(lapply(repertoires, `[[`, "group_key"), use.names = FALSE)
  occupancy <- table(all_keys)

  # A library drops out of a locus when none of its records survived
  # build_sharing_units(). The two donor-level means are published only when no
  # library dropped out, because a mean over whichever pairs happen to be
  # comparable is a mean over a different, support-selected set of libraries at
  # every locus, and would compare donors on unequal terms. This is a scientific
  # inclusion rule, not a convenience: it is why a donor can have pairwise rows
  # and still no mean. Once every library of the donor is nonempty, every pair
  # has a non-empty union and positive read support on both sides, so neither
  # index can be missing.
  every_library_contributed <- length(repertoires) == n_libraries_in_donor &&
    nrow(pair_rows) > 0L

  data.frame(
    donor_id = donor_id,
    locus = locus,
    group_unit = group_unit,
    n_libraries_in_donor = n_libraries_in_donor,
    n_libraries_nonempty = length(repertoires),
    n_pairs_total = as.integer(choose(n_libraries_in_donor, 2L)),
    n_pairs_available = sum(pair_rows$pair_status == "available"),
    donor_mean_status = if (!nrow(pair_rows)) "no_pairs" else
      if (!every_library_contributed) "empty_library" else
      if (any(pair_rows$pair_status != "available")) "available_below_minimum_support" else "available",
    n_groups_union = length(occupancy),
    n_groups_shared_2plus = sum(occupancy >= 2L),
    n_groups_shared_all = sum(occupancy == n_libraries_in_donor),
    mean_jaccard =
      if (every_library_contributed) mean(pair_rows$jaccard) else NA_real_,
    mean_morisita_horn =
      if (every_library_contributed) mean(pair_rows$morisita_horn) else NA_real_,
    stringsAsFactors = FALSE
  )
}

# ------------------------------------------------------------ which clones ---

#' Describe each clone by its best-supported observation.
#'
#' `records` are every one of a donor's records at one locus, carrying `group_key`,
#' the `consensus_count` that observation was seen at, the full `sequence` and a
#' `record_id` that break a tie deterministically, and the
#' SHARED_CLONE_DESCRIPTOR_COLUMNS. Returns one row
#' per key (`group_key` followed by those columns) in key order. Equal support is
#' broken by the descriptor content before `record_id`, because IGH record IDs
#' follow TRUST4 numbering, which changes between runs.
#'
#' One clone can have multiple observations across libraries and sequences.
#' Inferred IGH clone members can differ in their junctions and contigs.
#' Rather than build a consensus that was never assembled, the row shows the
#' observation the most reads stood behind and says so. For the non-IGH loci the
#' junction is identical across the members by construction, so only the contig
#' and the displayed call sets can vary.
summarise_shared_clone_descriptors <- function(records) {
  columns <- c("group_key", SHARED_CLONE_DESCRIPTOR_COLUMNS)
  if (nrow(records) == 0L) {
    return(make_empty_overlap_tables()$shared_groups[, columns, drop = FALSE])
  }
  best <- observed_representative_rows(records$group_key, records$consensus_count,
    records$representative_junction, records$sequence, records$v_gene_key,
    records$j_gene_key, records$record_id)
  descriptors <- records[best, columns, drop = FALSE]
  descriptors[order(descriptors$group_key), , drop = FALSE]
}

#' Build the released rows naming the clones a donor's libraries have in common.
#'
#' `repertoires` is the same named list of sharing-unit frames the two summaries
#' above read; `descriptors` comes from summarise_shared_clone_descriptors();
#' `group_unit` names what a key identifies; `library_ids` are all of the donor's
#' in-scope libraries and `tissue_by_sample` maps each to its tissue.
#'
#' A clone is reported when two or more contributing libraries detected it, which
#' is the same occupancy that `n_groups_shared_2plus` counts, so the two tables
#' agree by construction rather than by argument. Each reported clone then gets
#' one row for every one of the donor's libraries, including those that did not
#' detect it, which carry zero support. The reader is comparing what each library
#' recovered, so the libraries that did not recover the clone have to be visible
#' beside those that did. Leaving those rows out would make non-detection
#' indistinguishable from a library that was never sequenced.
#'
#' `support_fraction` is the clone's share of that library's read support at this
#' locus, the same denominator grouped_consensus_count in diversity.tsv
#' uses. It is NA, not zero, for a library that contributed no clone at all at
#' this locus: there is no repertoire for the clone to be a fraction of.
#'
#' Returns a frame with the released shared-clone columns, empty when no clone
#' was detected twice. Rows run by the clone's total support across the donor's
#' libraries, largest first, then group_key, then library_id, so every library's
#' view of one clone sits together.
build_donor_locus_shared_clone_rows <- function(repertoires, descriptors,
                                                donor_id, locus, group_unit,
                                                library_ids, tissue_by_sample) {
  empty <- make_empty_overlap_tables()$shared_groups
  all_keys <- unlist(lapply(repertoires, `[[`, "group_key"), use.names = FALSE)
  occupancy <- table(all_keys)
  shared <- sort(names(occupancy)[occupancy >= 2L])
  if (!length(shared)) return(empty)

  # One lookup per library, reused across every shared clone, so a donor with
  # thousands of shared clones does not rescan its repertoires per clone.
  support_by_sample <- lapply(repertoires, function(units) {
    setNames(units$consensus_count, units$group_key)
  })
  totals <- vapply(repertoires, function(units) sum(units$consensus_count), numeric(1))

  described <- descriptors[match(shared, descriptors$group_key), , drop = FALSE]
  rows <- lapply(library_ids, function(library_id) {
    support <- unname(support_by_sample[[library_id]][shared])
    if (is.null(support)) support <- rep(NA_real_, length(shared))
    support[is.na(support)] <- 0
    total <- if (library_id %in% names(totals)) totals[[library_id]] else 0

    data.frame(
      library_id = library_id,
      tissue_id = unname(tissue_by_sample[library_id]),
      donor_id = donor_id,
      locus = locus,
      group_unit = group_unit,
      group_key = shared,
      described[, SHARED_CLONE_DESCRIPTOR_COLUMNS, drop = FALSE],
      n_libraries_in_donor = length(library_ids),
      n_libraries_detected = as.integer(occupancy[shared]),
      consensus_count = support,
      support_fraction = if (total > 0) support / total else NA_real_,
      stringsAsFactors = FALSE
    )
  })

  out <- do.call(rbind, rows)
  clone_total <- tapply(out$consensus_count, out$group_key, sum)[out$group_key]
  out <- out[order(-clone_total, out$group_key, out$library_id),
             OVERLAP_SCHEMAS$shared_groups, drop = FALSE]
  rownames(out) <- NULL
  out
}

#' One locus's records in the shape summarise_shared_clone_descriptors() reads.
#' The two arms name their key, their tie-break identifier, their sequence
#' column and their V and J entity columns differently; canonicalising here is
#' what lets one released table carry both. Both arms describe V and J by the
#' analytical entity, read from `gene_columns` (v_gene and j_gene in the IGH
#' lineage records, v_gene_key and j_gene_key in the clonotype tables), so a
#' descriptor means the same at every locus, and vj_entity_display_label()
#' names each entity through `entity_labels`. `representative_column` holds the
#' library-qualified sequence_id of the record the sequence comes from: the
#' record itself at IGH, and the clonotype's representative elsewhere. The
#' `sequence` ranks tied observations and is not released.
as_descriptor_records <- function(df, key_column, id_column,
                                  sequence_column = "sequence",
                                  representative_column = id_column,
                                  gene_columns = c("v_gene", "j_gene"),
                                  entity_labels = character()) {
  v_gene <- as.character(df[[gene_columns[1]]])
  j_gene <- as.character(df[[gene_columns[2]]])
  data.frame(
    group_key = as.character(df[[key_column]]),
    consensus_count = as_read_support(df$consensus_count, "consensus_count"),
    record_id = as.character(df[[id_column]]),
    sequence = as.character(df[[sequence_column]]),
    representative_junction = as.character(df$junction),
    representative_junction_aa = as.character(df$junction_aa),
    v_gene_key = v_gene,
    v_gene_names = vj_entity_display_label(v_gene, entity_labels),
    j_gene_key = j_gene,
    j_gene_names = vj_entity_display_label(j_gene, entity_labels),
    representative_sequence_id = as.character(df[[representative_column]]),
    stringsAsFactors = FALSE
  )
}
