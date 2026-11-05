"""Bind supplied measurements with explicit ordering and bounded collection memory."""

import csv
import heapq
from collections.abc import Collection, Mapping, Sequence
from contextlib import ExitStack
from operator import itemgetter
from pathlib import Path

# Every released table that carries these identifiers leads with them, in this
# order. Donor-level and pair tables carry only the first, or the first and the
# last, and keep that order too.
IDENTITY_ORDER = ("scope_id", "library_id", "tissue_id", "donor_id")


def identity_first(columns) -> list[str]:
    """Move the identifiers present in `columns` to the front, in IDENTITY_ORDER."""
    columns = list(columns)
    return ([name for name in IDENTITY_ORDER if name in columns]
            + [name for name in columns if name not in IDENTITY_ORDER])


def read_entity_labels(path: Path) -> dict[str, str]:
    """Read the readable label of each composite V/J entity from reference metadata.

    `path` is the reference allele metadata. Returns `identity_label` keyed by
    `identity_id`, as read_vj_entity_labels() in gene_usage.R does, and stops when
    one identifier carries two labels.
    """
    labels: dict[str, str] = {}
    with Path(path).open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            identity = row["identity_id"]
            if identity in ("", "NA"):
                continue
            if labels.setdefault(identity, row["identity_label"]) != row["identity_label"]:
                raise ValueError(f"Reference metadata gives {identity} more than one identity_label: {path}")
    return labels


def entity_label(entity: str, labels: dict[str, str]) -> str:
    """Name one V/J entity as vj_entity_display_label() in gene_usage.R does.

    A composite entity takes its reference label, the sorted member genes of its
    exact sequence class; a single gene is its own label, and an empty entity
    stays empty.
    """
    return labels.get(entity, entity)


def collect_scope_metric(manifest_path: Path, inputs: list[Path], output: Path,
                         context_columns: tuple[str, ...],
                         quoting: int = csv.QUOTE_NONE,
                         sort_columns: tuple[str, ...] = (),
                         label_columns: dict[str, str] | None = None,
                         entity_labels: dict[str, str] | None = None) -> None:
    """Add the manifest context columns while streaming each retained library table.

    The identifiers lead the output in IDENTITY_ORDER, followed by the source
    columns. Rows keep their source order, unless `sort_columns` is given: the
    rows are then held in memory and stably sorted by those columns as text,
    which suits a table that combines a per-library arm with a scope arm.
    Manifest values replace any source cells with the same context column names.
    `label_columns` maps a V/J entity column to the readable label column written
    after it, named through `entity_labels` from read_entity_labels().
    Inputs must share a header; header-only tables retain that header. A
    successful collection atomically replaces the output, so published hard links
    keep their previous contents.
    """
    label_columns = label_columns or {}
    if label_columns and entity_labels is None:
        raise ValueError("Label columns need the reference entity labels")
    with manifest_path.open(newline="") as handle:
        manifest = list(csv.DictReader(handle, delimiter="\t"))
    libraries = {row["library_id"]: row for row in manifest}
    if len(libraries) != len(manifest) or not libraries:
        raise ValueError("Scope manifest must contain unique libraries")
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_suffix(output.suffix + ".tmp")
    try:
        with temporary.open("w", newline="") as destination:
            writer = None
            columns = None
            held = []
            for path in inputs:
                with path.open(newline="") as handle:
                    reader = csv.DictReader(handle, delimiter="\t")
                    if columns is None:
                        columns = reader.fieldnames
                        if not columns or "library_id" not in columns:
                            raise ValueError(f"Missing library identifier in {path}")
                        added = [name for name in context_columns if name not in columns]
                        fields = identity_first(columns + added)
                        for source, label in label_columns.items():
                            if source not in fields or label in fields:
                                raise ValueError(f"Cannot label {source} as {label} in {path}")
                            fields.insert(fields.index(source) + 1, label)
                        writer = csv.DictWriter(destination, fields, delimiter="\t",
                                                lineterminator="\n", quoting=quoting)
                        writer.writeheader()
                    elif reader.fieldnames != columns:
                        raise ValueError(f"Inconsistent library table header: {path}")
                    for row in reader:
                        library = row["library_id"]
                        if library not in libraries:
                            raise ValueError(f"Library {library} in {path} is outside the scope")
                        row.update({name: libraries[library][name]
                                    for name in context_columns})
                        row.update({label: entity_label(row[source], entity_labels)
                                    for source, label in label_columns.items()})
                        if sort_columns:
                            held.append(row)
                        else:
                            writer.writerow(row)
            if writer is None:
                raise ValueError("No library tables supplied")
            held.sort(key=lambda row: tuple(row[name] for name in sort_columns))
            writer.writerows(held)
        temporary.replace(output)
    finally:
        temporary.unlink(missing_ok=True)


def collect_donor_tables(inputs: list[Path], output: Path, keys: list[str],
                         value_orders: Mapping[str, Sequence[str]] | None = None,
                         gated_columns: Mapping[str, str] | None = None,
                         output_columns: Mapping[str, str] | None = None,
                         dropped_values: Mapping[str, Collection[str]] | None = None,
                         reconstruction_statuses: Mapping[str, str] | None = None) -> None:
    """Combine donor tables already sorted by their key columns.

    A key column named in `value_orders` sorts by the position of its value in
    that list, and a value the list omits stops the collection. Other key
    columns sort as text. `gated_columns` maps a measurement column to the
    eligibility column that must read TRUE for the output to keep its value;
    otherwise the output cell is empty, and the donor table keeps the value.
    `dropped_values` maps a column to values whose rows the output leaves out.
    `output_columns` maps each output column, in output order, to the donor
    table column it copies; without it the output keeps every donor column.
    Keys, gates and dropped values name donor table columns, and the donor
    tables keep every row and column.
    Inputs must share a header and contain unique keys across all tables. The
    output leads with its identifiers in IDENTITY_ORDER.
    Inputs with headers but no rows produce an output with the same header.
    Keep one pending row per input in memory. Replace the output only after
    all rows pass the checks, so a failure leaves the previous file intact.
    """
    ranks = {column: {value: rank for rank, value in enumerate(values)}
             for column, values in (value_orders or {}).items()}
    if not set(ranks).issubset(keys):
        raise ValueError("Ordered donor table columns must be key columns")
    gated_columns = dict(gated_columns or {})
    if set(gated_columns) & set(keys):
        raise ValueError("Gated donor table columns cannot be key columns")
    dropped = {column: set(values) for column, values in (dropped_values or {}).items()}
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_suffix(output.suffix + ".tmp")
    try:
        with ExitStack() as stack:
            readers = [csv.DictReader(stack.enter_context(path.open(newline="")),
                                      delimiter="\t") for path in inputs]
            columns = readers[0].fieldnames if readers else None
            if not columns or not set(keys).issubset(columns):
                raise ValueError("Donor tables must have the declared key columns")
            if any(reader.fieldnames != columns for reader in readers):
                raise ValueError("Donor table headers differ")
            if not {*gated_columns, *gated_columns.values()}.issubset(columns):
                raise ValueError("Donor tables must have the gated and eligibility columns")
            if not set(dropped).issubset(columns):
                raise ValueError("Donor tables must have the columns whose values drop rows")
            projection = dict(output_columns) if output_columns else {name: name for name in columns if not name.startswith("sort_")}
            if not set(projection.values()).issubset(columns):
                raise ValueError("Donor tables must have every column the output copies")

            def keyed_rows(reader, path):
                for row in reader:
                    key = []
                    for column in keys:
                        value = row[column]
                        if column in ("sort_group_support", "sort_pooled_support"):
                            value = int(value)
                        elif column in ranks:
                            if value not in ranks[column]:
                                raise ValueError(f"Donor table {path} has {column} {value!r}, "
                                                 "which its declared order omits")
                            value = ranks[column][value]
                        key.append(value)
                    yield tuple(key), row

            destination = stack.enter_context(temporary.open("w", newline=""))
            fields = identity_first(projection)
            writer = csv.DictWriter(destination, fields, delimiter="\t",
                                    lineterminator="\n", quoting=csv.QUOTE_NONE)
            writer.writeheader()
            previous = None
            for key, row in heapq.merge(*(keyed_rows(reader, path)
                                          for reader, path in zip(readers, inputs)),
                                        key=itemgetter(0)):
                if previous is not None and key <= previous:
                    raise ValueError("Unsorted or duplicate donor table key: "
                                     f"{tuple(row[column] for column in keys)}")
                previous = key
                if any(row[column] in values for column, values in dropped.items()):
                    continue
                for column, eligibility in gated_columns.items():
                    if row[eligibility] != "TRUE":
                        row[column] = ""
                if reconstruction_statuses is not None:
                    row["germline_reconstruction_status"] = reconstruction_statuses.get(
                        row["sequence_id"], row["germline_reconstruction_status"])
                writer.writerow({name: row[projection[name]] for name in fields})
        temporary.replace(output)
    finally:
        temporary.unlink(missing_ok=True)


# The header of each collected overlap product when the scope has no donor
# shard, matching OVERLAP_SCHEMAS in overlap_metrics.R with scope_id first.
OVERLAP_HEADERS = {
    'donor_overlap_pairs': (
        'scope_id',
        'donor_id',
        'locus',
        'group_unit',
        'library_id_1',
        'tissue_id_1',
        'grouped_consensus_count_1',
        'n_groups_1',
        'library_id_2',
        'tissue_id_2',
        'grouped_consensus_count_2',
        'n_groups_2',
        'pair_status',
        'min_consensus_count_per_library',
        'n_groups_shared',
        'shared_consensus_count_1',
        'shared_consensus_count_2',
        'shared_support_fraction_1',
        'shared_support_fraction_2',
        'jaccard',
        'morisita_horn',
        'bray_curtis_proportion',
        'n_shared_higher_in_1',
        'n_shared_higher_in_2',
        'n_shared_equal',
    ),
    'donor_overlap_summary': (
        'scope_id',
        'donor_id',
        'locus',
        'group_unit',
        'n_libraries_in_donor',
        'n_libraries_nonempty',
        'n_pairs_total',
        'n_pairs_available',
        'donor_mean_status',
        'n_groups_union',
        'n_groups_shared_2plus',
        'n_groups_shared_all',
        'mean_jaccard',
        'mean_morisita_horn',
    ),
    'shared_groups': (
        'scope_id',
        'library_id',
        'tissue_id',
        'donor_id',
        'locus',
        'group_unit',
        'group_key',
        'representative_junction',
        'representative_junction_aa',
        'v_gene_key',
        'v_gene_names',
        'j_gene_key',
        'j_gene_names',
        'representative_sequence_id',
        'n_libraries_in_donor',
        'n_libraries_detected',
        'consensus_count',
        'support_fraction',
    ),
}


def _rows(paths: list[Path], columns: list[str]):
    """Read each required shard once, preserving its header and textual cells."""
    for path in paths:
        with path.open(newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if reader.fieldnames != columns:
                raise ValueError(f"Inconsistent table header: {path}")
            yield from reader


def _header(paths: list[Path]) -> list[str]:
    if not paths:
        raise ValueError("No required tables supplied")
    with paths[0].open(newline="") as handle:
        columns = csv.DictReader(handle, delimiter="\t").fieldnames
    if not columns:
        raise ValueError(f"Missing table header: {paths[0]}")
    return columns


def _write_atomic(output: Path, columns: list[str], rows, quoting: int = csv.QUOTE_NONE) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_suffix(output.suffix + ".tmp")
    try:
        with temporary.open("w", newline="") as destination:
            writer = csv.DictWriter(destination, columns, delimiter="\t",
                                    lineterminator="\n", quoting=quoting)
            writer.writeheader()
            writer.writerows(rows)
        temporary.replace(output)
    finally:
        temporary.unlink(missing_ok=True)


def collect_scope_selection(inputs: list[Path], output: Path) -> None:
    """Stable library ordering retains each producer's canonical isotype order."""
    columns = _header(inputs)
    rows = list(_rows(inputs, columns))
    rows.sort(key=lambda row: row["library_id"])
    _write_atomic(output, identity_first(columns), rows)


def collect_scope_overlap(inputs: list[Path], output: Path, scope_id: str,
                          product: str) -> None:
    """Stream one overlap product in donor order, including scopes with no pairs.

    The output leads with its identifiers in IDENTITY_ORDER.
    """
    columns = _header(inputs) if inputs else list(OVERLAP_HEADERS[product])
    fields = identity_first(["scope_id"] + [name for name in columns if name != "scope_id"])

    def scoped_rows():
        for row in _rows(inputs, columns):
            row["scope_id"] = scope_id
            yield row

    _write_atomic(output, fields, scoped_rows())


def collect_scope_group_abundance(manifest_path: Path, library_inputs: list[Path],
                                  scope_input: Path, output: Path, scope_id: str) -> None:
    """Merge the per-library and scope abundance tables by library, then locus.

    `library_inputs` hold one library each, in library order, and `scope_input`
    holds every library in the scope. Rows in each argument must already be
    sorted by library_id, then locus, as text. Each library and locus must come
    from one table, and keeps that table's row order. Manifest donor and tissue
    IDs and `scope_id` are added. Memory holds one pending row per argument, and
    the output is replaced only after every row passes the checks.
    """
    with manifest_path.open(newline="") as handle:
        records = list(csv.DictReader(handle, delimiter="\t"))
    manifest = {row["library_id"]: row for row in records}
    if len(manifest) != len(records):
        raise ValueError("Scope manifest has a duplicate library_id")
    columns = _header([*library_inputs, scope_input])
    fields = identity_first(["scope_id", "donor_id", "tissue_id"] + columns)

    def keyed_rows(paths):
        for path in paths:
            for row in _rows([path], columns):
                yield (row["library_id"], row["locus"]), path, row

    def merged_rows():
        previous = previous_path = None
        for key, path, row in heapq.merge(keyed_rows(library_inputs), keyed_rows([scope_input]),
                                          key=itemgetter(0)):
            if previous is not None and key < previous:
                raise ValueError("Clone abundance rows must be sorted by library_id, then locus: "
                                 f"{path} has {key} after {previous} in {previous_path}")
            if key == previous and path != previous_path:
                raise ValueError(f"Clone abundance for {key} must come from one table: "
                                 f"found in {previous_path} and {path}")
            library = row["library_id"]
            if library not in manifest:
                raise ValueError(f"Library {library} is outside the scope")
            row.update(scope_id=scope_id, donor_id=manifest[library]["donor_id"],
                       tissue_id=manifest[library]["tissue_id"])
            previous, previous_path = key, path
            yield row

    _write_atomic(output, fields, merged_rows())
