"""Read supplied display metadata and attach it without changing library rows.

Tissue and donor context comes from a declared tissue table. Library context
comes from a declared library table. Only report preparation reads these files;
scientific tables and AIRR exports retain their core identities and measurements.
"""

import csv
import math
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Mapping, Sequence


IDENTIFIERS = {"library_id", "tissue_id", "donor_id", "scope_id"}
TIMEPOINT_TYPES = {"timepoint": "category", "timepoint_order": "number"}


@dataclass(frozen=True)
class MetadataTable:
    """Typed metadata indexed by tissue or library, with report descriptors."""

    key: str
    rows: dict[str, dict[str, Any]]
    descriptors: list[dict[str, str]]


def _kind(name: str, values: Sequence[str]) -> str:
    """Infer numbers from the whole source column; preserve labels and padded codes."""
    if name in TIMEPOINT_TYPES:
        return TIMEPOINT_TYPES[name]
    observed = [value.strip() for value in values if value.strip() not in {"", "NA"}]
    if not observed:
        return "category"
    for value in observed:
        # Leading zeroes in integer codes are meaningful labels, such as 0460.
        if re.fullmatch(r"[+-]?0[0-9]+", value):
            return "category"
        try:
            if not math.isfinite(float(value)):
                return "category"
        except ValueError:
            return "category"
    return "number"


def _value(value: str, kind: str, context: str) -> str | float | None:
    value = value.strip()
    if value == "" or (kind == "number" and value == "NA"):
        return None
    if kind == "category":
        return value
    try:
        number = float(value)
    except ValueError as error:
        raise ValueError(f"{context}: invalid numeric value {value!r}") from error
    if not math.isfinite(number):
        raise ValueError(f"{context}: numeric value must be finite or empty/NA")
    return number


def _load(path: Path, columns: Sequence[str] | None, key: str) -> MetadataTable:
    path = Path(path)
    if path.suffix.lower() not in {".csv", ".tsv"}:
        raise ValueError(f"Metadata file must be CSV or TSV: {path}")
    if columns is not None:
        if not isinstance(columns, (list, tuple)) or not columns:
            raise ValueError("Metadata columns must be a nonempty list of unique column names")
        if (any(not isinstance(name, str) or not name.strip() for name in columns)
                or len(set(columns)) != len(columns)):
            raise ValueError("Metadata columns must be a nonempty list of unique column names")
    rows = {}
    with path.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle, delimiter="\t" if path.suffix.lower() == ".tsv" else ",")
        header = reader.fieldnames or []
        if any(not name.strip() for name in header):
            raise ValueError(f"{path}: empty column name")
        if len(set(header)) != len(header):
            raise ValueError(f"{path}: duplicate column names")
        if columns is None:
            selected = [name for name in header if name not in IDENTIFIERS]
        else:
            selected = list(columns)
        if set(selected) & IDENTIFIERS:
            raise ValueError("Core identifiers cannot be selected as metadata columns")
        missing = {key, *selected} - set(header)
        if missing:
            raise ValueError(f"{path}: missing selected columns {sorted(missing)}")
        if key == "library_id" and set(selected) & set(TIMEPOINT_TYPES):
            raise ValueError("timepoint and timepoint_order must be tissue-level metadata")
        for line, row in enumerate(reader, start=2):
            if None in row or None in row.values():
                raise ValueError(f"{path}:{line}: row width differs from the header")
            identity = row[key].strip()
            if not identity:
                raise ValueError(f"{path}:{line}: empty {key}")
            if identity in rows:
                raise ValueError(f"{path}:{line}: duplicate {key} {identity!r}")
            values = {key: identity}
            for identifier in ("tissue_id", "donor_id"):
                if identifier in row and identifier != key:
                    values[identifier] = row[identifier].strip() or None
            values.update({"meta_" + name: row[name] for name in selected})
            rows[identity] = values
    descriptors = []
    for name in selected:
        column = "meta_" + name
        kind = _kind(name, [row[column] for row in rows.values()])
        descriptors.append({"column": column, "label": name, "type": kind,
                            "source_level": "tissue" if key == "tissue_id" else "library"})
        for identity, row in rows.items():
            row[column] = _value(row[column], kind, f"{path}, {key}={identity!r}, {name}")
    return MetadataTable(key, rows, descriptors)


def load_tissue_metadata(path: Path, columns: Sequence[str] | None = None) -> MetadataTable:
    """Read selected tissue fields, or all non-identifier columns when omitted.

    Optional donor IDs validate joins. The source column determines its display
    type before the scope join excludes extra tissues. Headers remain unchanged.
    """
    return _load(path, columns, "tissue_id")


def load_library_metadata(path: Path, columns: Sequence[str] | None = None) -> MetadataTable:
    """Read selected library fields, or all non-identifier fields when omitted.

    A samplesheet may be an explicit source. Empty cells become missing values.
    Category fields retain literal NA; numeric fields treat NA as missing.
    """
    return _load(path, columns, "library_id")


def _validate_timepoints(libraries: Sequence[Mapping[str, Any]]) -> None:
    """Require one-to-one visit labels and orders among this scope's known pairs.

    Missing labels/orders remain missing. Extra source tissues cannot create a
    conflict in a scope that does not contain them. Orders carry no clinical meaning.
    """
    orders = {}
    labels = {}
    for row in libraries:
        label, order = row.get("meta_timepoint"), row.get("meta_timepoint_order")
        if label is None or order is None:
            continue
        tissue = row["tissue_id"]
        if label in orders and orders[label][0] != order:
            previous, previous_tissue = orders[label]
            raise ValueError(f"Timepoint {label!r} has conflicting timepoint_order values in this scope: "
                             f"{previous:g} (tissue_id={previous_tissue!r}) and {order:g} (tissue_id={tissue!r}). "
                             "Use one order per timepoint label.")
        if order in labels and labels[order][0] != label:
            previous, previous_tissue = labels[order]
            raise ValueError(f"timepoint_order {order:g} has different timepoint labels in this scope: "
                             f"{previous!r} (tissue_id={previous_tissue!r}) and {label!r} (tissue_id={tissue!r}). "
                             "Use distinct orders for distinct timepoint labels.")
        orders.setdefault(label, (order, tissue))
        labels.setdefault(order, (label, tissue))


def join_scope_metadata(libraries: Sequence[Mapping[str, Any]],
                        tissue_metadata: MetadataTable | None = None,
                        library_metadata: MetadataTable | None = None,
                        ) -> tuple[list[dict[str, Any]], list[dict[str, str]]]:
    """Left-join metadata onto one report row per library and return descriptors.

    Keep row order, measurements and every library. Missing values stay missing;
    annotations are never filled from another tissue. Repeated libraries from
    one tissue remain separate rows. Validate supplied visit order within the scope.
    """
    tables = [("tissue_id", tissue_metadata), ("library_id", library_metadata)]
    descriptors = [descriptor.copy() for _, table in tables if table is not None
                   for descriptor in table.descriptors]
    names = [descriptor["column"] for descriptor in descriptors]
    if len(set(names)) != len(names):
        raise ValueError("Tissue and library metadata columns collide")
    for key, table in tables:
        if table is not None and table.key != key:
            raise ValueError(f"Expected metadata keyed by {key}, got {table.key}")
    seen_libraries = set()
    tissue_donors = {}
    joined = []
    for library in libraries:
        for key in ("library_id", "tissue_id", "donor_id"):
            if not isinstance(library.get(key), str) or not library[key].strip():
                raise ValueError(f"Scope rows require nonempty string {key}")
        identity, tissue, donor = (library[key] for key in ("library_id", "tissue_id", "donor_id"))
        if identity in seen_libraries:
            raise ValueError(f"Duplicate library_id in scope rows: {identity!r}")
        seen_libraries.add(identity)
        if tissue in tissue_donors and tissue_donors[tissue] != donor:
            raise ValueError(f"Tissue {tissue!r} maps to more than one donor")
        tissue_donors[tissue] = donor
        collisions = set(library) & set(names)
        if collisions:
            raise ValueError(f"Metadata would overwrite existing columns: {sorted(collisions)}")
        row = dict(library)
        for key, table in tables:
            if table is None:
                continue
            metadata = table.rows.get(library[key], {})
            for identifier in ("tissue_id", "donor_id"):
                supplied = metadata.get(identifier)
                if supplied is not None and supplied != library[identifier]:
                    raise ValueError(f"Library {identity!r}: supplied {identifier} {supplied!r} "
                                     f"does not match {library[identifier]!r}")
            for descriptor in table.descriptors:
                name = descriptor["column"]
                value = metadata.get(name)
                row[name] = value
        joined.append(row)
    _validate_timepoints(joined)
    return joined, descriptors
