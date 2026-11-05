"""Assemble the curated library sheet from authoritative scope measurements.

One row is one library within one scope. Start with the complete manifest and
copy measurements by their declared keys, including unavailable values. The
ratios formed here describe whole-V SHM eligibility from supplied counts. All
scientific estimates and availability reasons come upstream.
"""

import csv
import math
from pathlib import Path
from typing import Any

IDENTITY_COLUMNS = ("scope_id", "library_id", "tissue_id", "donor_id")
LOCI = ("igh", "igk", "igl", "tra", "trb", "trd", "trg")
YIELD_COLUMNS = tuple(f"{locus}_yield_per_million_read_pairs" for locus in LOCI)
# Copied from qc.tsv for IGH and TRB.
QC_FIELDS = ("common_qc_fraction_of_annotated", "grouping_eligible_fraction_of_common_qc")
# Copied from diversity.tsv for IGH and TRB, support first and each status before
# the estimate it qualifies.
DIVERSITY_FIELDS = (
    "grouped_consensus_count", "richness", "diversity_status",
    "top_5_support_fraction", "rarefaction_target", "rarefaction_status",
    "rarefied_inverse_simpson",
)
CLASSES = ("IgM", "IgD", "IgG", "IgA", "IgE")
SHM_GROUPS = ("all_IGH", "IgM", "IgG", "IgA")
SHM_FIELDS = (
    "n_clones_eligible", "mean_clone_median_shm", "fraction_clones_zero_median_shm",
    "n_members_eligible", "n_members_measured",
)
# Copied from the whole-V rows of shm.tsv, so each headline name ends in whole_v.
# n_clones_eligible there counts clones with an SHM-eligible member, which the
# headline name says because the headline also holds grouping-eligible counts.
SHM_NAMES = {"n_clones_eligible": "n_shm_eligible_clones"}
# The column order is the export contract: docs/data_dictionary.csv lists
# the headline rows in this order. Topic tables hold every measurement the
# headline leaves out. The QC and diversity block covers IGH and TRB whatever loci
# params.rarefaction.targets configures.
HEADLINE_COLUMNS = (
    *IDENTITY_COLUMNS, "input_read_pairs", *YIELD_COLUMNS,
    *(f"{locus}_{field}" for locus in ("igh", "trb") for field in DIVERSITY_FIELDS
      if field != "diversity_status"),
    "class_callable_consensus_count",
    *(f"class_{isotype.lower()}_fraction" for isotype in CLASSES),
    "log2_igg_iga",
    "igh_n_shm_eligible_clones_whole_v",
    "igh_shm_eligible_fraction_of_measured_members_whole_v",
    "igh_mean_clone_median_shm_whole_v",
    "igh_fraction_clones_zero_median_shm_whole_v",
    *(f"{group.lower()}_{field}_whole_v" for group in SHM_GROUPS[1:]
      for field in ("n_shm_eligible_clones", "mean_clone_median_shm")),
)
assert len(HEADLINE_COLUMNS) == 41


def _read(path: Path, keys: tuple[str, ...], columns: tuple[str, ...],
          manifest: dict[str, dict[str, str]] | None = None) -> dict[tuple, dict]:
    """Read TSV cells unchanged and check identities before joining a table."""
    path = Path(path)
    result = {}
    with path.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        header = reader.fieldnames or []
        if not header or len(header) != len(set(header)) or any(not name for name in header):
            raise ValueError(f"{path.name}: missing or repeated column names")
        missing = set(keys + columns) - set(header)
        if missing:
            raise ValueError(f"{path.name}: missing required columns {sorted(missing)}")
        for row in reader:
            if None in row or None in row.values():
                raise ValueError(f"{path.name}:{reader.line_num}: malformed row")
            key = tuple(row[column] for column in keys)
            if any(not value or value != value.strip() for value in key):
                raise ValueError(f"{path.name}:{reader.line_num}: empty or padded key")
            if key in result:
                raise ValueError(f"{path.name}: duplicate key {key!r}")
            if manifest is not None:
                identity = manifest.get(row["library_id"])
                if identity is None:
                    raise ValueError(f"{path.name}: out-of-scope library {row['library_id']!r}")
                for column in IDENTITY_COLUMNS:
                    if column in row and row[column] != identity[column]:
                        raise ValueError(f"{path.name}: conflicting {column} for {row['library_id']!r}")
            result[key] = row
    return result


def _required(rows: dict, key: tuple, source: Path) -> dict:
    if key not in rows:
        raise ValueError(f"{Path(source).name}: missing required key {key!r}")
    return rows[key]


def _measurement(row: dict, column: str) -> str | None:
    """Validate numeric cells without changing their source precision.

    An empty cell is the pipeline null. A measurement column has no legitimate
    string value, so a cell spelling the two letters NA is read as missing here
    as well, although table_io.R keeps that spelling distinguishable in the
    record tables.
    """
    value = row[column]
    if value in ("", "NA"):
        return None
    if not math.isfinite(float(value)):
        raise ValueError(f"{row['library_id']!r}: {column} must be finite or missing")
    return value


def _copy(row: dict, fields: tuple[str, ...], prefix: str = "") -> dict:
    copied = {}
    for field in fields:
        if field.endswith("_status"):
            value = row[field]
            if value in ("", "NA"):
                raise ValueError(f"{row['library_id']!r}: {field} requires an availability reason")
        else:
            value = _measurement(row, field)
        copied[prefix + field] = value
    return copied


def build_headline(manifest: Path, yields: Path, qc: Path, diversity: Path,
                   isotype: Path, isotype_summary: Path,
                   shm: Path, *, scope_id: str | None = None) -> list[dict[str, Any]]:
    """Return the fixed scientific sheet in export order.

    Numeric source cells retain their TSV spelling; missing values become None.
    Each source table must carry a row for every manifest library and key; a
    missing row stops the build. An unavailable cell is never read as an empty
    biological repertoire.
    The manifest supplies the identity columns, and every identity column a
    source carries must agree with it.
    """
    manifest_rows = _read(manifest, ("library_id",), IDENTITY_COLUMNS)
    libraries = {key[0]: {column: row[column] for column in IDENTITY_COLUMNS}
                 for key, row in manifest_rows.items()}
    scopes = {row["scope_id"] for row in libraries.values()}
    if len(scopes) != 1 or "" in scopes or (scope_id is not None and scopes != {scope_id}):
        raise ValueError("Headline requires one nonempty scope matching the manifest")
    if any(not value or value != value.strip() for row in libraries.values() for value in row.values()):
        raise ValueError("Manifest identities must be nonempty strings without surrounding whitespace")

    yield_fields = ("input_read_pairs", "yield_status", *YIELD_COLUMNS)
    yield_rows = _read(yields, ("library_id",), yield_fields, libraries)
    qc_rows = _read(qc, ("library_id", "locus"), QC_FIELDS, libraries)
    diversity_rows = _read(diversity, ("library_id", "locus"), DIVERSITY_FIELDS, libraries)
    class_rows = _read(isotype, ("library_id", "isotype_level", "isotype"), ("support_fraction",), libraries)
    composition_fields = ("class_callable_consensus_count", "igg_iga_status", "log2_igg_iga")
    composition_rows = _read(isotype_summary, ("library_id",),
                             (*composition_fields, "common_qc_consensus_count"), libraries)
    shm_rows = _read(shm, ("library_id", "isotype_group", "v_region"), SHM_FIELDS, libraries)

    rows = []
    for library_id, identity in libraries.items():
        signal = _required(yield_rows, (library_id,), yields)
        row = {**identity, **_copy(signal, yield_fields)}
        for locus in ("IGH", "TRB"):
            copied = {**_copy(_required(qc_rows, (library_id, locus), qc), QC_FIELDS),
                      **_copy(_required(diversity_rows, (library_id, locus), diversity), DIVERSITY_FIELDS)}
            row.update({f"{locus.lower()}_{field}": value for field, value in copied.items()})
        for group in CLASSES:
            source = _required(class_rows, (library_id, "class", group), isotype)
            row[f"class_{group.lower()}_fraction"] = _measurement(source, "support_fraction")
        composition = _required(composition_rows, (library_id,), isotype_summary)
        row.update(_copy(composition, composition_fields))
        total = _measurement(composition, "common_qc_consensus_count")
        callable_support = row["class_callable_consensus_count"]
        if total is None or callable_support is None or not 0 <= float(callable_support) <= float(total):
            raise ValueError(f"{library_id!r}: isotype callability requires valid authoritative counts")
        for group in SHM_GROUPS:
            source = _required(shm_rows, (library_id, group, "whole_v"), shm)
            row.update({f"{('igh' if group == 'all_IGH' else group.lower())}_{SHM_NAMES.get(field, field)}_whole_v": value
                        for field, value in _copy(source, SHM_FIELDS).items()})
        for field, numerator, denominator in (
                ("member", "igh_n_members_eligible_whole_v", "igh_n_members_measured_whole_v"),
                ("clone", "igh_n_shm_eligible_clones_whole_v", "igh_richness")):
            eligible, total = row[numerator], row[denominator]
            fraction = None
            if eligible is not None and total is not None:
                if not 0 <= float(eligible) <= float(total):
                    raise ValueError(f"{library_id!r}: invalid SHM {field} eligibility counts")
                if float(total) > 0:
                    fraction = float(eligible) / float(total)
            row["igh_shm_eligible_fraction_of_measured_members_whole_v" if field == "member"
                else "igh_shm_eligible_clone_fraction_whole_v"] = fraction
        rows.append({column: row[column] for column in HEADLINE_COLUMNS})
    return rows


def write_headline(rows: list[dict[str, Any]], output: Path) -> None:
    """Atomically write the scientific sheet, with empty cells for missing values."""
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_suffix(".tmp")
    try:
        with temporary.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=HEADLINE_COLUMNS, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            for row in rows:
                writer.writerow({column: "" if value is None else value
                                 for column, value in row.items()})
        temporary.replace(output)
    finally:
        temporary.unlink(missing_ok=True)
