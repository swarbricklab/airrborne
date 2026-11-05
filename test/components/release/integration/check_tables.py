#!/usr/bin/env python3
"""Validate scientific source tables and their delivered data dictionaries.

Tables must match their dictionary header order, finite numeric types, declared
vocabularies, nullability and primary keys. Declared references must resolve
within each release. Unrecognised files fail unless
workflow.lib.release_bundle.is_exempt() explicitly excludes them.

Each delivered table must retain its source scientific columns and dictionary
contract. check_release.py checks exported observations and links. The record-QC
reference calculations check source retention, support transfer and exact
inference-sequence uniqueness.

Run from the repository root:

    python test/components/release/integration/check_tables.py \\
        --results results/airrborne --dictionary docs/data_dictionary.csv

Value checks read every row. Use a batch job for a full cohort.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import math
import sys
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path

# Read the workflow's filename-to-table mapping so these checks use the same
# scientific table names as publication.
sys.path.insert(0, str(Path(__file__).resolve().parents[4]))
from workflow.lib.release_bundle import (  # noqa: E402
    is_exempt,
    is_header_only,
    table_for,
)
from workflow.lib.scope_qc_audit import COLUMNS as QC_AUDIT_SHARD_COLUMNS  # noqa: E402

# The dictionary spells a null as this token, because an empty cell in the CSV
# means "no vocabulary declared" rather than "the empty string is allowed".
EMPTY_TOKEN = "empty"


@dataclass(frozen=True)
class Reference:
    """Where one column's values must appear, on rows meeting a condition.

    `condition` is a (column, operator, value) triple, with operator = or !=, or
    None when the reference holds on every row.
    """

    column: str
    table: str
    target: str
    condition: tuple[str, str, str] | None = None


@dataclass
class TableContract:
    """What the dictionary says one released table must look like."""

    columns: list[str] = field(default_factory=list)
    # Column name -> the values it may take. Columns with no declared vocabulary
    # are absent.
    options: dict[str, set[str]] = field(default_factory=dict)
    primary_key: list[str] = field(default_factory=list)
    types: dict[str, str] = field(default_factory=dict)
    # Columns the dictionary declares nullable=FALSE, whose cells must not be empty.
    required: set[str] = field(default_factory=set)
    references: list[Reference] = field(default_factory=list)


def open_table(path: Path, mode: str = "rt"):
    return gzip.open(path, mode, newline="") if path.suffix == ".gz" else path.open(mode, newline="")


def table_paths(root: Path) -> list[Path]:
    return sorted([*root.rglob("*.tsv"), *root.rglob("*.tsv.gz")])


def read_condition(text: str) -> tuple[str, str, str] | None:
    """Read one references_when entry, such as group_unit=clone or locus!=IGH."""
    if not text:
        return None
    for operator in ("!=", "="):
        column, found, value = text.partition(operator)
        if found:
            return column, operator, value
    raise ValueError(f"Unreadable references_when entry: {text!r}")


def read_references(row: dict[str, str]) -> list[Reference]:
    """Read a row's references. Alternatives are |-separated, in step across fields."""
    if not row["references_table"]:
        return []
    tables = row["references_table"].split("|")
    targets = row["references_column"].split("|")
    conditions = row["references_when"].split("|") if row["references_when"] else [""] * len(tables)
    if not len(tables) == len(targets) == len(conditions):
        raise ValueError(f"{row['table']}.{row['column']}: reference fields differ in length")
    return [Reference(row["column"], table, target, read_condition(condition))
            for table, target, condition in zip(tables, targets, conditions)]


def load_dictionary(path: Path) -> dict[str, TableContract]:
    """Read the dictionary fields the released tables are checked against."""
    contracts: dict[str, TableContract] = defaultdict(TableContract)
    with path.open(newline="") as handle:
        for row in csv.DictReader(handle):
            contract = contracts[row["table"]]
            column = row["column"]
            if column in contract.columns:
                raise ValueError(f"{path}: repeated {row['table']}.{column}")
            contract.columns.append(column)
            contract.types[column] = row["type"]
            if row["options"]:
                contract.options[column] = {
                    "" if value == EMPTY_TOKEN else value
                    for value in row["options"].split("|")
                }
            if row["is_primary_key"] == "TRUE":
                contract.primary_key.append(column)
            if row["nullable"] not in {"TRUE", "FALSE"}:
                raise ValueError(f"{path}: {row['table']}.{column} has nullable {row['nullable']!r}")
            if row["nullable"] == "FALSE":
                contract.required.add(column)
            contract.references.extend(read_references(row))
    return dict(contracts)


# Columns a scope step adds to a library table: readable V/J labels on
# clonotypes, as LIBRARY_METRIC_LABELS in workflow/rules/scope.smk declares them,
# and the member assignments annotate_scope_airr adds to each AIRR export.
SCOPE_ADDED_COLUMNS = {"clonotypes": {"v_gene_names", "j_gene_names"},
                       "exported_AIRR": {"clone_id"}}


def expected_columns(path: Path, contract: TableContract, table: str) -> list[str]:
    """The columns this particular file owes, which is not always the whole table.

    Scope collectors append specimen and donor context to library measurements,
    readable V/J labels to clonotypes and clone_id to the AIRR exports. Clone abundance also receives scope
    identity when its two arms are collected. Library and donor QC audit shards
    keep the column order qc_audit.R writes, and scope_qc_audit.py leads the
    collected qc_audit.tsv with the identifiers.
    """
    if "_qc_audit_" in path.name:
        if set(QC_AUDIT_SHARD_COLUMNS) != set(contract.columns):
            raise ValueError("QC audit shard columns differ from the qc_audit dictionary")
        return list(QC_AUDIT_SHARD_COLUMNS)
    if "libraries" in path.parts:
        added = {"scope_id", "tissue_id", "donor_id"} | SCOPE_ADDED_COLUMNS.get(table, set())
        return [column for column in contract.columns if column not in added]
    return contract.columns


def check_header(path: Path, expected: list[str], table: str) -> list[str]:
    with open_table(path) as handle:
        observed = next(csv.reader(handle, delimiter="\t"), [])
    if observed != expected:
        return [f"{path}: header differs from {table}\n"
                f"  expected={expected}\n  observed={observed}"]
    return []


def check_values(path: Path, contract: TableContract,
                 expected: list[str]) -> list[str]:
    """Check vocabularies, nullability and primary-key uniqueness over every row.

    Only columns this file owes are checked, so a table whose per-library form is
    narrower than its aggregate is not asked for columns it never carried.
    """
    vocabularies = {column: values for column, values in contract.options.items()
                    if column in expected}
    key = [column for column in contract.primary_key if column in expected]
    required = [column for column in expected if column in contract.required]
    empty_required: dict[str, int] = defaultdict(int)
    offending: dict[str, set[str]] = defaultdict(set)
    wrong_types: dict[str, set[str]] = defaultdict(set)
    malformed = 0
    seen_keys: set[tuple[str, ...]] = set()
    duplicate_keys: set[tuple[str, ...]] = set()
    with open_table(path) as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if None in row or None in row.values():
                malformed += 1
                continue
            for column in expected:
                value, kind = row[column], contract.types[column]
                if not value or kind not in {"numeric", "integer"}:
                    continue
                try:
                    number = float(value)
                    valid = math.isfinite(number) and (kind != "integer" or number.is_integer())
                except ValueError:
                    valid = False
                if not valid:
                    wrong_types[column].add(value)
            for column, values in vocabularies.items():
                value = row[column]
                if value not in values:
                    offending[column].add(value)
            for column in required:
                if not row[column]:
                    empty_required[column] += 1
            if key:
                key_values = tuple(row[column] for column in key)
                if key_values in seen_keys:
                    duplicate_keys.add(key_values)
                seen_keys.add(key_values)

    failures = []
    if malformed:
        failures.append(f"{path}: {malformed} malformed row(s)")
    for column, values in sorted(wrong_types.items()):
        failures.append(f"{path}: {column} contains invalid {contract.types[column]} values: {sorted(values)[:5]}")
    for column, count in sorted(empty_required.items()):
        failures.append(f"{path}: {column} is declared not nullable but has {count} empty cell(s)")
    for column in sorted(offending):
        examples = sorted(offending[column], key=str)[:5]
        failures.append(
            f"{path}: {column} holds {len(offending[column])} value(s) outside "
            f"its declared vocabulary\n"
            f"  declared={sorted(vocabularies[column])}\n  found={examples}"
        )
    if duplicate_keys:
        examples = sorted(duplicate_keys)[:5]
        failures.append(
            f"{path}: primary key {key} repeats for "
            f"{len(duplicate_keys)} value(s)\n  examples={examples}"
        )
    return failures


def validate_tables(results: Path, dictionary: dict[str, TableContract]) -> int:
    checked = 0
    failures = []
    unrecognised = []
    for path in table_paths(results):
        if is_exempt(path):
            continue
        # Check source table values once here. The publisher normally hardlinks
        # those tables and can copy them when linking is unavailable. Bundle
        # headers and the shipped dictionary are checked by
        # validate_release_bundles(); check_release.py checks the
        # actual handover and can compare it with source results.
        if "release" in path.parts or "intermediate" in path.parts:
            continue
        table = table_for(path)
        if table is None:
            unrecognised.append(path)
            continue
        if table not in dictionary:
            failures.append(f"{path}: table {table} is absent from data dictionary")
            continue
        contract = dictionary[table]
        expected = expected_columns(path, contract, table)
        header_failures = check_header(path, expected, table)
        failures.extend(header_failures)
        # The value checks read the file by column name, so they mean nothing
        # until the header is the one the dictionary declares.
        if not header_failures and not is_header_only(path):
            failures.extend(check_values(path, contract, expected))
        checked += 1
    if unrecognised:
        listed = "\n".join(f"  {path}" for path in unrecognised[:20])
        more = ("\n  ..." if len(unrecognised) > 20 else "")
        failures.append(
            f"{len(unrecognised)} released TSV(s) match no dictionary table and "
            f"are not listed as exempt:\n{listed}{more}\n"
            f"Give each one a table in the data dictionary and a rule in "
            f"table_for(), or a reason in EXEMPT_SUFFIXES."
        )
    if failures:
        raise SystemExit("\n".join(failures))
    return checked


def check_references(bundle: Path, contracts: dict[str, TableContract]) -> list[str]:
    """Check that each declared reference resolves within one release.

    A target outside the release's tables, such as Repertoire.json, is not
    checked here. Target values are held in memory, so a full cohort needs a
    batch job.
    """
    paths: dict[str, list[Path]] = defaultdict(list)
    for path in table_paths(bundle):
        table = table_for(path)
        if table is not None:
            paths[table].append(path)
    known: dict[tuple[str, str], set[str]] = {}

    def target_values(table: str, column: str) -> set[str]:
        if (table, column) not in known:
            values: set[str] = set()
            for path in paths[table]:
                with open_table(path) as handle:
                    values.update(row[column] for row in csv.DictReader(handle, delimiter="\t"))
            known[(table, column)] = values
        return known[(table, column)]

    def applies(row: dict[str, str], reference: Reference) -> bool:
        if reference.condition is None:
            return True
        column, operator, value = reference.condition
        return (row[column] == value) == (operator == "=")

    failures = []
    for table, files in sorted(paths.items()):
        contract = contracts.get(table, TableContract())
        references = [reference for reference in contract.references
                      if reference.table in paths]
        if not references:
            continue
        targets = [target_values(reference.table, reference.target) for reference in references]
        missing: dict[Reference, set[str]] = defaultdict(set)
        for path in files:
            with open_table(path) as handle:
                for row in csv.DictReader(handle, delimiter="\t"):
                    for reference, values in zip(references, targets):
                        value = row[reference.column]
                        if value and value not in values and applies(row, reference):
                            missing[reference].add(value)
        for reference, values in missing.items():
            failures.append(f"{bundle}: {table}.{reference.column} holds {len(values)} value(s) absent "
                            f"from {reference.table}.{reference.target}: {sorted(values)[:5]}")
    return failures


def validate_release_bundles(results: Path,
                             dictionary: dict[str, TableContract]) -> int:
    """Check delivered tables' declared scientific columns and their references."""
    failures = []
    bundles = sorted(results.glob("scopes/*/release"))
    for bundle in bundles:
        shipped = bundle / "data_dictionary.csv"
        shipped_contracts = load_dictionary(shipped)
        for path in table_paths(bundle):
            table = table_for(path)
            if table is None:
                failures.append(f"{path}: unrecognised shipped table")
                continue
            if table not in shipped_contracts:
                failures.append(f"{path}: table {table} is absent from the shipped dictionary")
                continue
            contract = shipped_contracts[table]
            errors = check_header(path, contract.columns, table)
            failures.extend(errors)
            # Every released table retains the scientific source contract.
            if table not in dictionary or contract != dictionary[table]:
                failures.append(f"{path}: shipped dictionary differs from its source contract")
        failures.extend(check_references(bundle, shipped_contracts))
    if failures:
        raise SystemExit("\n".join(failures))
    return len(bundles)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results", required=True, type=Path)
    parser.add_argument("--dictionary", required=True, type=Path)
    args = parser.parse_args()
    dictionary = load_dictionary(args.dictionary)
    tables = validate_tables(args.results, dictionary)
    bundles = validate_release_bundles(args.results, dictionary)
    print(f"Validated {tables} scientific source tables against {len(dictionary)} dictionary tables")
    print(f"Validated scientific columns and dictionaries in {bundles} scope releases")


if __name__ == "__main__":
    main()
