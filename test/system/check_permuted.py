#!/usr/bin/env python3
"""Check that renumbering TRUST4 CDR3 variants leaves every released value unchanged.

TIE_A and TIE_B hold the same records under different variant indices (see
test/data_builders/permuted_indices.py). Their rows in every released table
and AIRR export must agree once library, donor, tissue and sequence IDs are
mapped across and IGH clones are named by their members, since SCOPer's
numeric labels follow input order.
"""

import argparse
import csv
import gzip
from pathlib import Path
import re
import sys

import yaml

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from test.data_builders.permuted_indices import (
    DONORS, IDENTITY, PERMUTED, SOURCE, TIED, read_tsv, renumber, variant_counts)
from test.support.reference_tables import MISSING, TOLERANCE

POOLED = re.compile(r"^(TIE_[AB])(\W)(assemble\d+)_(\d+)$")
BACKBONE = re.compile(r"^(TIE_[AB])(\W)(assemble\d+)$")
SOURCE_ID = re.compile(r"^(assemble\d+)_(\d+)$")
CLONE = re.compile(r"^cohort\W(DA|DB)\W\d+$")
ID_COLUMNS = ("library_id", "donor_id", "tissue_id")
REPERTOIRE_LINKS = ("repertoire_id", "sample_processing_id", "data_processing_id")


class Normaliser:
    """Map TIE_B values to their TIE_A equivalents and name clones by their members."""

    def __init__(self, counts: dict[str, int]):
        self.counts = counts
        self.exact = {PERMUTED: IDENTITY, DONORS[PERMUTED]: DONORS[IDENTITY],
                      f"{PERMUTED}_tissue": f"{IDENTITY}_tissue"}
        self.clones: dict[str, str] = {}

    def sequence(self, value: str, permuted: bool) -> str:
        if match := POOLED.match(value):
            _, sep, contig, index = match.groups()
            variant = f"{contig}_{index}"
            return f"{IDENTITY}{sep}{renumber(variant, self.counts) if permuted else variant}"
        if match := SOURCE_ID.match(value):
            return renumber(value, self.counts) if permuted else value
        return value

    def learn_clones(self, members: list[dict[str, str]]) -> None:
        grouped: dict[str, list[str]] = {}
        for row in members:
            permuted = row["library_id"] == PERMUTED
            grouped.setdefault(row["clone_id"], []).append(self.sequence(row["sequence_id"], permuted))
        self.clones = {clone: f"clone[{min(ids)}]" for clone, ids in grouped.items()}

    def value(self, value: str, permuted: bool) -> str:
        if value in self.exact:
            return self.exact[value]
        if CLONE.match(value):
            if value not in self.clones:
                raise AssertionError(f"clone {value!r} has no members in clone_members")
            return self.clones[value]
        if BACKBONE.match(value):
            return BACKBONE.sub(rf"{IDENTITY}\2\3", value)
        normalised = self.sequence(value, permuted)
        if normalised == value and value.startswith(f"{PERMUTED}_"):
            return IDENTITY + value[len(PERMUTED):]
        return normalised


def same(left: str, right: str) -> bool:
    if left == right:
        return True
    if left in MISSING or right in MISSING:
        return False
    try:
        a, b = float(left), float(right)
    except ValueError:
        return False
    return abs(a - b) <= TOLERANCE * max(1.0, abs(a), abs(b))


def compare(name: str, header: list[str], first: list[list[str]], second: list[list[str]]) -> None:
    first, second = sorted(first), sorted(second)
    if len(first) != len(second):
        raise AssertionError(f"{name}: {len(first)} {IDENTITY} rows against {len(second)} {PERMUTED} rows")
    for a, b in zip(first, second):
        for column, left, right in zip(header, a, b):
            if not same(left, right):
                raise AssertionError(
                    f"{name}, {column}: {left!r} for {IDENTITY} against {right!r} for {PERMUTED}\n"
                    f"  {IDENTITY}: {a}\n  {PERMUTED}: {b}")


def side(row: dict[str, str]) -> str | None:
    values = {row.get(column) for column in ID_COLUMNS}
    if values & {IDENTITY, DONORS[IDENTITY], f"{IDENTITY}_tissue"}:
        return IDENTITY
    if values & {PERMUTED, DONORS[PERMUTED], f"{PERMUTED}_tissue"}:
        return PERMUTED
    return None


def read_rows(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return list(reader.fieldnames or []), list(reader)


def check_tie(release: Path) -> None:
    """Fail if the constructed variants no longer tie on every ranked criterion."""
    _, rows = read_rows(release / "airr" / f"{IDENTITY}_rearrangements.tsv.gz")
    tied = [row for row in rows if row["source_sequence_id"] in TIED]
    if len(tied) != 2:
        raise AssertionError(f"{IDENTITY} export lacks the constructed variants {TIED}")
    first, second = tied
    for label, measure in (("support", lambda row: row["consensus_count"]),
                           ("V identity", lambda row: row["v_identity"]),
                           ("alignment length", lambda row: len(row["sequence_alignment"]))):
        if measure(first) != measure(second):
            raise AssertionError(f"The constructed variants differ in {label}, so they do not tie")
    if first["junction"] == second["junction"]:
        raise AssertionError("The constructed variants share a junction, so content cannot decide")


def check_case(work: Path) -> None:
    config = yaml.safe_load((work / "config.yaml").read_text())
    release = Path(config["outs"]["results_dir"]) / "scopes" / "cohort" / "release"
    check_tie(release)
    _, cdr3 = read_tsv(Path("test/data/post_trust4/libraries") / SOURCE / "trust4" /
                       f"{SOURCE}_cdr3.out", header=False)
    normaliser = Normaliser(variant_counts(cdr3))
    members, = release.rglob("clone_members.tsv")
    normaliser.learn_clones(read_rows(members)[1])

    compared = 0
    for path in sorted(release.rglob("*.tsv")):
        header, rows = read_rows(path)
        split: dict[str, list[list[str]]] = {IDENTITY: [], PERMUTED: []}
        for row in rows:
            owner = side(row)
            if owner is not None:
                permuted = owner == PERMUTED
                split[owner].append([normaliser.value(row[column], permuted) for column in header])
        if split[IDENTITY] or split[PERMUTED]:
            compare(str(path.relative_to(release)), header, split[IDENTITY], split[PERMUTED])
            compared += 1
    exports = {}
    for library in (IDENTITY, PERMUTED):
        header, rows = read_rows(release / "airr" / f"{library}_rearrangements.tsv.gz")
        # These Repertoire links are derived from the library ID, so they differ by design.
        header = [column for column in header if column not in REPERTOIRE_LINKS]
        exports[library] = [[normaliser.value(row[column], library == PERMUTED) for column in header]
                            for row in rows]
    compare("airr exports", header, exports[IDENTITY], exports[PERMUTED])
    print(f"Renumbered variants gave identical results in {compared} tables and the AIRR exports.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--work-dir", required=True, type=Path)
    args = parser.parse_args()
    check_case(args.work_dir.resolve())
