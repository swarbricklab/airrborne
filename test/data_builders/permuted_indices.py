#!/usr/bin/env python3
"""Construct two libraries that differ only in TRUST4 CDR3 variant numbering.

TRUST4 numbers the CDR3 variants of each contig with a trailing index that can
change between runs. TIE_A copies one public library and gives two variants of
its first IGH contig equal read support, so a ranking that falls back to the
sequence ID would choose differently. TIE_B holds the same records with the
indices of every contig reversed and the rows in reverse order. Each library is
its own donor, so donor-level results can be compared too. Neither library has
assembly producer history, because the edited files no longer match it.
"""

import argparse
import csv
from pathlib import Path
import sys

import yaml

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from test.system.stage_inputs import prepare

SOURCE = "SYN_D1_T2"
IDENTITY = "TIE_A"
PERMUTED = "TIE_B"
DONORS = {IDENTITY: "DA", PERMUTED: "DB"}
# Two variants of one IGH contig with junctions of equal length and sequences of
# equal length. TIE_A raises the second to the first one's support.
TIED = ("assemble0_0", "assemble0_1")


def variant_counts(cdr3: list[list[str]]) -> dict[str, int]:
    """Count the CDR3 variants of each contig from cdr3.out, which lists every one.

    The AIRR table omits partial CDR3 calls, so it can hold a subset of them.
    """
    indices: dict[str, list[int]] = {}
    for row in cdr3:
        indices.setdefault(row[0], []).append(int(row[1]))
    for contig, seen in indices.items():
        if sorted(seen) != list(range(len(seen))):
            raise ValueError(f"cdr3.out indices of {contig} are not 0 to n - 1: {seen}")
    return {contig: len(seen) for contig, seen in indices.items()}


def renumber(sequence_id: str, counts: dict[str, int]) -> str:
    """Reverse the variant index within its contig. Applied twice, it is the identity."""
    contig, index = sequence_id.rsplit("_", 1)
    return f"{contig}_{counts[contig] - 1 - int(index)}"


def read_tsv(path: Path, header: bool) -> tuple[list[str] | None, list[list[str]]]:
    with path.open(newline="") as handle:
        rows = list(csv.reader(handle, delimiter="\t"))
    return (rows[0], rows[1:]) if header else (None, rows)


def write_tsv(path: Path, header: list[str] | None, rows: list[list[str]]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        if header is not None:
            writer.writerow(header)
        writer.writerows(rows)


def write_library(root: Path, library: str, header: list[str], airr: list[list[str]],
                  cdr3: list[list[str]], final: str, depth: str) -> None:
    trust4 = root / "libraries" / library / "trust4"
    metrics = root / "libraries" / library / "metrics"
    trust4.mkdir(parents=True)
    metrics.mkdir()
    write_tsv(trust4 / f"{library}_airr.tsv", header, airr)
    write_tsv(trust4 / f"{library}_cdr3.out", None, cdr3)
    (trust4 / f"{library}_final.out").write_text(final)
    (metrics / f"{library}_seq_depth.tsv").write_text(depth.replace(SOURCE, library))


def prepare_case(work: Path, resources: Path | None = None) -> None:
    config_path = work / "config.yaml"
    config = prepare(Path("test/system/config.yaml"), work, config_path, resources=resources)
    root = Path(config["outs"]["results_dir"])
    config.pop("scope_metadata", None)
    config.pop("membership_csv", None)
    samples = work / "samplesheet.csv"
    samples.write_text("library_id,tissue_id,donor_id,path\n" + "".join(
        f"{library},{library}_tissue,{donor},{library}\n" for library, donor in DONORS.items()))
    config["deps"]["samplesheet"] = str(samples.resolve())
    config_path.write_text(yaml.safe_dump(config, sort_keys=False))

    source = Path("test/data/post_trust4/libraries") / SOURCE
    header, airr = read_tsv(source / "trust4" / f"{SOURCE}_airr.tsv", header=True)
    _, cdr3 = read_tsv(source / "trust4" / f"{SOURCE}_cdr3.out", header=False)
    final = (source / "trust4" / f"{SOURCE}_final.out").read_text()
    depth = (source / "metrics" / f"{SOURCE}_seq_depth.tsv").read_text()

    sequence_id = header.index("sequence_id")
    support = header.index("consensus_count")
    by_id = {row[sequence_id]: row for row in airr}
    if not all(tied in by_id for tied in TIED):
        raise ValueError(f"{SOURCE} no longer carries the variants {TIED}")
    by_id[TIED[1]][support] = by_id[TIED[0]][support]
    write_library(root, IDENTITY, header, airr, cdr3, final, depth)

    counts = variant_counts(cdr3)
    permuted_airr = [row[:sequence_id] + [renumber(row[sequence_id], counts)] + row[sequence_id + 1:]
                     for row in reversed(airr)]
    # cdr3.out keys each variant by contig (column 1) and index (column 2).
    permuted_cdr3 = [[row[0], renumber(f"{row[0]}_{row[1]}", counts).rsplit("_", 1)[1], *row[2:]]
                     for row in reversed(cdr3)]
    write_library(root, PERMUTED, header, permuted_airr, permuted_cdr3, final, depth)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--work-dir", required=True, type=Path)
    parser.add_argument("--resources", type=Path,
                        help="Prepared resources directory to use instead of the configuration's.")
    args = parser.parse_args()
    prepare_case(args.work_dir.resolve(), args.resources)
