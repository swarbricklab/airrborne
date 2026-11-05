#!/usr/bin/env python3
"""Validate public construction evidence, paired reads and native assembly links.

These checks use the standard library and explicit sequence/count relationships.
They do not call pipeline QC, grouping or metric functions. Frozen tool outputs
are observed evidence; construction facts remain separate from recovered support.
"""

import argparse
import csv
import gzip
import hashlib
import json
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1] / "data"


def table(path: Path) -> list[dict]:
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def reverse_complement(sequence: str) -> str:
    return sequence.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def fastq(path: Path) -> list[tuple[str, str, str]]:
    with gzip.open(path, "rt") as handle:
        lines = handle.read().splitlines()
    assert len(lines) % 4 == 0, f"Incomplete FASTQ record: {path}"
    result = []
    for index in range(0, len(lines), 4):
        name, bases, separator, quality = lines[index:index + 4]
        assert name.startswith("@") and separator == "+", f"FASTQ record format: {path}"
        assert len(bases) == len(quality) == 150, f"Read length: {path}"
        assert set(bases) <= set("ACGT") and set(quality) == {"I"}, f"Read bases/quality: {path}"
        result.append((name[1:], bases, quality))
    return result


def check_pair(first: tuple, second: tuple, fragment: dict, transcript: str) -> None:
    assert first[0] == fragment["fragment_id"] + "/1", "R1 name differs from fragment ledger"
    assert second[0] == fragment["fragment_id"] + "/2", "R2 name differs from fragment ledger"
    start, end = int(fragment["start"]), int(fragment["end"])
    assert end - start == 250 and 0 <= start < end <= len(transcript), "Fragment coordinates"
    assert first[1] == transcript[start:start + 150], "R1 sequence differs from fragment"
    assert second[1] == reverse_complement(transcript[end - 150:end]), "R2 sequence differs from fragment"


def backbones(lines: list[str]) -> dict[str, str]:
    assert len(lines) % 6 == 0, "Native backbone requires six lines"
    result = {}
    for index in range(0, len(lines), 6):
        assert lines[index].startswith(">"), "Native backbone header"
        name = lines[index][1:].split()[0]
        sequence = lines[index + 1]
        assert name not in result and set(sequence) <= set("ACGTN"), "Native backbone identity/bases"
        for line in lines[index + 2:index + 6]:
            votes = [int(value) for value in line.split()]
            assert len(votes) == len(sequence), "Vote array length differs from backbone"
            assert min(votes) >= 0, "Negative native base votes"
        result[name] = sequence
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=ROOT)
    args = parser.parse_args()
    root = args.root
    construction = json.loads((root / "construction/construction.json").read_text())
    transcripts = json.loads((root / "construction/transcripts.json").read_text())
    reference = json.loads((root / "construction/reference_alleles.json").read_text())
    sequences = {row["transcript_id"]: row["sequence"] for row in transcripts}
    assert len(sequences) == len(transcripts) == 11
    alleles = {row["name"]: row["sequence"] for entries in reference.values() for row in entries}
    trees = [{row["node_id"]: row for row in table(root / "construction" / f"family_{index}.tsv")} for index in (0, 1)]
    for row in transcripts:
        pieces = []
        for segment in "vdj":
            sequence = alleles[row[f"{segment}_call"]]
            begin = row["trims"][f"{segment}_trim_5"]
            end = len(sequence) - row["trims"][f"{segment}_trim_3"]
            pieces.append(sequence[begin:end])
            if segment != "j":
                pieces.append(row["np1" if segment == "v" else "np2"])
        assert "".join(pieces) == row["ancestor"], f"Ancestor reconstruction: {row['transcript_id']}"
        edits = [dict(position=i, before=a, after=b) for i, (a, b) in
                 enumerate(zip(row["ancestor"], row["variable"])) if a != b]
        assert edits == row["edits"], f"Mutation coordinates: {row['transcript_id']}"
        if row["lineage_node"]:
            clone, node = row["lineage_node"].split(":")
            assert trees[int(clone)][node]["sequence"].upper() == row["variable"], "Lineage node sequence"
            assert trees[int(clone)]["0"]["sequence"].upper() == row["ancestor"], "Lineage ancestor"
        assert row["variable"][row["junction_start"]:row["junction_end"]] == row["junction"], "Junction coordinates"
        assert row["variable"][-1] == row["constant"][0] == "G", "J-to-C shared nucleotide"
        assert row["c_start"] == len(row["variable"]) - 1, "Constant start coordinate"
        assert row["sequence"] == row["variable"] + row["constant"][1:], "Transcript join"
        assert len(row["constant"]) == row["constant_reference_end"], "Constant segment length"
    with gzip.open(root / "construction/fragments.tsv.gz", "rt") as handle:
        fragments = list(csv.DictReader(handle, delimiter="\t"))
    assert len({row["fragment_id"] for row in fragments}) == len(fragments), "Repeated fragment identifier"
    sample_rows = list(csv.DictReader((root / "samplesheet.csv").open()))
    libraries = {row["library_id"] for row in sample_rows}
    assert libraries == set(construction["representation"]), "Construction/library membership"
    assert libraries == {row["library_id"] for row in fragments}, "Fragment/library membership"
    counts = Counter(row["library_id"] for row in fragments)
    native_count = 0
    recovery = table(root / "construction/native_recovery.tsv")
    for library in sorted(libraries):
        for lane in (1, 2):
            expected = [row for row in fragments if row["library_id"] == library and int(row["lane"]) == lane]
            reads = [fastq(root / "reads" / library / f"{library}_L00{lane}_R{mate}_001.fastq.gz") for mate in (1, 2)]
            assert len(reads[0]) == len(reads[1]) == len(expected), "Lane fragment count"
            for first, second, fragment in zip(*reads, expected):
                check_pair(first, second, fragment, sequences[fragment["transcript_id"]])
        target = root / "post_trust4" / "libraries" / library
        depth = table(target / "metrics" / f"{library}_seq_depth.tsv")
        assert len(depth) == 1 and depth[0]["gussid"] == library, "Depth library identity"
        for field, value in [("n_fastq_pairs", 2), ("total_read_pairs", counts[library]),
                             ("total_reads", counts[library] * 2), ("total_bases", counts[library] * 300)]:
            assert int(depth[0][field]) == value, f"Sequencing depth {library}/{field}"
        assembly = target / "trust4"
        lines = (assembly / f"{library}_final.out").read_text().splitlines()
        native_backbones = backbones(lines)
        cdr3 = [line.split("\t") for line in (assembly / f"{library}_cdr3.out").read_text().splitlines()]
        cdr3_ids = {f"{row[0]}_{row[1]}" for row in cdr3}
        native = table(assembly / f"{library}_airr.tsv")
        assert len({row["sequence_id"] for row in native}) == len(native), "Repeated native record ID"
        assert {row["locus"] for row in native} == {"IGH", "TRB"}, "Native IGH/TRB recovery"
        expected_junctions = {row["junction"] for row in transcripts
                              if row["transcript_id"] in construction["representation"][library]}
        assert expected_junctions == {row["junction"] for row in native}, "Constructed/native junction sets differ"
        observed = {(row["sequence_id"], row["locus"], int(row["consensus_count"]),
                     row["productive"], row["c_call"], len(row["sequence"])) for row in native}
        recorded = {(row["sequence_id"], row["locus"], int(row["native_support"]),
                     row["native_productive"], row["native_c_call"], int(row["assembled_bases"]))
                    for row in recovery if row["library_id"] == library}
        for row in native:
            assert row["sequence_id"] in cdr3_ids, "Native CDR3 link"
            assert row["sequence_id"].rsplit("_", 1)[0] in native_backbones, "Native backbone link"
            assert int(row["consensus_count"]) > 0, "Native recovered support"
        assert observed == recorded, "Native recovery ledger differs"
        producer = json.loads((assembly / f"{library}_assembly_provenance.json").read_text())
        assert producer["library_id"] == library and producer["stage"] == "assembly", "Producer identity"
        assert producer["inputs"]["reference"]["sha256"] == construction["reference_sha256"], "Producer reference"
        for entry in producer["outputs"].values():
            path = assembly / entry["name"]
            assert hashlib.sha256(path.read_bytes()).hexdigest() == entry["sha256"], "Producer output checksum"
        native_count += len(native)
    manifest = json.loads((root / "checksums.json").read_text())
    data_files = {str(path.relative_to(root)) for path in root.rglob("*")
                  if path.is_file() and path.suffix != ".md" and path.name != "checksums.json"}
    assert set(manifest) == data_files, "Checksum catalogue differs from data files"
    reviewed = table(root / "expected/expected_records.tsv")
    assert {(row["library_id"], row["sequence_id"]) for row in reviewed} == {
        (row["library_id"], row["sequence_id"]) for row in recovery}, "Reviewed/native record IDs differ"
    for name, expected in manifest.items():
        assert hashlib.sha256((root / name).read_bytes()).hexdigest() == expected, f"Frozen artifact changed: {name}"
    print(f"Validated {len(transcripts)} transcripts, {len(fragments)} read pairs and {native_count} native records across {len(libraries)} libraries.")


if __name__ == "__main__":
    main()
