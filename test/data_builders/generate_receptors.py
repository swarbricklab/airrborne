#!/usr/bin/env python3
"""Build the public receptor panel and paired reads in a new output directory.

GenAIRR supplies rearrangements and related IGH sequences. The declared library
representation and read tiling are controlled examples for pipeline testing.
Coordinates are zero-based and ends are exclusive. Run with GenAIRR 2.3.1 and
Biopython 1.85. Normal acceptance uses the frozen files and needs neither package.
"""

import argparse
import csv
import gzip
import hashlib
import io
import json
from pathlib import Path

import GenAIRR as ga
from Bio import SeqIO
from Bio.Seq import Seq

SEED = 1337
REPRESENTATION = {
    "SYN_D1_T1": {"IGH_root": 14, "IGH_early": 10, "IGH_other": 8, "TRB_1": 4, "TRB_2": 2},
    "SYN_D1_T2": {"IGH_early": 12, "IGH_late": 10, "IGH_distinct": 8, "IGH_short_c": 3, "TRB_1": 2, "TRB_3": 4},
    "SYN_D2_T1": {"IGH_root": 16, "IGH_other_late": 12, "TRB_1": 3, "TRB_2": 5},
    "SYN_D2_T2": {"IGH_other": 14, "IGH_distinct": 12, "IGH_stop": 3, "TRB_2": 4, "TRB_3": 2},
}


def write_json(path: Path, value) -> None:
    path.write_text(json.dumps(value, indent=2) + "\n")


def matched_reference(locus: str, reference: dict[str, str]):
    """Retain exact shared alleles and GenAIRR's curated anchor coordinates."""
    by_sequence = {}
    for name, gapped in reference.items():
        if name.startswith(locus):
            by_sequence.setdefault(gapped.replace(".", ""), []).append(name)
    alias = "human_igh" if locus == "IGH" else "human_tcrb"
    original = ga.Experiment.on(alias).curate_refdata("functional_anchors_only").recombine().compile().refdata
    ref = ga.RefDataConfig.vdj()
    catalogue = []
    for segment in "vdj":
        selected = set()
        for index in range(getattr(original, f"{segment}_pool_size")()):
            allele = getattr(original, f"{segment}_allele")(index)
            sequence = allele.seq().decode().upper().replace(".", "")
            names = by_sequence.get(sequence, [])
            if not names:
                continue
            name = sorted(names)[0]
            gene = name.split("*")[0]
            if gene in selected:
                continue
            anchor = allele.anchor
            if segment == "v":
                assert str(Seq(sequence[anchor:anchor + 3]).translate()) == "C"
                assert reference[name][:309].replace(".", "") == sequence[:anchor]
            if segment == "j":
                assert str(Seq(sequence[anchor:anchor + 3]).translate()) in {"F", "W"}
            kwargs = {"anchor": anchor} if segment != "d" else {}
            getattr(ref, f"add_{segment}_allele")(name, gene, sequence.encode(), **kwargs)
            catalogue.append(dict(segment=segment, name=name, source_alias=allele.name,
                                  sequence=sequence, anchor=anchor))
            selected.add(gene)
            if len(selected) == 2:
                break
        assert len(selected) == 2
    ref.validate_strict()
    return ref, catalogue


def experiment(ref):
    lengths = [(i, 1.0) for i in range(7)]
    trims = [(i, 1.0) for i in range(4)]
    return (ga.Experiment.on(ref).recombine(np1_lengths=lengths, np2_lengths=lengths)
            .trim(v_3=trims, d_5=trims, d_3=trims, j_5=trims).productive_only())


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    evidence = args.output / "construction"
    evidence.mkdir()
    reference = {r.id: str(r.seq).upper() for r in SeqIO.parse(args.reference, "fasta")}
    refs, catalogues, baselines = {}, {}, {}
    for locus in ("IGH", "TRB"):
        refs[locus], catalogues[locus] = matched_reference(locus, reference)
        baselines[locus] = experiment(refs[locus]).run_records(
            n=200, seed=SEED, strict=True, expose_provenance=True, validate_records=True)
    family = experiment(refs["IGH"]).clonal_lineage(
        n_clones=2, max_generations=4, n_max=100, n_sample=12, rate=0.01,
        lambda_base=2.0, selection_strength=0.0).run_records(
            seed=SEED, expose_provenance=True, validate_records=True)
    assert family.validate_families()
    trees = []
    for index, tree in enumerate(family.lineage_trees):
        tree.validate()
        text = tree.to_node_table_tsv()
        (evidence / f"family_{index}.tsv").write_text(text)
        trees.append({int(r["node_id"]): r for r in csv.DictReader(io.StringIO(text), delimiter="\t")})
    truth = []

    def transcript(name, record, c_call, variable=None, ancestor=None, node=None, c_length=270):
        variable = (variable or record["sequence"]).upper()
        ancestor = (ancestor or variable).upper()
        constant = reference[c_call].replace(".", "")[:c_length]
        assert len(variable) % 3 == 1 and variable[-1] == constant[0] == "G"
        sequence = variable + constant[1:]
        amino_acids = str(Seq(sequence).translate())
        if name != "IGH_stop":
            assert "*" not in amino_acids
        else:
            assert "*" in amino_acids
        junction_start, junction_end = record["junction_start"], record["junction_end"]
        edits = [dict(position=i, before=a, after=b) for i, (a, b) in
                 enumerate(zip(ancestor, variable)) if a != b]
        truth.append(dict(transcript_id=name, locus=record["locus"],
                          generator_sequence_id=record["sequence_id"], lineage_node=node,
                          v_call=record["truth_v_call"], d_call=record["truth_d_call"],
                          j_call=record["truth_j_call"], c_call=c_call,
                          ancestor=ancestor, variable=variable, edits=edits,
                          junction_start=junction_start, junction_end=junction_end,
                          junction=variable[junction_start:junction_end],
                          constant=constant, constant_reference_start=0,
                          constant_reference_end=c_length, join_overlap=1,
                          c_start=len(variable) - 1, sequence=sequence,
                          amino_acids=amino_acids,
                          trims={f"{s}_trim_{end}": record[f"{s}_trim_{end}"]
                                 for s in "vdj" for end in (5, 3)},
                          np1=record["np1"].upper(), np2=record["np2"].upper()))

    for name, clone, node, c_call in [
        ("IGH_root", 0, 0, "IGHM*01"), ("IGH_early", 0, 5, "IGHG1*01"),
        ("IGH_late", 0, 18, "IGHA1*01"), ("IGH_other", 1, 23, "IGHG1*01"),
        ("IGH_other_late", 1, 25, "IGHG1*01"), ("IGH_stop", 0, 17, "IGHM*01"),
    ]:
        root = trees[clone][0]["sequence"].upper()
        record = dict(next(r for r in family.records if r["clone_id"] == clone))
        record["sequence"] = root
        record["sequence_id"] = f"clone{clone}_node{node}"
        record["np1"] = root[record["v_sequence_end"]:record["d_sequence_start"]]
        record["np2"] = root[record["d_sequence_end"]:record["j_sequence_start"]]
        transcript(name, record, c_call, trees[clone][node]["sequence"], root, f"{clone}:{node}")
    transcript("IGH_distinct", baselines["IGH"].records[1], "IGHG1*01")
    transcript("IGH_short_c", baselines["IGH"].records[2], "IGHM*01", c_length=30)
    selected = []
    for record in baselines["TRB"].records:
        if record["truth_d_call"] != "TRBD1*01" or not record["truth_j_call"].startswith("TRBJ1-"):
            continue
        if len(selected) < 2 and record["truth_v_call"] in {r["truth_v_call"] for r in selected}:
            continue
        selected.append(record)
        if len(selected) == 3:
            break
    for index, record in enumerate(selected, 1):
        transcript(f"TRB_{index}", record, "TRBC1*02")
    write_json(evidence / "transcripts.json", truth)
    write_json(evidence / "reference_alleles.json", catalogues)
    write_json(evidence / "construction.json", dict(
        seed=SEED, genairr_version=ga.__version__,
        reference_sha256=hashlib.sha256(args.reference.read_bytes()).hexdigest(),
        read_length=150, fragment_length=250, stride=5, phred_quality=40,
        representation=REPRESENTATION))
    sequences = {r["transcript_id"]: r["sequence"] for r in truth}
    fragments = []
    for library, panel in REPRESENTATION.items():
        directory = args.output / "reads" / library
        directory.mkdir(parents=True)
        streams = {(lane, mate): io.StringIO() for lane in (1, 2) for mate in (1, 2)}
        number = 0
        for name, replicates in panel.items():
            sequence = sequences[name]
            starts = list(range(0, len(sequence) - 250 + 1, 5))
            if starts[-1] != len(sequence) - 250:
                starts.append(len(sequence) - 250)
            for start in starts:
                for replicate in range(replicates):
                    number += 1
                    fragment_id = f"{library}_{number:06d}"
                    lane = 1 + (number % 2)
                    for mate, bases in [(1, sequence[start:start + 150]),
                                        (2, str(Seq(sequence[start + 100:start + 250]).reverse_complement()))]:
                        streams[lane, mate].write(f"@{fragment_id}/{mate}\n{bases}\n+\n{'I' * 150}\n")
                    fragments.append(dict(library_id=library, lane=lane, fragment_id=fragment_id,
                                          transcript_id=name, start=start, end=start + 250,
                                          replicate=replicate))
        for (lane, mate), stream in streams.items():
            path = directory / f"{library}_L00{lane}_R{mate}_001.fastq.gz"
            with path.open("wb") as handle:
                with gzip.GzipFile(fileobj=handle, mode="wb", filename="", mtime=0) as zipped:
                    zipped.write(stream.getvalue().encode())
    with (evidence / "fragments.tsv.gz").open("wb") as handle:
        with gzip.GzipFile(fileobj=handle, mode="wb", filename="", mtime=0) as zipped:
            buffer = io.StringIO()
            writer = csv.DictWriter(buffer, fieldnames=fragments[0], delimiter="\t")
            writer.writeheader()
            writer.writerows(fragments)
            zipped.write(buffer.getvalue().encode())
    print(f"Built {len(truth)} transcripts and {len(fragments)} read pairs in {args.output}")


if __name__ == "__main__":
    main()
