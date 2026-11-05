#!/usr/bin/env python3
"""Construct an empty library and a library containing only a failing receptor.

Two public libraries provide unchanged controls. Constructed inputs have positive
sequencing depth and deliberately lack assembly producer history.
"""

import argparse
import csv
import gzip
from pathlib import Path
import re
import shutil
import sys

import yaml

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from test.system.stage_inputs import prepare

CONTROLS = ("SYN_D1_T1", "SYN_D1_T2")
CONSTRUCTED = ("EMPTY", "FAIL")


def prepare_case(work: Path, resources: Path | None = None) -> None:
    config_path = work / "config.yaml"
    config = prepare(Path("test/system/config.yaml"), work, config_path, resources=resources)
    root = Path(config["outs"]["results_dir"])
    config.pop("scope_metadata", None)
    samples = work / "samplesheet.csv"
    samples.write_text("library_id,tissue_id,donor_id,path\n" + "".join(
        f"{library},{library}_tissue,D1,{library}\n" for library in (*CONTROLS, *CONSTRUCTED)))
    membership = work / "membership.csv"
    membership.write_text("project,library_id,usage\nempty_only,EMPTY,analysed\nempty_only,FAIL,analysed\n")
    config["deps"]["samplesheet"] = str(samples.resolve())
    config["membership_csv"] = str(membership.resolve())
    # Light-chain targets reach empty loci here; the R execution test supplies
    # nonempty IGK/IGL vectors. Positive minima reach the actual workflow here.
    config["params"]["rarefaction"]["targets"].update(IGK=4, IGL=8)
    config["params"]["diversity"]["min_consensus_count"].update(IGH=2000, TRB=400)
    config["params"]["overlap"]["min_consensus_count"].update(IGH=1900, TRB=300)
    config_path.write_text(yaml.safe_dump(config, sort_keys=False))
    frozen = Path("test/data/post_trust4/libraries")
    for library in CONTROLS:
        shutil.copytree(frozen / library, root / "libraries" / library)
    source = frozen / "SYN_D2_T2" / "trust4"
    airr = (source / "SYN_D2_T2_airr.tsv").read_text().splitlines(keepends=True)
    scores = (source / "SYN_D2_T2_cdr3.out").read_text().splitlines(keepends=True)
    votes = re.split(r"(?=^>)", (source / "SYN_D2_T2_final.out").read_text(), flags=re.M)
    for library in CONSTRUCTED:
        directory = root / "libraries" / library
        (directory / "trust4").mkdir(parents=True)
        (directory / "metrics").mkdir()
        prefix = directory / "trust4" / library
        Path(str(prefix) + "_airr.tsv").write_text(airr[0] + ("".join(
            line for line in airr[1:] if line.startswith("assemble3_0\t")) if library == "FAIL" else ""))
        Path(str(prefix) + "_cdr3.out").write_text("".join(
            line for line in scores if line.startswith("assemble3\t")) if library == "FAIL" else "")
        Path(str(prefix) + "_final.out").write_text("".join(
            block for block in votes if block.startswith(">assemble3 ")) if library == "FAIL" else "")
        (directory / "metrics" / f"{library}_seq_depth.tsv").write_text(
            "gussid\tn_fastq_pairs\ttotal_read_pairs\ttotal_reads\ttotal_bases\n"
            f"{library}\t1\t1000\t2000\t300000\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--work-dir", required=True, type=Path)
    parser.add_argument("--resources", type=Path,
                        help="Prepared resources directory to use instead of the configuration's.")
    args = parser.parse_args()
    prepare_case(args.work_dir.resolve(), args.resources)
