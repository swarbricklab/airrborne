#!/usr/bin/env python3
"""Assemble the frozen paired reads and retain native files with producer history.

This is a curation command. Its output directory must be new. A normal full test
runs the workflow's assembly rule and leaves these frozen inputs unchanged.
"""

import argparse
import csv
import json
import subprocess
import sys
from pathlib import Path

import yaml

REPO = Path(__file__).resolve().parents[2]
ROOT = REPO / "test/data"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    config = yaml.safe_load((REPO / "test/system/config.yaml").read_text())
    image = config["deps"]["trust4_container"].removeprefix("docker://")
    reference = args.reference.resolve()
    prefix = [str(REPO / "workflow/scripts/container_exec.sh"), "docker://" + image, image]
    settings = config["params"]["trust4"]
    for record in csv.DictReader((ROOT / "samplesheet.csv").open()):
        library = record["library_id"]
        target = args.output.resolve() / "libraries" / library
        assembly = target / "trust4"
        assembly.mkdir(parents=True)
        reads = {mate: [path.relative_to(REPO) for path in
                        sorted((ROOT / "reads" / library).glob(f"*_R{mate}_*.fastq.gz"))]
                 for mate in (1, 2)}
        read_args = [argument for mate in (1, 2) for path in reads[mate] for argument in (f"-{mate}", str(path))]
        command = prefix + ["run-trust4", "-f", str(reference), "--ref", str(reference), *read_args,
                            "-t", "2", "-k", str(settings["k_mer_size"]), "--contigMinCov",
                            str(settings["contig_min_cov"]), "-o", library, "--od", str(assembly), "--clean", "0"]
        with (target / "assembly.log").open("w") as log:
            subprocess.run(command, cwd=REPO, check=True, stdout=log, stderr=subprocess.STDOUT)
            context = dict(library_id=library,
                           image=dict(format="oci", digest="sha256:" + image.rsplit("sha256:", 1)[1],
                                      configured_reference=image),
                           settings={**settings, "clean": 0},
                           read_pairs=[dict(r1=str(first), r2=str(second))
                                       for first, second in zip(reads[1], reads[2])])
            captured = subprocess.check_output([sys.executable, str(REPO / "workflow/scripts/capture_assembly_context.py"),
                                               "--repo", str(REPO), "--context", json.dumps(context)], text=True)
            subprocess.run(prefix + ["perl", str(REPO / "workflow/scripts/record_assembly_provenance.pl"),
                "--context", captured, "--reference", str(reference),
                "--output", str(assembly / f"{library}_assembly_provenance.json"),
                "--airr", str(assembly / f"{library}_airr.tsv"), "--cdr3", str(assembly / f"{library}_cdr3.out"),
                "--final", str(assembly / f"{library}_final.out")], cwd=REPO, check=True, stdout=log, stderr=subprocess.STDOUT)
            metrics = target / "metrics"
            metrics.mkdir()
            subprocess.run(["bash", str(REPO / "workflow/scripts/count_seq_depth.sh"),
                str(metrics / f"{library}_seq_depth.tsv"), library,
                "--r1", *map(str, reads[1]), "--r2", *map(str, reads[2])], cwd=REPO,
                check=True, stdout=log, stderr=subprocess.STDOUT)
        print(f"Assembled {library}", flush=True)


if __name__ == "__main__":
    main()
