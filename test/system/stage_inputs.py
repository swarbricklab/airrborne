#!/usr/bin/env python3
"""Prepare a fresh fixture run, retaining assembly evidence when supplied."""

import argparse
import csv
import shutil
from pathlib import Path

import yaml


def prepare(config_path: Path, work_dir: Path, out_config: Path,
            seed_results: Path | None = None, resources: Path | None = None) -> dict:
    """Write isolated output paths and copy every manifest library's native bundle.

    Omitting seed_results prepares a FASTQ run with no retained assemblies.
    Input paths are relative to the repository working directory, as they are
    when Snakemake reads the original config. resources replaces deps.resources
    with a directory that setup.sh already prepared, so a staged run reuses
    its images, receipts and reference.
    """
    work_dir = work_dir.resolve()
    actual = work_dir / "results" / "airrborne"
    logs = work_dir / "logs" / "airrborne"
    if actual.exists() or logs.exists():
        raise ValueError(f"Fixture outputs already exist under {work_dir}")
    if out_config.resolve() == config_path.resolve():
        raise ValueError("The isolated config must differ from the input config")
    config = yaml.safe_load(config_path.read_text())
    for key in ("samplesheet", "data", "resources"):
        config["deps"][key] = str(Path(config["deps"][key]).resolve())
    if resources is not None:
        config["deps"]["resources"] = str(Path(resources).resolve())
    if config.get("membership_csv"):
        config["membership_csv"] = str(Path(config["membership_csv"]).resolve())
    for sources in config.get("scope_metadata", {}).values():
        for source in sources.values():
            source["file"] = str(Path(source["file"]).resolve())
    config["outs"]["results_dir"] = str(actual)
    config["outs"]["log_dir"] = str(logs)

    seeds = []
    if seed_results is not None:
        seed_results = seed_results.resolve()
        with Path(config["deps"]["samplesheet"]).open(newline="") as handle:
            libraries = [row["library_id"].strip() for row in csv.DictReader(handle)]
        if not libraries or len(set(libraries)) != len(libraries):
            raise ValueError("The fixture samplesheet must name distinct libraries")
        for library in libraries:
            if not library or library in (".", "..") or "/" in library or "\\" in library:
                raise ValueError(f"Invalid fixture library_id: {library!r}")
            library_dir = Path("libraries") / library
            seeds.extend(library_dir / "trust4" / f"{library}_{suffix}" for suffix in (
                "airr.tsv", "cdr3.out", "final.out", "assembly_provenance.json"))
            seeds.append(library_dir / "metrics" / f"{library}_seq_depth.tsv")
        missing = [str(path) for path in seeds if not (seed_results / path).is_file()]
        if missing:
            raise ValueError(f"Incomplete fixture under {seed_results}:\n  " + "\n  ".join(missing))

    # Validate the whole bundle before writing any seed. copy2 retains the
    # original producer record verbatim and keeps the frozen inputs independent.
    actual.mkdir(parents=True)
    for relative in seeds:
        destination = actual / relative
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(seed_results / relative, destination)
    out_config.parent.mkdir(parents=True, exist_ok=True)
    out_config.write_text(yaml.safe_dump(config, sort_keys=False))
    print(f"Staged {len(seeds)} native seed files in {actual}")
    print(f"Wrote isolated config to {out_config}")
    return config


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--seed-results", type=Path,
                        help="Native bundles for post-TRUST4 execution; omit for a FASTQ run.")
    parser.add_argument("--work-dir", required=True, type=Path)
    parser.add_argument("--out-config", required=True, type=Path)
    parser.add_argument("--resources", type=Path,
                        help="Prepared resources directory to use instead of the configuration's.")
    args = parser.parse_args()
    try:
        prepare(args.config, args.work_dir, args.out_config, args.seed_results, args.resources)
    except (OSError, ValueError, KeyError) as error:
        parser.exit(1, f"{error}\n")


if __name__ == "__main__":
    main()
