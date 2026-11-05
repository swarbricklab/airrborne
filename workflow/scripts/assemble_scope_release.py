"""Publish the complete scope handover through one recoverable transaction.

Collects the release allow-list from the scope directory and each library's AIRR
export with this scope's member assignments, and requires every library to share
one processing identity. The staged
folder's airr/validation.json is written from the staged files before the
folder replaces the release. The receipt .complete.json is written last, so its
presence marks a finished publication. It holds the summary
build_release_bundle() returns, including the hash of manifest/release_files.tsv.
"""

import csv
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(snakemake.input.module).resolve().parents[2]))
from workflow.lib.release_bundle import (
    build_release_bundle, release_path, scope_release_source_files, table_for)
from workflow.lib.alignment_export import read_reference
from workflow.lib.repertoire_metadata import validate_scope_repertoires

scope = str(snakemake.wildcards.scope)
source = Path(snakemake.params.scope_dir)
paths = scope_release_source_files(scope, bool(snakemake.params.membership))
generated = {"provenance.json", "effective_config.yaml"}
files = [relative for relative in paths if relative != "report" and relative not in generated]
airr = {library: {Path(path): f"airr/{library}_rearrangements.tsv.gz"}
        for library, path in zip(snakemake.params.libraries, snakemake.input.airr)}
tables = {release_path(relative): table_for(Path(release_path(relative)))
          for relative in files if relative.endswith(".tsv")}
tables.update({relative: "exported_AIRR" for sources in airr.values() for relative in sources.values()})
identities = [json.loads(Path(path).read_text()) for path in snakemake.input.identities]
if {row["library_id"] for row in identities} != set(airr):
    raise ValueError("Processing records differ from the released libraries")
processing_ids = {row["processing_identity"] for row in identities}
if len(processing_ids) != 1:
    raise ValueError("Released libraries have different processing settings")
processing_identity = processing_ids.pop()
with (source / "manifest/scope_manifest.tsv").open(newline="") as handle:
    manifest = list(csv.DictReader(handle, delimiter="\t"))


def validate_staged(staged: Path) -> None:
    """Validate the actual staged AIRR links before replacing the prior release."""
    record = json.loads((staged / "airr/Repertoire.json").read_text())
    files = {library: next(iter(sources.values())) for library, sources in airr.items()}
    observed = validate_scope_repertoires(record, staged, manifest, files, processing_identity,
                                          member_assignments=True,
                                          reference=read_reference(Path(snakemake.input.reference)))
    (staged / "airr/validation.json").write_text(
        json.dumps(observed, indent=2, sort_keys=True, allow_nan=False) + "\n")


summary = build_release_bundle(
    source, Path(snakemake.params.bundle),
    scope_files=files, scope_directories=["report"], airr_files=airr,
    dictionary=Path(snakemake.input.dictionary), table_names=tables,
    provenance=source / "provenance.json", effective_config=source / "effective_config.yaml",
    calibration=source / f"clones/{scope}_clone_threshold.tsv",
    readme_template=Path(snakemake.input.readme),
    validate_staged=validate_staged)
receipt = Path(snakemake.output.receipt)
temporary = receipt.with_suffix(".tmp")
temporary.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
temporary.replace(receipt)
Path(snakemake.log[0]).write_text(json.dumps(summary, sort_keys=True) + "\n")
