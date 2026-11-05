"""Literal release context for publication and rendering contract fixtures."""

import csv
import json
from pathlib import Path

import yaml


def write_release_context(scope: Path) -> tuple[dict, Path, Path, Path]:
    """Write default settings and constructed provenance for a populated scope.

    The manifest supplies the scope and library IDs. The image digests and
    fallback threshold are fixture values, not evidence of a workflow run.
    """
    with (scope / "manifest/scope_manifest.tsv").open(newline="") as handle:
        libraries = list(csv.DictReader(handle, delimiter="\t"))
    scope_id = libraries[0]["scope_id"]
    params = yaml.safe_load((Path(__file__).resolve().parents[2] / "config/airrborne/config.yaml").read_text())["params"]
    config = {"schema_version": "3.2", "scope_id": scope_id, "seed": 1337, "params": params}
    provenance = {
        "schema_version": "3.2", "scope_id": scope_id,
        "source": {"revision": "fixture", "dirty": False},
        "images": {"analysis": {"format": "docker", "digest": "sha256:" + "a" * 64},
                   "assembly": {"format": "docker", "digest": "sha256:" + "b" * 64}},
        "software": {"analysis": {"image": "analysis", "tools": []},
                     "configured_assembly": {"image": "assembly", "tools": []}},
        "references": {},
        # provenance.py names the scope working-tree paths of these files.
        "files": {"manifest": "manifest/scope_manifest.tsv", "effective_config": "effective_config.yaml",
                  "threshold": f"clones/{scope_id}_clone_threshold.tsv",
                  "donor_thresholds": f"clones/{scope_id}_donor_clone_thresholds.tsv",
                  "distance_diagnostics": f"clones/{scope_id}_clone_distance_diagnostics.tsv"},
        "assembly_producers": {"libraries": {}, "producers": {},
                               "unknown": {"libraries": [row["library_id"] for row in libraries],
                                           "reason": "No fixture producer record was supplied."}},
    }
    ppath, cpath, threshold = scope / "provenance.json", scope / "effective_config.yaml", scope / "threshold.tsv"
    ppath.write_text(json.dumps(provenance))
    cpath.write_text(yaml.safe_dump(config))
    threshold.write_text("scope_id\tclone_threshold\tthreshold_source\n" + scope_id + "\t0.1\tfallback\n")
    return provenance, ppath, cpath, threshold
