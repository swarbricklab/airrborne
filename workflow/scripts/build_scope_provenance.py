"""Write one scope's provenance.json and effective_config.yaml.

Reads the scope manifest, the measured analysis and assembly tool versions, the
curated reference, its auxiliary data and reference receipt, and each library's
assembly outputs with its producer record when one exists. A library without a
producer record is recorded with unknown history. The pipeline source state,
its patches and its untracked files are kept under the source-evidence
directory named by the source-state hash. effective_config.yaml holds the
configured settings and, under fixed_methods, the scientific choices the code
hard-codes.
"""

import csv
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(snakemake.input.module).resolve().parents[2]))
from workflow.lib.provenance import (build_scope_provenance,
    effective_scope_config, load_assembly_provenance, write_scope_provenance, write_source_evidence)
from workflow.lib.release_methods import FIXED_METHODS

with open(snakemake.input.manifest, newline="") as handle:
    libraries = list(csv.DictReader(handle, delimiter="\t"))
versions = {role: json.loads(Path(path).read_text()) for role, path in
            (("analysis", snakemake.input.analysis_versions),
             ("assembly", snakemake.input.assembly_versions))}
assembly = {}
for library in libraries:
    name = library["library_id"]
    directory = Path(snakemake.params.results_dir) / "libraries" / name / "trust4"
    sidecar = directory / (name + "_assembly_provenance.json")
    outputs = {key: directory / (name + suffix) for key, suffix in
               (("airr", "_airr.tsv"), ("cdr3", "_cdr3.out"), ("final_out", "_final.out"))}
    assembly[name] = load_assembly_provenance(name, sidecar if sidecar.exists() else None, outputs)
record = build_scope_provenance(
    snakemake.wildcards.scope, libraries, snakemake.params.effective_config,
    source=snakemake.params.source,
    images={role: value["image"] for role, value in versions.items()},
    tool_versions=versions,
    references={"imgt_reference": Path(snakemake.input.reference),
                "reference_metadata": Path(snakemake.input.reference_metadata),
                "auxiliary_data": Path(snakemake.input.auxiliary)},
    reference_receipt=json.loads(Path(snakemake.input.reference_receipt).read_text()),
    assembly_records=assembly)
write_source_evidence(Path(snakemake.params.repo), Path(snakemake.params.source_evidence),
                      snakemake.params.source)
effective = {**effective_scope_config(snakemake.params.effective_config, str(snakemake.wildcards.scope)),
             "fixed_methods": FIXED_METHODS}
write_scope_provenance(record, effective,
                       Path(snakemake.output.provenance),
                       Path(snakemake.output.config))
Path(snakemake.log[0]).write_text(f"Recorded provenance for {len(libraries)} libraries\n")
