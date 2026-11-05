"""Write airr/Repertoire.json linking this scope's shared library AIRR exports.

Reads the scope manifest, the scope provenance and effective configuration, and
each library's processing record. Every library must share one processing
identity. The log states the strict metadata result and the number of metadata
gaps. The publisher writes airr/validation.json from the staged release.
"""

import csv
import json
import sys
import yaml
from pathlib import Path

sys.path.insert(0, str(Path(snakemake.input.module).resolve().parents[2]))
from workflow.lib.airr_exports import AIRR_VERSION
from workflow.lib.repertoire_metadata import build_scope_repertoires, validate_repertoire_metadata


def software_description(entries: list[dict[str, str | None]]) -> str:
    """Name measured versions and explain each unavailable version."""
    return "; ".join(f"{entry['name']} {entry['version']}" if entry["status"] == "measured"
                     else f"{entry['name']} (version unavailable)" for entry in entries)


with open(snakemake.input.manifest, newline="") as handle:
    libraries = list(csv.DictReader(handle, delimiter="\t"))
identities = {}
for path in snakemake.input.identities:
    record = json.loads(Path(path).read_text())
    library_id = record["library_id"]
    if library_id in identities:
        raise ValueError(f"Duplicate processing record for library {library_id}: {path}")
    identities[library_id] = record
if set(identities) != {row["library_id"] for row in libraries}:
    raise ValueError("Processing records differ from scope libraries")
processing_ids = {record["processing_identity"] for record in identities.values()}
if len(processing_ids) != 1:
    raise ValueError("Scope library annotations were produced with different processing settings")
processing_identity = processing_ids.pop()
files = {library: f"airr/{library}_rearrangements.tsv.gz" for library in identities}
provenance = json.loads(Path(snakemake.input.provenance).read_text())
# The analysis image annotates the records; the release environment validates the
# exports against the pinned AIRR schema with its own airr-python.
software = (software_description(provenance["software"]["analysis"]["tools"])
            + f"; airr-python {AIRR_VERSION} (export validation)")
history = provenance["assembly_producers"]
library_processing = {}
for library in libraries:
    library_id = library["library_id"]
    producer_id = history["libraries"].get(library_id)
    producer = history["producers"][producer_id] if producer_id is not None else None
    assembly_software = (software_description(producer["software"]) if producer is not None else "")
    library_processing[library_id] = {
        "software_versions": "; ".join(value for value in (
            assembly_software or "TRUST4 (producer version unavailable)", software) if value),
        "paired_reads_assembly": (
            "TRUST4 assembled paired-end reads; recorded producer settings are in provenance.json."
            if producer is not None else
            "Reused TRUST4 outputs; assembly version and settings are unavailable because no producer record was supplied."),
    }
effective = yaml.safe_load(Path(snakemake.input.effective_config).read_text())
record = build_scope_repertoires(
    str(snakemake.wildcards.scope), libraries, files, processing_identity,
    library_processing_fields=library_processing,
    processing_fields={
        "quality_thresholds": "No separate FASTQ base-quality filter is run by airrborne. Upstream filtering is not documented in the inputs.",
        "primer_match_cutoffs": "No primer sequences or match cutoffs are supplied to airrborne; no separate primer trimming is run.",
        "collapsing_method": (
            "Every common-QC TRUST4 record is exported as its own row. Grouping-eligible records of one "
            "IGH, IGK or IGL assembled contig collapse to one contig representative "
            "(is_contig_representative); at TRA, TRB, TRD and TRG each eligible record represents itself. "
            "Contig representatives with an identical full sequence within a library and locus merge into "
            "one member (is_member_representative), whose member_consensus_count sums the eligible support "
            "transferred to it; every other row carries zero there. clone_id names the member's IGH clone, "
            "or its clonotype at the other loci, in this scope, and is empty on every row that is not a "
            "member representative."),
        "data_processing_protocols": (
            "Imputed TRUST4 CDR3 records are excluded before IgBLAST/MakeDb annotation. "
            "Common QC requires productive, locus-consistent records with positive read support, "
            "junction length divisible by three and alignment N fraction at most "
            + str(effective["params"]["airr_qc"]["max_alignment_n_fraction"])
            + ". A record enters grouping when its junction is valid, its V and J calls each resolve to one "
            "germline entity and its constant-region call is consistent with its locus. At IGH, IGK and IGL "
            "the dominant common-QC variant of the record's contig, chosen by original read support, "
            "must also meet these requirements, and the record's support then transfers to that "
            "contig representative. grouping_status gives each record's "
            "reason. effective_config.yaml records the settings."),
        "germline_database": ("IMGT/GENE-DB download dated "
                              + provenance["references"]["imgt_reference"]["imgt_snapshot"]["date"]
                              + " in the analysis image; curated reference SHA256 "
                              + provenance["references"]["imgt_reference"]["sha256"])})
output = Path(snakemake.output.repertoire)
output.parent.mkdir(parents=True, exist_ok=True)
validation = validate_repertoire_metadata(record, libraries, files, processing_identity)
temporary = output.with_suffix(".tmp")
temporary.write_text(json.dumps(record, indent=2, ensure_ascii=True, allow_nan=False) + "\n")
temporary.replace(output)
Path(snakemake.log[0]).write_text(
    f"Linked {len(libraries)} libraries; strict metadata valid: {validation['strict_valid']}; "
    f"strict metadata gaps: {len(validation['strict_metadata_gaps'])}\n")
