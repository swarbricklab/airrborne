"""Write and validate every common-QC AIRR record for one library.

read_library_annotations() joins the grouping and qualified-call decisions from
their producers, and add_makedb_fields() the MakeDb germline coordinates and
V-domain regions. export_library_airr() preserves source records and adds
library-qualified links; it does not select inference members. Processing
identity comes from library_processing_identity() using the declared code,
settings, image and resources. validate_export_content() checks the written
export before this script atomically replaces the identity record.
"""

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(snakemake.input.export_module).resolve().parents[2]))
from workflow.lib.airr_exports import (
    add_makedb_fields, export_library_airr, read_library_annotations, validate_export_content)
from workflow.lib.alignment_export import add_alignment_evidence, read_reference, validate_alignment
from workflow.lib.airr_exports import iter_library_airr
from workflow.lib.provenance import library_processing_identity

library = str(snakemake.wildcards.library_id)
image = json.loads(Path(snakemake.input.analysis_image).read_text())
root = Path(snakemake.input.export_module).resolve().parents[2]
identity = library_processing_identity(
    {str(Path(path).resolve().relative_to(root)): Path(path) for path in snakemake.input.code},
    snakemake.params.settings, image, {"imgt_reference": Path(snakemake.input.reference), "auxiliary_data": Path(snakemake.input.auxiliary),
     "reference_metadata": Path(snakemake.input.reference_metadata),
     "airr_schema": Path(snakemake.input.schema)})
annotations, fields = read_library_annotations(
    [Path(snakemake.input.igh_grouping), Path(snakemake.input.nonigh_grouping)],
    Path(snakemake.input.isotype_calls), library)
fields = add_makedb_fields(annotations, fields, Path(snakemake.input.makedb))
fields = add_alignment_evidence(annotations, fields, Path(snakemake.input.igblast_hits))
# The native MakeDb CIGAR predates the producer I/D adapter.
import csv
import gzip
with gzip.open(snakemake.input.makedb, "rt", newline="") as handle:
    for row in csv.DictReader(handle, delimiter="\t"):
        if row["sequence_id"] in annotations:
            annotations[row["sequence_id"]]["changeo_d_cigar"] = row.get("d_cigar", "")
ids = export_library_airr(Path(snakemake.input.airr), Path(snakemake.output.airr),
                         library, identity, annotations=annotations, annotation_fields=fields)
count = validate_export_content(Path(snakemake.input.airr), Path(snakemake.output.airr),
                               library, identity, annotations=annotations, annotation_fields=fields)
reference = read_reference(Path(snakemake.input.reference))
for row in iter_library_airr(Path(snakemake.output.airr), library, identity, reference=reference):
    pass
record = {"library_id": library, "processing_identity": identity, "records": count, **ids}
output = Path(snakemake.output.identity)
temporary = output.with_suffix(".tmp")
temporary.write_text(json.dumps(record, indent=2, sort_keys=True) + "\n")
temporary.replace(output)
Path(snakemake.log[0]).write_text(f"Validated all {count} source records and export links\n")
