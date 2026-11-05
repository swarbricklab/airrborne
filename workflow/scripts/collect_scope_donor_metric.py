"""Merge the declared sorted donor measurement tables with bounded memory.

collect_donor_tables() in scope_tables.py checks the rule's key_columns, in the
value order key_orders lists for a column or as text otherwise, and rejects
repeated or unsorted keys. It leaves out rows holding dropped_values, empties
each gated_columns cell whose eligibility column is not TRUE and writes the
output_columns when the rule lists them. Producers define the row units,
estimates and eligibility. Header-only inputs produce a header-only output.
Collection writes a temporary file and replaces the destination only after all
checks pass.
"""

import csv
import gzip
import sys
from pathlib import Path

sys.path.insert(0, str(Path(snakemake.input.collection_module).resolve().parents[2]))
from workflow.lib.scope_tables import collect_donor_tables

statuses = None
if snakemake.wildcards.donor_metric == "clone_member_shm_coverage":
    statuses = {}
    for path in snakemake.input.exports:
        with gzip.open(path, "rt", newline="") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                if row["locus"] == "IGH" and row["germline_reconstruction_status"] != "available":
                    statuses[row["sequence_id"]] = row["germline_reconstruction_status"]
collect_donor_tables([Path(path) for path in snakemake.input.tables],
                     Path(snakemake.output.table), list(snakemake.params.key_columns),
                     dict(snakemake.params.key_orders), dict(snakemake.params.gated_columns),
                     dict(snakemake.params.output_columns), dict(snakemake.params.dropped_values), statuses)
Path(snakemake.log[0]).write_text(f"Merged {len(snakemake.input.tables)} donor shards\n")
