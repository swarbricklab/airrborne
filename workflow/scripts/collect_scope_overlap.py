"""Collect each supplied donor overlap product independently.

collect_scope_overlap() in scope_tables.py adds scope_id and keeps supplied
donor and row order. Producers own the pair, locus-summary and shared-clone
measurements. With no donor inputs, each product still receives its declared
header.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(snakemake.input.collection_module).resolve().parents[2]))
from workflow.lib.scope_tables import collect_scope_overlap

for product in ("donor_overlap_pairs", "donor_overlap_summary", "shared_groups"):
    collect_scope_overlap([Path(path) for path in snakemake.input[product]],
                           Path(snakemake.output[product]), str(snakemake.params.scope_id),
                           product)
Path(snakemake.log[0]).write_text("Collected the three donor overlap products\n")
