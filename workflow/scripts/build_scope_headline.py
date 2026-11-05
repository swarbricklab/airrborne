"""Assemble one scientific headline row per manifest library.

build_headline() in headline.py joins the authoritative measurements and checks
required source rows. Unavailable estimates remain missing. write_headline()
atomically writes the table; project annotations enter only report preparation.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(snakemake.input.headline_module).resolve().parents[2]))
from workflow.lib.headline import build_headline, write_headline
rows = build_headline(
    snakemake.input.manifest, snakemake.input.yields, snakemake.input.qc,
    snakemake.input.diversity, snakemake.input.isotype, snakemake.input.isotype_summary, snakemake.input.shm,
    scope_id=snakemake.wildcards.scope)
write_headline(rows, Path(snakemake.output.headline))
Path(snakemake.log[0]).write_text(f"Wrote {len(rows)} library headline rows\n")
