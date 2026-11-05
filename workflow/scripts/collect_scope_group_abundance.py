"""Collect retained group support from the IGH and non-IGH repertoire arms.

collect_scope_group_abundance() in scope_tables.py adds manifest context and
merges the per-library non-IGH tables with the scope's IGH table with bounded
memory. It writes rows by library, then locus, keeps source order within each
locus and retains a header for header-only tables. The producer tables define
the groups and their recovered read support.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(snakemake.input.collection_module).resolve().parents[2]))
from workflow.lib.scope_tables import collect_scope_group_abundance

collect_scope_group_abundance(Path(snakemake.input.manifest),
    [Path(path) for path in snakemake.input.nonigh], Path(snakemake.input.igh),
    Path(snakemake.output.group_abundance), str(snakemake.params.scope_id))
Path(snakemake.log[0]).write_text("Collected library clone-abundance vectors\n")
