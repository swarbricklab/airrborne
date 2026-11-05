"""Collect supplied library-by-isotype selection summaries.

collect_scope_selection() in scope_tables.py sorts rows by textual library_id
and keeps each producer's isotype order within a library. Producers own the
estimates and their availability. Header-only inputs retain the output schema.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(snakemake.input.collection_module).resolve().parents[2]))
from workflow.lib.scope_tables import collect_scope_selection

collect_scope_selection([Path(path) for path in snakemake.input.selection],
                         Path(snakemake.output.selection))
Path(snakemake.log[0]).write_text("Collected donor selection summaries\n")
