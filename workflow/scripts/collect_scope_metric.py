"""Stream supplied library measurements and attach manifest context.

collect_scope_metric() in scope_tables.py applies the rule's context_columns,
quoting and optional sort_columns and label_columns settings without
recalculating measurements. Label columns name V/J entities from the reference
metadata input. It preserves source rows and textual identifiers, accepts
header-only tables and rejects missing tables or inconsistent headers before
replacing the output.
"""

import csv
import sys
from pathlib import Path

sys.path.insert(0, str(Path(snakemake.input.collection_module).resolve().parents[2]))
from workflow.lib.scope_tables import collect_scope_metric, read_entity_labels

label_columns = dict(getattr(snakemake.params, "label_columns", {}) or {})
reference = list(getattr(snakemake.input, "reference_metadata", []) or [])
collect_scope_metric(Path(snakemake.input.manifest),
                     [Path(path) for path in snakemake.input.metric],
                     Path(snakemake.output.metric),
                     tuple(snakemake.params.context_columns), snakemake.params.quoting,
                     tuple(getattr(snakemake.params, "sort_columns", ())),
                     label_columns, read_entity_labels(Path(reference[0])) if reference else None)
Path(snakemake.log[0]).write_text(f"Collected {len(snakemake.input.metric)} measurement shards\n")
