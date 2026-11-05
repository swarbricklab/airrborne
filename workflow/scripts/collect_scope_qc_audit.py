"""Collect the scope's supplied population counts and recovered read support.

collect_scope_qc_audit() in scope_qc_audit.py validates producer counts, joins
manifest context and writes the audit, QC and QC trail tables. The producer
records supplied by the rule identify known assembly history; other manifest
libraries receive unknown_history. The QC table reads locus, support and
grouping status from each library's AIRR export. This script applies no record
QC.
"""
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(snakemake.input.module).resolve().parents[2]))
from workflow.lib.scope_qc_audit import EXPORT_SUFFIX, collect_scope_qc_audit

known = [json.loads(Path(path).read_text())["library_id"] for path in snakemake.input.producers]
exports = {Path(path).name.removesuffix(EXPORT_SUFFIX): Path(path) for path in snakemake.input.exports}
collect_scope_qc_audit(
    Path(snakemake.input.manifest), list(map(Path, snakemake.input.shards)),
    list(map(Path, snakemake.input.diversity)), Path(snakemake.input.selection),
    Path(snakemake.output.audit), known,
    exports=exports, qc_output=Path(snakemake.output.qc), trail_output=Path(snakemake.output.trail),
    value_orders=snakemake.params.value_orders, coverage_file=Path(snakemake.input.coverage))
Path(snakemake.log[0]).write_text("Validated scope population counts and canonical-locus totals\n")
