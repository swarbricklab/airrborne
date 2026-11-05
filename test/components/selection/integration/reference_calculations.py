"""Reconstruct selection results independently from retained scientific inputs."""

from pathlib import Path
from test.support.reference_tables import check_fields, number, read_table, index_rows
import math
from test.components.isotype.integration.reference_calculations import resolve_call


def verify_selection_estimates(rows: list[dict[str, str]]) -> int:
    """Require usable estimates and intervals wherever selection reports success."""
    for row in rows:
        for region in ("cdr", "fwr"):
            context = (row["library_id"], row["isotype_group"], region)
            used = int(row[f"n_clone_consensus_sequences_used_{region}"])
            assert 0 <= used <= int(row["n_clone_consensus_sequences"]), context
            values = [number(row[column]) for column in
                      (f"sigma_{region}", f"sigma_{region}_ci_low", f"sigma_{region}_ci_high")]
            if row[f"sigma_{region}_status"] == "estimated":
                assert used > 0 and all(value is not None and math.isfinite(value)
                                        for value in values), (context, "missing selection estimate")
                estimate, lower, upper = values
                assert lower <= estimate <= upper, (context, "invalid selection interval")
            else:
                assert all(value is None for value in values), (context, "unavailable selection estimate")
    return len(rows)


def verify_selection_audit(root: Path, scope: Path) -> int:
    """Rebuild selection eligibility populations and both support measures from members.

    The donor member shards keep the qualified C calls and alignments the
    released member table leaves out.
    """
    members = [row for path in sorted((scope / "intermediate").glob("*/*_clone_members.tsv"))
               for row in read_table(path)]
    audit = index_rows([row for row in read_table(scope / "metrics/qc_audit.tsv")
                        if row["stage"] in {"selection_input", "selection_input_source_support"}],
                       ("library_id", "isotype_group", "stage"))
    expected_keys = set()
    for context in read_table(scope / "manifest/scope_manifest.tsv"):
        library = context["library_id"]
        directory = root / "libraries" / library
        sources = {library + "~" + row["sequence_id"]: row for row in
                   read_table(directory / "airr" / f"{library}_airr_qc_pass.tsv")}
        local = [row for row in members if row["library_id"] == library]
        for group in ("all_IGH", "IgM", "IgG", "IgA"):
            def belongs(row):
                if group == "all_IGH":
                    return True
                return resolve_call(row["qualified_c_call"])[0] == group
            before = [row for row in local if belongs(row)]
            after = [row for row in before if all(row[column] != ""
                     for column in ("clone_id", "sequence_alignment", "germline_alignment_d_mask"))]
            for stage in ("selection_input", "selection_input_source_support"):
                key = (library, group, stage)
                expected_keys.add(key)
                def support(population):
                    return sum(int((sources[row["sequence_id"]] if stage.endswith("source_support")
                                    else row)["consensus_count"]) for row in population)
                check_fields(audit[key], dict(input_unit="inference_representative",
                    output_unit="selection_input_record", n_units_in=len(before), n_units_out=len(after),
                    consensus_count_in=support(before), consensus_count_out=support(after)), (scope.name, key))
    assert set(audit) == expected_keys
    return len(audit)
