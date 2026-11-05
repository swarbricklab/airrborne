"""Reconstruct isotype results independently from retained scientific inputs."""

from pathlib import Path
from test.support.reference_tables import check_fields, number, read_table, index_rows
from collections import Counter
import math
import re
CLASSES = ("IgM", "IgD", "IgG", "IgA", "IgE")
SUBCLASSES = ("IgM", "IgD", "IgG1", "IgG2", "IgG3", "IgG4", "IgA1", "IgA2", "IgE")


def resolve_call(call: str) -> tuple[str | None, str | None]:
    genes = {token.strip().split("*", 1)[0] for token in call.split(",") if token.strip()}
    genes = {"IGHG4" if gene == "IGHG4A" else gene for gene in genes}
    labels = {"IGH" + label[2:]: label for label in SUBCLASSES}
    if not genes or not genes <= labels.keys():
        return None, None
    subclasses = {labels[gene] for gene in genes}
    classes = {re.sub(r"[1-4]$", "", label) for label in subclasses}
    return (next(iter(classes)) if len(classes) == 1 else None,
            next(iter(subclasses)) if len(subclasses) == 1 else None)


def verify_isotypes(root: Path, scope: Path) -> tuple[int, int, int]:
    metrics = scope / "metrics"
    composition = index_rows(read_table(metrics / "isotype_usage.tsv"), ("library_id", "isotype_level", "isotype"))
    summaries = index_rows(read_table(metrics / "isotype_summary.tsv"), ("library_id",))
    expected_keys = set()
    class_only_total = 0
    for library, in summaries:
        directory = root / "libraries" / library
        records = [row for row in read_table(directory / "airr" / f"{library}_airr_qc_pass.tsv") if row["locus"] == "IGH"]
        calls = index_rows(read_table(directory / "metrics" / f"{library}_isotype_calls.tsv"), ("sequence_id",))
        assert set(calls) == {(row["sequence_id"],) for row in records}
        support = {"class": Counter(), "subclass": Counter()}
        total = 0
        class_only = 0
        for row in records:
            count = int(row["consensus_count"])
            call = calls[(row["sequence_id"],)]
            major, minor = resolve_call(call["qualified_c_call"])
            # The sidecar names the same isotype for the AIRR export.
            check_fields(call, {"isotype_class": major, "isotype_subclass": minor},
                         (library, row["sequence_id"]))
            total += count
            if major is not None:
                support["class"][major] += count
                if minor is None:
                    class_only += count
            if minor is not None:
                support["subclass"][minor] += count
        for resolution, levels in (("class", CLASSES), ("subclass", SUBCLASSES)):
            callable_count = sum(support[resolution].values())
            for label in levels:
                key = (library, resolution, label)
                expected_keys.add(key)
                count = support[resolution][label]
                check_fields(composition[key], {
                    "consensus_count": count, "callable_consensus_count": callable_count,
                    "support_fraction": count / callable_count if callable_count else None}, key)
        major = support["class"]
        callable_count = sum(major.values())
        igg, iga = major["IgG"], major["IgA"]
        status = ("no_class_callable_evidence" if not callable_count else
                  "both_zero" if not igg and not iga else "zero_igg" if not igg else
                  "zero_iga" if not iga else "available")
        check_fields(summaries[(library,)], {
            "common_qc_consensus_count": total,
            "class_callable_consensus_count": callable_count,
            "class_uncallable_consensus_count": total - callable_count,
            "subclass_callable_consensus_count": sum(support["subclass"].values()),
            "class_only_consensus_count": class_only,
            "igm_consensus_count": major["IgM"], "igg_consensus_count": igg,
            "iga_consensus_count": iga, "igg_iga_status": status,
            "log2_igg_iga": math.log2(igg) - math.log2(iga) if igg and iga else None}, library)
        class_only_total += class_only
    assert set(composition) == expected_keys
    return len(composition), len(summaries), class_only_total
