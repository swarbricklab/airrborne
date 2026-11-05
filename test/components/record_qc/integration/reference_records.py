"""Check reviewed record decisions and the complete early QC audit independently.

The small connected fixture has no imputed calls, exact duplicates or annotation
loss. Boundary and real-tool tests exercise those branches separately. The
reviewed table anchors decisions here instead of reusing production predicates.
"""

import csv
from collections import Counter
from itertools import combinations
from pathlib import Path

from test.support.reference_tables import check_fields, read_table

LOCI = ("IGH", "IGK", "IGL", "TRA", "TRB", "TRD", "TRG")
CRITERIA = ("locus_match", "d_locus_match", "low_n", "productive",
            "junction_mod3", "positive_support")


def keyed(rows: list[dict], columns: tuple[str, ...]) -> dict:
    result = {tuple(row[column] for column in columns): row for row in rows}
    assert len(result) == len(rows), f"Repeated keys: {columns}"
    return result


def gene_set(call: str) -> set[str]:
    return {token.split("*", 1)[0] for token in call.split(",") if token}


def expected_early_rows(records: list[dict], cdr3: list[list[str]],
                        max_n_fraction: float) -> list[dict]:
    """Specify populations, units and statuses from reviewed fixture decisions."""
    result = []

    def add(stage, index, locus, incoming, outgoing, *, kind="checkpoint",
            basis="native", analysis="common", unit="rearrangement",
            conserved=True, status="applied", support=True):
        result.append(dict(stage=stage, stage_index=index, locus=locus,
            locus_basis=basis, analysis=analysis, grain="library", kind=kind,
            input_unit=unit, output_unit=unit, decision_status=status,
            support_conserved=conserved, isotype_group="", v_region="",
            n_units_in=len(incoming), n_units_out=len(outgoing),
            consensus_count_in=sum(int(row["consensus_count"]) for row in incoming) if support else None,
            consensus_count_out=sum(int(row["consensus_count"]) for row in outgoing) if support else None))

    contigs = set(row[0] for row in cdr3)
    for stage, index, population, unit in (
            ("trust4_contigs_with_cdr3_call", 1, contigs, "contig"),
            ("trust4_cdr3_calls", 2, cdr3, "cdr3_call")):
        add(stage, index, "all", population, population, basis="not_applicable",
            analysis="assembly", unit=unit, conserved=None, support=False)
    for locus in (*LOCI, "all"):
        population = [row for row in records if locus == "all" or row["locus"] == locus]
        for stage, index, kind, conserved in (
                ("trust4_reported_airr", 3, "checkpoint", True),
                ("imputed_cdr3_filter", 4, "eligibility", False),
                ("igblast_makedb", 7, "eligibility", False),
                ("igblast_locus_agreement", 8, "diagnostic", False)):
            add(stage, index, locus, population, population, kind=kind,
                conserved=conserved, basis="annotated" if index == 8 else "native")
        # Each criterion uses the same annotated population. Failures can overlap,
        # so these rows are not sequential losses.
        for index, criterion in enumerate(CRITERIA, 11):
            kept = [row for row in population if criterion not in row["qc_failures"].split(",")]
            add("qc_" + criterion, index, locus, population, kept,
                basis="annotated", kind="criterion", conserved=False,
                status="disabled" if criterion == "low_n" and max_n_fraction == 1 else "applied")
        kept = [row for row in population if not row["qc_failures"]]
        add("qc_pass_total", 20, locus, population, kept, basis="annotated",
            kind="eligibility", conserved=False)
    add("igblast_locus_routing", 6, "all", records, records, kind="eligibility", conserved=False)
    return result


def check_audit(actual: list[dict], expected: list[dict], context: str) -> int:
    observed = keyed(actual, ("stage", "locus"))
    wanted = keyed(expected, ("stage", "locus"))
    assert observed.keys() == wanted.keys(), f"{context}: incomplete or unexpected early audit stages/loci"
    for key, row in wanted.items():
        check_fields(observed[key], row, (context, key))
    return len(wanted)


def check_clone_partitions(members: dict, expected: dict, context: str) -> None:
    """Compare donor-specific equivalence classes without fixing clone labels."""
    assert members.keys() == expected.keys(), (context, "clone member set")
    for key, row in members.items():
        check_fields(row, {"member_consensus_count": expected[key][1]}, (context, key))
    for first, second in combinations(members, 2):
        if members[first]["donor_id"] == members[second]["donor_id"]:
            assert (members[first]["clone_id"] == members[second]["clone_id"]) == (
                expected[first][0] == expected[second][0]), (context, first, second, "clone partition")


def verify_records(root: Path, expectations: Path, populations: dict[str, set[str]],
                   config: dict) -> dict[str, int]:
    truth = read_table(expectations)
    libraries = set().union(*populations.values())
    assert libraries <= {row["library_id"] for row in truth}, "Missing reviewed libraries"
    expected_audits = {}
    expected_members = {}
    totals = Counter()
    for library in sorted(libraries):
        records = [row for row in truth if row["library_id"] == library]
        expected = keyed(records, ("sequence_id",))
        directory = root / "libraries" / library
        native = keyed(read_table(directory / "trust4" / f"{library}_airr.tsv"), ("sequence_id",))
        annotated = keyed(read_table(directory / "airr" / f"{library}_igblast.tsv"), ("sequence_id",))
        passed = keyed(read_table(directory / "airr" / f"{library}_airr_qc_pass.tsv"), ("sequence_id",))
        assert native.keys() == expected.keys() == annotated.keys(), (library, "native/annotation record set")
        assert passed.keys() == {key for key, row in expected.items() if not row["qc_failures"]}, (library, "common QC record set")
        assert len({row["sequence"] for row in native.values()}) == len(native), "Review exact duplicate outcomes"
        with (directory / "trust4" / f"{library}_cdr3.out").open() as handle:
            cdr3 = list(csv.reader(handle, delimiter="\t"))
        assert all(abs(float(row[9]) - .01) >= 1e-6 for row in cdr3), "Review imputed outcomes"
        for key, wanted in expected.items():
            for actual in (native[key], annotated[key], *([passed[key]] if key in passed else [])):
                check_fields(actual, {"locus": wanted["locus"], "consensus_count": int(wanted["consensus_count"])}, (library, key))
            for segment in ("v", "j"):
                assert gene_set(annotated[key][segment + "_call"]) == {wanted[segment + "_gene"]}, (library, key, segment)
            assert (annotated[key]["productive"] in ("T", "TRUE")) == ("productive" not in wanted["qc_failures"]), (library, key, "productivity")
        calls = keyed(read_table(directory / "metrics" / f"{library}_isotype_calls.tsv"), ("sequence_id",))
        assert calls.keys() == {key for key in passed if expected[key]["locus"] == "IGH"}
        for key, call in calls.items():
            gene = expected[key]["qualified_c_gene"]
            assert gene_set(call["qualified_c_call"]) == ({gene} if gene else set()), (library, key, "qualified C")
            check_fields(call, {"c_call_callable": bool(gene), "c_call_enough_evidence": bool(gene),
                                "c_call_present": True, "c_call_locus_match": True}, (library, key))
        mapping = keyed([row for branch in ("igh", "nonigh") for row in read_table(
            directory / "metrics" / f"{library}_{branch}_grouping.tsv")], ("sequence_id",))
        assert mapping.keys() == passed.keys(), (library, "grouping record set")
        support = Counter()
        for key in passed:
            support[expected[key]["representative"]] += int(expected[key]["consensus_count"])
        for key, row in mapping.items():
            representative = expected[key]["representative"]
            is_representative = key[0] == representative
            check_fields(row, {"grouping_eligible": True, "grouping_status": "eligible",
                "backbone_representative_sequence_id": representative,
                "is_backbone_representative": is_representative,
                "inference_representative_sequence_id": representative,
                "is_inference_representative": is_representative,
                "inference_consensus_count": support[key[0]] if is_representative else 0}, (library, key))
            if is_representative and expected[key]["locus"] == "IGH":
                expected_members[(library, library + "~" + key[0])] = (
                    expected[key]["clone_partition"], support[key[0]])
        rows = expected_early_rows(records, cdr3, config["params"]["airr_qc"]["max_alignment_n_fraction"])
        expected_audits[library] = rows
        actual = [row for shard in ("reannotate", "qc") for row in read_table(
            directory / "qc" / f"{library}_qc_audit_{shard}.tsv")]
        shard_expected = [dict(row, library_id=library, scope_id="", donor_id="", tissue_id="") for row in rows]
        totals["early_audit_rows"] += check_audit(actual, shard_expected, library)
        totals["reviewed_records"] += len(records)
    for name, included in populations.items():
        scope = root / "scopes" / name
        manifest = keyed(read_table(scope / "manifest/scope_manifest.tsv"), ("library_id",))
        audit = read_table(scope / "metrics/qc_audit.tsv")
        for library in included:
            context = manifest[(library,)]
            wanted = [dict(row, library_id=library, scope_id=name,
                           donor_id=context["donor_id"], tissue_id=context["tissue_id"])
                      for row in expected_audits[library]]
            actual = [row for row in audit if row["library_id"] == library and 0 < int(row["stage_index"]) <= 20]
            totals["scope_early_audit_rows"] += check_audit(actual, wanted, (name, library))
        members = keyed(read_table(scope / "metrics/clone_members.tsv"), ("library_id", "sequence_id"))
        wanted = {key: value for key, value in expected_members.items() if key[0] in included}
        # Compare equivalence classes within each donor; numeric clone labels
        # are incidental and can change when a scope omits a library.
        check_clone_partitions(members, wanted, name)
        totals["reviewed_clone_members"] += len(members)
    return dict(totals)
