"""Reconstruct record qc results independently from retained scientific inputs."""

from pathlib import Path
from test.support.reference_tables import check_fields, number, read_table, index_rows
from collections import Counter, defaultdict


# The grouping_status tokens each qc.tsv reason fraction sums, besides any
# v_call_ or j_call_ status.
REASON_TOKENS = {"ineligible_junction": "junction_ineligible", "c_locus_conflict": "c_call_ineligible",
                 "c_call_incomplete_candidates": "c_call_ineligible",
                 "ineligible_backbone_representative": "contig_representative_ineligible"}
REASONS = ("v_call_unresolved", "j_call_unresolved", "junction_ineligible", "c_call_ineligible",
           "contig_representative_ineligible")


def grouping_reasons(status: str) -> set[str]:
    if status == "eligible":
        return set()
    return {"v_call_unresolved" if token.startswith("v_call_") else
            "j_call_unresolved" if token.startswith("j_call_") else REASON_TOKENS[token]
            for token in status.split(";")}


def verify_qc(root: Path, scope: Path) -> int:
    """Compare audit, trail and qc.tsv grouping populations with the two scientific mappings."""
    manifest = read_table(scope / "manifest/scope_manifest.tsv")
    table = read_table(scope / "metrics/qc_audit.tsv")
    stages = {"grouping_eligible", "assembly_backbone_collapse", "exact_sequence_collapse", "inference_representatives"}
    tail = index_rows([row for row in table if row["stage"] in stages],
                      ("library_id", "locus", "stage"))
    trail = index_rows([row for row in read_table(scope / "metrics/qc_trail.tsv")
                        if row["stage"] in {"grouping_eligible", "contig_collapse", "exact_sequence_collapse"}],
                       ("library_id", "locus", "stage"))
    qc = index_rows(read_table(scope / "metrics/qc.tsv"), ("library_id", "locus"))
    expected_keys = set()
    for context in manifest:
        library = context["library_id"]
        directory = root / "libraries" / library
        records = index_rows(read_table(directory / "airr" / f"{library}_airr_qc_pass.tsv"), ("sequence_id",))
        members = index_rows([row for branch in ("igh", "nonigh")
                              for row in read_table(directory / "metrics" / f"{library}_{branch}_grouping.tsv")],
                             ("sequence_id",))
        assert set(records) == set(members)
        transferred = Counter()
        distinct_sequences = set()
        for key, mapping in members.items():
            if mapping["grouping_eligible"] != "TRUE":
                continue
            backbone = (mapping["backbone_representative_sequence_id"],)
            final = (mapping["inference_representative_sequence_id"],)
            assert backbone in records and final in records
            assert members[backbone]["is_backbone_representative"] == "TRUE"
            assert members[final]["is_inference_representative"] == "TRUE"
            assert members[backbone]["inference_representative_sequence_id"] == final[0]
            assert records[backbone]["sequence"] == records[final]["sequence"], (key, "exact-sequence link")
            transferred[final] += int(records[key]["consensus_count"])
            if mapping["is_inference_representative"] == "TRUE":
                sequence_key = (records[key]["locus"], records[key]["sequence"])
                assert sequence_key not in distinct_sequences, (library, "duplicate inference sequence")
                distinct_sequences.add(sequence_key)
        for key, mapping in members.items():
            assert int(mapping["inference_consensus_count"]) == transferred[key], (library, key, "transferred support")
        shard = index_rows([row for branch in ("igh", "nonigh")
                            for row in read_table(directory / "qc" / f"{library}_qc_audit_{branch}_grouping.tsv")
                            if row["stage"] in stages], ("library_id", "locus", "stage"))
        for locus in ("IGH", "IGK", "IGL", "TRA", "TRB", "TRD", "TRG", "all"):
            selected = [key for key, row in records.items() if locus == "all" or row["locus"] == locus]
            eligible = [key for key in selected if members[key]["grouping_eligible"] == "TRUE"]
            backbones = [key for key in selected if members[key]["is_backbone_representative"] == "TRUE"]
            representatives = [key for key in selected if members[key]["is_inference_representative"] == "TRUE"]
            assert set(representatives) <= set(backbones) <= set(eligible)
            incoming = sum(int(records[key]["consensus_count"]) for key in selected)
            retained = sum(int(records[key]["consensus_count"]) for key in eligible)
            collapsed = sum(int(members[key]["inference_consensus_count"]) for key in representatives)
            assert retained == collapsed
            expected = {
                "grouping_eligible": ("eligibility", 30, "rearrangement", "rearrangement", len(selected), len(eligible), incoming, retained),
                "exact_sequence_collapse": ("grouping", 32, "backbone_representative", "inference_representative", len(backbones), len(representatives), retained, collapsed),
                "inference_representatives": ("checkpoint", 33, "inference_representative", "inference_representative", len(representatives), len(representatives), collapsed, collapsed),
            }
            # The released trail renames the units and numbers the stages from 1.
            released = {
                "grouping_eligible": ("eligibility", 15, "airr_record", "airr_record", *expected["grouping_eligible"][4:]),
                "exact_sequence_collapse": ("grouping", 17, "contig", "member", *expected["exact_sequence_collapse"][4:]),
            }
            if locus in ("IGH", "IGK", "IGL"):
                expected["assembly_backbone_collapse"] = ("grouping", 31, "rearrangement", "backbone_representative", len(eligible), len(backbones), retained, retained)
                released["contig_collapse"] = ("grouping", 16, "airr_record", "contig", *expected["assembly_backbone_collapse"][4:])
            for stage, values in expected.items():
                key = (library, locus, stage)
                expected_keys.add(key)
                fields = dict(zip(("kind", "stage_index", "input_unit", "output_unit", "n_units_in", "n_units_out", "consensus_count_in", "consensus_count_out"), values))
                check_fields(tail[key], fields, (scope.name, key))
                if locus != "all":
                    check_fields(shard[key], fields, key)
            for stage, values in released.items():
                fields = dict(zip(("kind", "stage_order", "input_unit", "output_unit", "n_units_in", "n_units_out", "consensus_count_in", "consensus_count_out"), values))
                check_fields(trail[(library, locus, stage)], fields, (scope.name, library, locus, stage))
            if locus != "all":
                # Each reason sums the original support of the records whose
                # grouping status names it, over common-QC support.
                reasons = Counter()
                for key in selected:
                    for reason in grouping_reasons(members[key]["grouping_status"]):
                        reasons[reason] += int(records[key]["consensus_count"])
                check_fields(qc[(library, locus)], {
                    "common_qc_consensus_count": incoming, "grouped_consensus_count": retained,
                    "grouping_eligible_fraction_of_common_qc": retained / incoming if incoming else None,
                    **{reason + "_fraction": reasons[reason] / incoming if incoming else None for reason in REASONS},
                }, (scope.name, library, locus))
        assert set(shard) == {key for key in expected_keys if key[0] == library and key[1] != "all"}
        for locus in ("IGK", "IGL", "TRA", "TRB", "TRD", "TRG"):
            groups = read_table(directory / "clonotypes" / f"{library}_{locus}_clonotypes.tsv")
            contributors = defaultdict(set)
            representatives = defaultdict(set)
            group_support = Counter()
            for key, mapping in members.items():
                if records[key]["locus"] != locus or mapping["grouping_eligible"] != "TRUE":
                    continue
                target = records[(mapping["inference_representative_sequence_id"],)]
                group = (target["junction"].strip().upper(), target["v_gene"], target["j_gene"])
                contributors[group].add(mapping["assembly_backbone_id"])
                representatives[group].add(mapping["inference_representative_sequence_id"])
                group_support[group] += int(records[key]["consensus_count"])
            observed = index_rows(groups, ("junction", "v_gene_key", "j_gene_key"))
            assert observed.keys() == contributors.keys(), (library, locus, "exact group population")
            for key, group in observed.items():
                check_fields(group, {"n_contigs": len(contributors[key]),
                                     "n_members": len(representatives[key]),
                                     "consensus_count": group_support[key]}, (library, locus, key))
    assert set(tail) == expected_keys
    assert set(qc) == {(context["library_id"], locus) for context in manifest
                       for locus in ("IGH", "IGK", "IGL", "TRA", "TRB", "TRD", "TRG")}
    assert not any("relax" in row["stage"] or row["stage"] == "informative_alignment_ge_200" for row in table)
    return len(tail)
