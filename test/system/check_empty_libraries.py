#!/usr/bin/env python3
"""Check publication and scientific availability for empty and all-failing libraries."""

import argparse
import csv
import gzip
from pathlib import Path
import re
import shutil
import sys

import yaml

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from test.support.reference_tables import check_fields, read_table
from test.components.record_qc.integration.reference_records import check_audit, expected_early_rows
from test.data_builders.empty_libraries import CONTROLS, CONSTRUCTED


def check_case(work: Path) -> None:
    config = yaml.safe_load((work / "config.yaml").read_text())
    root = Path(config["outs"]["results_dir"])
    stop, = [row for row in read_table(Path("test/data/expected/expected_records.tsv"))
             if row["transcript_junction"] == "IGH_stop"]
    for library in CONSTRUCTED:
        directory = root / "libraries" / library
        assert read_table(directory / "airr" / f"{library}_airr_qc_pass.tsv") == []
        annotated = read_table(directory / "airr" / f"{library}_igblast.tsv")
        assert len(annotated) == (1 if library == "FAIL" else 0)
        with gzip.open(directory / "export" / f"{library}_rearrangements.tsv.gz", "rt") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            assert "sequence_id" in reader.fieldnames and list(reader) == []
        with (directory / "trust4" / f"{library}_cdr3.out").open() as handle:
            cdr3 = list(csv.reader(handle, delimiter="\t"))
        expected = expected_early_rows([stop] if library == "FAIL" else [], cdr3, .1)
        actual = [row for shard in ("reannotate", "qc") for row in read_table(
            directory / "qc" / f"{library}_qc_audit_{shard}.tsv")]
        check_audit(actual, expected, library)
    for name in ("cohort", "empty_only"):
        release = root / "scopes" / name / "release"
        assert (release / ".complete.json").is_file()
        metrics = release / "metrics"
        headline = {row["library_id"]: row for row in read_table(metrics / "headline.tsv")}
        assert headline.keys() == set(CONSTRUCTED) | (set(CONTROLS) if name == "cohort" else set())
        for library in CONSTRUCTED:
            check_fields(headline[library], {"input_read_pairs": 1000,
                "igh_grouped_consensus_count": 0, "trb_grouped_consensus_count": 0,
                "igh_richness": 0, "trb_richness": 0,
                "igh_rarefied_inverse_simpson": None, "trb_rarefied_inverse_simpson": None,
                "igh_rarefaction_status": "empty_repertoire", "trb_rarefaction_status": "empty_repertoire",
                "igh_yield_per_million_read_pairs": 174000 if library == "FAIL" else 0,
                "class_callable_consensus_count": 0, "igh_mean_clone_median_shm_whole_v": None,
                "log2_igg_iga": None}, (name, library))
        # The working audit stays in the scope directory; the release has its trail.
        audit = read_table(root / "scopes" / name / "metrics" / "qc_audit.tsv")
        history = [row for row in audit if row["stage"] == "assembly_history" and row["library_id"] in CONSTRUCTED]
        assert len(history) == 2 and all(row["decision_status"] == "unknown_history" for row in history)
        qc = {(row["library_id"], row["locus"]): row for row in read_table(metrics / "qc.tsv")}
        for library in CONSTRUCTED:
            for locus in ("IGH", "TRB"):
                check_fields(qc[(library, locus)], {"common_qc_consensus_count": 0,
                    "grouped_consensus_count": 0, "grouping_eligible_fraction_of_common_qc": None,
                    "v_call_unresolved_fraction": None, "class_callable_fraction": None}, (name, library, locus))
        # FAIL's one annotated IGH record fails the productive criterion.
        check_fields(qc[("FAIL", "IGH")], {"igblast_consensus_count": 174, "nonproductive_fraction": 1,
                                            "common_qc_fraction_of_annotated": 0}, (name, "FAIL"))
        check_fields(qc[("EMPTY", "IGH")], {"igblast_consensus_count": 0, "nonproductive_fraction": None,
                                             "common_qc_fraction_of_annotated": None}, (name, "EMPTY"))
        for row in read_table(metrics / "selection_pressure.tsv"):
            if row["library_id"] in CONSTRUCTED:
                check_fields(row, {"n_clone_consensus_sequences": 0, "n_clone_consensus_sequences_used_cdr": 0, "n_clone_consensus_sequences_used_fwr": 0,
                    "n_clone_consensus_sequences_with_positions_cdr": 0, "n_clone_consensus_sequences_with_positions_fwr": 0,
                    "sigma_cdr": None, "sigma_fwr": None, "sigma_cdr_status": "no_eligible_clones",
                    "sigma_fwr_status": "no_eligible_clones"}, (name, row["library_id"]))
        for row in read_table(metrics / "shm.tsv"):
            if row["library_id"] in CONSTRUCTED:
                check_fields(row, {"n_members_eligible": 0, "n_clones_eligible": 0, "mean_clone_median_shm": None}, name)
        pairs = read_table(metrics / "donor_overlap_pairs.tsv")
        empty_pair = [row for row in pairs if {row["library_id_1"], row["library_id_2"]} == set(CONSTRUCTED)]
        assert len(empty_pair) == 2
        for row in empty_pair:
            check_fields(row, {"pair_status": "empty_both", "grouped_consensus_count_1": 0,
                "grouped_consensus_count_2": 0, "morisita_horn": None, "jaccard": None}, name)
        if name == "cohort":
            for library in CONTROLS:
                for locus in ("igh", "trb"):
                    check_fields(headline[library], {locus + "_rarefaction_status": "below_minimum_support",
                        locus + "_top_5_support_fraction": None}, library)
            control_pairs = [row for row in pairs if {row["library_id_1"], row["library_id_2"]} == set(CONTROLS)]
            assert len(control_pairs) == 2
            # Both endpoints are nonempty but below the minima, so the pair keeps
            # its values with a flag, and their diversity rows withhold the
            # inverse Simpson the report's Simpson difference reads.
            diversity = {(row["library_id"], row["locus"]): row for row in read_table(metrics / "diversity.tsv")}
            for row in control_pairs:
                check_fields(row, {"pair_status": "below_minimum_support"}, name)
                assert row["morisita_horn"] != "" and row["bray_curtis_proportion"] != ""
                for library in (row["library_id_1"], row["library_id_2"]):
                    check_fields(diversity[(library, row["locus"])], {"diversity_status": "below_minimum_support",
                        "inverse_simpson": None}, (library, row["locus"]))
            assert read_table(metrics / "rarefaction_curves.tsv"), "Minima must preserve curves"
    print("Constructed empty/all-failing libraries, empty scope and positive minima passed")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--work-dir", required=True, type=Path)
    args = parser.parse_args()
    check_case(args.work_dir.resolve())
