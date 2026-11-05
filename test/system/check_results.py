#!/usr/bin/env python3
"""Check both system modes with reviewed records and independent calculations.

SHM reconstruction reads measured features; it does not independently align
sequences. Pinned component tests protect alignment and numerical estimators.
"""

import argparse
import csv
from collections import Counter
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from test.support.reference_tables import check_fields, read_table, index_rows
from test.support.expectations import verify_expected_values
from test.components.diversity.integration import reference_calculations as DIVERSITY
from test.components.shm.integration.reference_calculations import verify_shm, REGIONS, SHM_GROUPS
from test.components.isotype.integration.reference_calculations import verify_isotypes
from test.components.overlap.integration.reference_calculations import verify_sharing
from test.components.inputs.integration.reference_calculations import verify_input_depth, verify_yield
from test.components.annotation.integration.reference_calculations import verify_reference_receipt
from test.components.record_qc.integration.reference_calculations import verify_qc
from test.components.record_qc.integration.reference_records import verify_records
from test.components.selection.integration.reference_calculations import verify_selection_estimates, verify_selection_audit
from test.components.igh_clones.integration.reference_calculations import verify_calibration, verify_clones

def expected_scope_libraries(config: dict) -> dict[str, set[str]]:
    """Read the fixture populations independently from its input CSV files."""
    with Path(config["deps"]["samplesheet"]).open(newline="") as handle:
        libraries = {row["library_id"] for row in csv.DictReader(handle)}
    scopes = {"cohort": libraries}
    if config.get("membership_csv"):
        with Path(config["membership_csv"]).open(newline="") as handle:
            for row in csv.DictReader(handle):
                scope = scopes.setdefault(row["project"], set())
                if row["usage"] != "fail_qc":
                    scope.add(row["library_id"])
    return scopes


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results", required=True, type=Path)
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--expected-values", type=Path,
                        help="committed expectations to check against this run")
    parser.add_argument("--expected-records", type=Path,
                        help="reviewed annotation, QC, grouping and clone outcomes")
    parser.add_argument("--scope", help="check only this configured scope")
    args = parser.parse_args()
    import yaml
    config = yaml.safe_load(args.config.read_text())
    populations = expected_scope_libraries(config)
    if args.scope:
        populations = {args.scope: populations[args.scope]}
    totals = Counter()
    checked_depth = set()
    totals.update(reference_databases=verify_reference_receipt(args.results))
    for name, libraries in sorted(populations.items()):
        scope = args.results / "scopes" / name
        manifest = index_rows(read_table(scope / "manifest/scope_manifest.tsv"), ("library_id",))
        assert {key[0] for key in manifest} == libraries, f"{name}: incomplete scope manifest"
        totals.update(input_depth_rows=verify_input_depth(args.results, libraries - checked_depth))
        checked_depth.update(libraries)
        for table in ("yield.tsv", "isotype_summary.tsv"):
            rows = index_rows(read_table(scope / "metrics" / table), ("library_id",))
            assert {key[0] for key in rows} == libraries, f"{name}: incomplete {table}"
        diversity = index_rows(read_table(scope / "metrics/diversity.tsv"),
                               ("library_id", "locus"))
        # Every library has a diversity row at all seven loci, with or without a target.
        diversity_loci = {"IGH", "IGK", "IGL", "TRA", "TRB", "TRD", "TRG"}
        assert set(diversity) == {(library, locus) for library in libraries
                                  for locus in diversity_loci}
        assert all(row["group_unit"] == ("clone" if row["locus"] == "IGH" else "clonotype")
                   for row in diversity.values())
        selection = index_rows(read_table(scope / "metrics/selection_pressure.tsv"),
                                ("library_id", "isotype_group"))
        assert set(selection) == {(library, group) for library in libraries
                                  for group in ("all_IGH", "IgM", "IgG", "IgA")}
        totals.update(selection_rows=verify_selection_estimates(list(selection.values())))
        totals.update(selection_audit_rows=verify_selection_audit(args.results, scope))
        shm = index_rows(read_table(scope / "metrics/shm.tsv"),
                          ("library_id", "isotype_group", "v_region"))
        assert set(shm) == {(library, group, region) for library in libraries
                            for group in SHM_GROUPS for region in ("whole_v", *REGIONS)}
        rows, points = DIVERSITY.verify_scope(scope)
        totals.update(diversity_rows=rows, curve_points=points)
        for row in read_table(scope / "metrics/diversity.tsv"):
            check_fields(row, {"min_consensus_count_for_diversity": config["params"]["diversity"]["min_consensus_count"].get(row["locus"], 0),
                               "rarefaction_target": config["params"]["rarefaction"]["targets"].get(row["locus"])}, row["library_id"])
        composition, ratios, class_only = verify_isotypes(args.results, scope)
        detail, shm = verify_shm(args.results, scope, config["params"]["shm"])
        pairs, shared_groups = verify_sharing(scope, config["params"])
        yields = verify_yield(args.results, scope)
        qc_rows = verify_qc(args.results, scope)
        donor_fits, failed_with_members = verify_calibration(
            args.results, scope, config["params"]["lineage"]["min_fit_distances"],
            config["params"]["lineage"]["fit_subsample"])
        clones = verify_clones(scope)
        totals.update(composition_rows=composition, isotype_summary_rows=ratios,
                      clone_shm_rows=detail, shm_rows=shm, overlap_pairs=pairs,
                      shared_groups=shared_groups,
                      yield_rows=yields, qc_grouping_rows=qc_rows,
                      donor_fits=donor_fits, unfitted_donors_assigned=failed_with_members,
                      clones=clones)
        print(f"{scope.name}: all metric reconstructions passed; class-only support={class_only}", flush=True)
    if args.expected_values:
        totals.update(expected_values=verify_expected_values(
            args.results, args.expected_values, args.scope,
            set().union(*populations.values())))
    if args.expected_records:
        totals.update(verify_records(args.results, args.expected_records, populations, config))
    print("Reconstructed: " + json.dumps(dict(totals), sort_keys=True))


if __name__ == "__main__":
    main()
