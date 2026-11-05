"""Reconstruct igh clones results independently from retained scientific inputs."""

from pathlib import Path
from test.support.reference_tables import check_fields, number, read_table, index_rows


FIT_STATUSES = {"fitted", "too_few_distances", "constant_distances", "fit_failed", "out_of_range"}
FIT_PARAMETERS = ("near_weight", "near_shape", "near_scale", "far_weight", "far_mean", "far_sd",
                  "absolute_log_likelihood", "fit_sensitivity")


def verify_calibration(root: Path, scope: Path, min_fit_distances: int,
                       fit_subsample: int | None) -> tuple[int, int]:
    """Check calibration status without treating an unusable fit as exclusion.

    Each donor's n_members and n_libraries must count the inference
    representatives its libraries' grouping retained, and a donor below
    `min_fit_distances` distances must report too_few_distances. An attempted fit
    must use the distances `fit_subsample` allows, and only a fitted donor holds
    mixture parameters.
    """
    manifest = read_table(scope / "manifest/scope_manifest.tsv")
    donors = {row["donor_id"] for row in manifest}
    fits = index_rows(read_table(scope / "clones" / f"{scope.name}_donor_clone_thresholds.tsv"), ("donor_id",))
    assert set(fits) == {(donor,) for donor in donors}
    diagnostics = read_table(scope / "clones" / f"{scope.name}_clone_distance_diagnostics.tsv")
    assert {row["donor_id"] for row in diagnostics} == donors
    threshold, = read_table(scope / "clones" / f"{scope.name}_clone_threshold.tsv")
    assert int(threshold["n_donors"]) == len(fits)
    assert int(threshold["n_donors_fitted"]) == sum(row["fit_status"] == "fitted" for row in fits.values())
    membership = read_table(scope / "metrics/clone_members.tsv")
    audit = index_rows([row for row in read_table(scope / "metrics/qc_audit.tsv")
                       if row["stage"] == "donor_clones"], ("donor_id",))
    failed_with_members = 0
    for donor in sorted(donors):
        fit = fits[(donor,)]
        assert fit["scope_id"] == scope.name
        assert fit["fit_status"] in FIT_STATUSES, (scope.name, donor, fit["fit_status"])
        assert (number(fit["fitted_threshold"]) is not None) == (fit["fit_status"] == "fitted")
        assert (int(fit["n_distances"]) < min_fit_distances) == (fit["fit_status"] == "too_few_distances")
        n_distances = int(fit["n_distances"])
        attempted = fit["fit_status"] not in ("too_few_distances", "constant_distances")
        expected_fitted = None if not attempted else n_distances if fit_subsample is None else min(n_distances, fit_subsample)
        assert number(fit["n_distances_fitted"]) == expected_fitted, (scope.name, donor)
        for column in FIT_PARAMETERS:
            assert (number(fit[column]) is not None) == (fit["fit_status"] == "fitted"), (scope.name, donor, column)
        assert sum(int(row["n_distances_in_bin"]) for row in diagnostics if row["donor_id"] == donor) == int(fit["n_distances"])
        expected = {}
        for context in manifest:
            if context["donor_id"] != donor:
                continue
            library = context["library_id"]
            mapping = read_table(root / "libraries" / library / "metrics" / f"{library}_igh_grouping.tsv")
            expected.update({(library, library + "~" + row["sequence_id"]): int(row["inference_consensus_count"])
                             for row in mapping if row["is_inference_representative"] == "TRUE"})
        assigned = [row for row in membership if row["donor_id"] == donor]
        actual = {(row["library_id"], row["sequence_id"]): int(row["member_consensus_count"])
                  for row in assigned}
        assert len(assigned) == len(expected) and actual == expected, (scope.name, donor)
        check_fields(fit, {"n_members": len(expected),
                           "n_libraries": len({library for library, _ in expected})}, (scope.name, donor))
        clones = len({row["clone_id"] for row in assigned})
        check_fields(audit[(donor,)], {"n_units_in": clones, "n_units_out": clones,
                                      "consensus_count_in": sum(expected.values()),
                                      "consensus_count_out": sum(expected.values())}, (scope.name, donor))
        if fit["fit_status"] != "fitted" and expected:
            assert clones > 0
            failed_with_members += 1
    return len(fits), failed_with_members


def verify_clones(scope: Path) -> int:
    """Rebuild each released clone's descriptors and counts from its members.

    The V and J entities and the junction length are the one value every member
    of the clone carries, `n_libraries` counts the members' distinct libraries
    and `consensus_count` sums their transferred support.
    """
    members: dict[str, list[dict[str, str]]] = {}
    for row in read_table(scope / "metrics/clone_members.tsv"):
        members.setdefault(row["clone_id"], []).append(row)
    clones = index_rows(read_table(scope / "metrics/clones.tsv"), ("clone_id",))
    assert {key[0] for key in clones} == set(members), scope.name
    for (clone_id,), clone in clones.items():
        rows = members[clone_id]
        for column in ("donor_id", "v_gene_key", "v_gene_names", "j_gene_key", "j_gene_names",
                       "junction_length"):
            values = {row[column] for row in rows}
            assert values == {clone[column]}, (scope.name, clone_id, column, values)
        check_fields(clone, {"n_members": len(rows),
                             "n_libraries": len({row["library_id"] for row in rows}),
                             "consensus_count": sum(int(row["member_consensus_count"]) for row in rows)},
                     (scope.name, clone_id))
    return len(clones)
