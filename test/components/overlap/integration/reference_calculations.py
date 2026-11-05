"""Reconstruct overlap results independently from retained scientific inputs."""

from pathlib import Path
from test.support.reference_tables import check_fields, number, read_table, index_rows
from collections import Counter, defaultdict
from fractions import Fraction
from itertools import combinations
import math


def pairwise(first: dict[str, int], second: dict[str, int], minimum: int = 0) -> dict:
    """Every released pair value for two support vectors keyed by group."""
    keys = first.keys() | second.keys()
    shared = first.keys() & second.keys()
    total_a, total_b = sum(first.values()), sum(second.values())
    shared_a = sum(first[key] for key in shared)
    shared_b = sum(second[key] for key in shared)
    status = ("empty_both" if not total_a and not total_b else
              "empty_first" if not total_a else "empty_second" if not total_b else "available")
    morisita = None
    if status == "available":
        numerator = 2 * math.fsum(first.get(key, 0) / total_a * second.get(key, 0) / total_b for key in keys)
        denominator = math.fsum((value / total_a) ** 2 for value in first.values())
        denominator += math.fsum((value / total_b) ** 2 for value in second.values())
        morisita = numerator / denominator
    pair_status = status if status != "available" else (
        "below_minimum_support" if min(total_a, total_b) < minimum else "available")
    result = {"n_groups_1": len(first), "n_groups_2": len(second),
              "grouped_consensus_count_1": total_a,
              "grouped_consensus_count_2": total_b,
              "n_groups_shared": len(shared),
              "shared_consensus_count_1": shared_a, "shared_consensus_count_2": shared_b,
              "shared_support_fraction_1": shared_a / total_a if total_a else None,
              "shared_support_fraction_2": shared_b / total_b if total_b else None,
              "jaccard": len(shared) / len(keys) if status == "available" else None,
              "morisita_horn": morisita, "min_consensus_count_per_library": minimum,
              "pair_status": pair_status, "bray_curtis_proportion": None, "n_shared_higher_in_2": None,
              "n_shared_higher_in_1": None, "n_shared_equal": None}
    if status == "available":
        p = {key: Fraction(first.get(key, 0), total_a) for key in keys}
        q = {key: Fraction(second.get(key, 0), total_b) for key in keys}
        differences = {key: q[key] - p[key] for key in keys}
        result["bray_curtis_proportion"] = float(sum(abs(value) for value in differences.values()) / 2)
        for side, predicate in (("higher_in_2", lambda x: x > 0), ("higher_in_1", lambda x: x < 0),
                                ("equal", lambda x: x == 0)):
            result[f"n_shared_{side}"] = sum(predicate(differences[key]) for key in shared)
    return result


def verify_sharing(scope: Path, params: dict) -> tuple[int, int]:
    """Rebuild the overlap pairs and the persistent clones from the support vector.

    A separate collector writes shared_groups.tsv, so recovering the same
    clones and detection counts from group_abundance.tsv checks that the
    two released tables describe one repertoire. A clone belongs to one donor, so
    a donor and a clone key together name it. Shared-clone rows must run by
    donor, locus, the clone's total support (largest first), group_key and
    library_id. Each donor and locus summary must carry the pair counts and the
    occupancy of these support vectors.
    """
    vectors = defaultdict(dict)
    occupancy = defaultdict(set)
    libraries = defaultdict(set)
    tissues = {}
    for row in read_table(scope / "manifest/scope_manifest.tsv"):
        libraries[row["donor_id"]].add(row["library_id"])
        tissues[row["library_id"]] = row["tissue_id"]
    for row in read_table(scope / "metrics/group_abundance.tsv"):
        vectors[(row["donor_id"], row["locus"], row["library_id"])][row["group_key"]] = int(row["consensus_count"])
        occupancy[(row["donor_id"], row["locus"], row["group_key"])].add(row["library_id"])
    pairs = read_table(scope / "metrics/donor_overlap_pairs.tsv")
    pair_keys = index_rows(pairs, ("donor_id", "locus", "library_id_1", "library_id_2"))
    assert set(pair_keys) == {
        (donor, locus, first, second)
        for donor, members in libraries.items()
        for locus in params["overlap"]["loci"]
        for first, second in combinations(
            sorted(members, key=lambda library: (tissues[library], library)), 2)
    }, "Overlap rows differ from the complete within-donor comparison set"
    for row in pairs:
        first = vectors[(row["donor_id"], row["locus"], row["library_id_1"])]
        second = vectors[(row["donor_id"], row["locus"], row["library_id_2"])]
        check_fields(row, pairwise(first, second,
                                  params["overlap"]["min_consensus_count"].get(row["locus"], 0)),
                     (row["donor_id"], row["locus"], row["library_id_1"], row["library_id_2"]))
        check_fields(row, {"group_unit": "clone" if row["locus"] == "IGH" else "clonotype"}, row["locus"])
    summaries = index_rows(read_table(scope / "metrics/donor_overlap_summary.tsv"), ("donor_id", "locus"))
    assert set(summaries) == {(donor, locus) for donor, members in libraries.items() if len(members) >= 2
                              for locus in params["overlap"]["loci"]}, \
        "Donor overlap rows differ from the complete donor and locus set"
    for (donor, locus), row in summaries.items():
        members = libraries[donor]
        counts = Counter(key for library in members for key in vectors[(donor, locus, library)])
        check_fields(row, {"group_unit": "clone" if locus == "IGH" else "clonotype",
                           "n_libraries_in_donor": len(members),
                           "n_libraries_nonempty": sum(bool(vectors[(donor, locus, library)]) for library in members),
                           "n_pairs_total": len(members) * (len(members) - 1) // 2,
                           "n_pairs_available": sum(pair["donor_id"] == donor and pair["locus"] == locus
                                                    and pair["pair_status"] == "available" for pair in pairs),
                           "n_groups_union": len(counts),
                           "n_groups_shared_2plus": sum(count >= 2 for count in counts.values()),
                           "n_groups_shared_all": sum(count == len(members) for count in counts.values())},
                     (donor, locus))
    # The shared-clone table carries one row per library of the donor, including
    # the libraries that recovered nothing, so the listed clones are its distinct
    # keys rather than its row count.
    shared_rows = read_table(scope / "metrics/shared_groups.tsv")
    shared = index_rows(shared_rows, ("donor_id", "locus", "group_key", "library_id"))
    listed = {key for key, members in occupancy.items()
              if len(members) >= 2 and key[1] in params["overlap"]["loci"]}
    assert set(shared) == {(*key, library) for key in listed for library in libraries[key[0]]}, \
        "Shared-clone rows differ from the complete clone-by-library set"
    for row in shared.values():
        key = (row["donor_id"], row["locus"], row["group_key"])
        vector = vectors[(key[0], key[1], row["library_id"])]
        count, total = vector.get(key[2], 0), sum(vector.values())
        check_fields(row, {"n_libraries_detected": len(occupancy[key]),
                           "n_libraries_in_donor": len(libraries[key[0]]),
                           "consensus_count": count,
                           "support_fraction": count / total if total else None}, key)
    totals = defaultdict(int)
    for key in listed:
        totals[key] = sum(vectors[(key[0], key[1], library)].get(key[2], 0) for library in libraries[key[0]])
    order = [(row["donor_id"], row["locus"], -totals[(row["donor_id"], row["locus"], row["group_key"])],
              row["group_key"], row["library_id"]) for row in shared_rows]
    assert order == sorted(order), "Shared-clone rows are not in clone-support order"
    return len(pairs), len(listed)
