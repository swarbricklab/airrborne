"""Independent reference for the released diversity and support-depth values.

This module reconstructs every value in diversity.tsv and
rarefaction_curves.tsv from the released
support vector in group_abundance.tsv, then compares the reconstruction
with what the pipeline published. It imports nothing from the pipeline, so the
values it computes do not depend on workflow/lib/diversity.R. Keep it that way:
everything below uses the standard library only.

Each count is recovered read support for one clone or single-chain clonotype.
It is not a cell or molecule count. Drawing these support units without
replacement describes the observed vector; it neither estimates unseen clones
nor replaces reassembling downsampled FASTQs.

At depth m, a group's sampled count X is hypergeometric. The three curves are
E[number of observed groups], exp(E[Shannon entropy]) and the reciprocal of
expected Simpson concentration. This is deterministic; there are no random
replicates. Hill-number interpretation: doi:10.2307/1934352. Finite-sample
richness: doi:10.2307/1934145. The reconstruction uses a probability recurrence
rather than the pipeline's R hypergeometric functions, and equal original counts
share one calculation.
"""

from bisect import bisect_right
from collections import Counter, defaultdict
from functools import lru_cache
import math
from pathlib import Path
from test.support.reference_tables import check_fields, number, read_table

MAX_CURVE_POINTS = 42  # Forty grid points plus exact target and endpoint.
OBSERVED_METRICS = (
    "shannon_entropy", "hill_shannon", "inverse_simpson", "clonality",
    "top_1_support_fraction", "top_5_support_fraction", "top_10_support_fraction",
)


def group_expectations(count: int, total: int, depth: int) -> tuple[float, float]:
    """Presence probability and expected entropy from a full hypergeometric PMF.

    Start at a mode with relative weight one and recur in each direction. The
    normalisation cancels its unknown probability. Weights never require large
    binomial coefficients or subtraction of nearly equal log-factorials.
    """
    lower = max(0, depth - (total - count))
    upper = min(count, depth)
    mode = min(upper, max(lower, (depth + 1) * (count + 1) // (total + 2)))
    weights = [0.0] * (upper - lower + 1)
    weights[mode - lower] = 1.0
    weight = 1.0
    for drawn in range(mode, lower, -1):
        weight *= (drawn * (total - count - depth + drawn)
                   / ((count - drawn + 1) * (depth - drawn + 1)))
        weights[drawn - 1 - lower] = weight
    weight = 1.0
    for drawn in range(mode, upper):
        weight *= ((count - drawn) * (depth - drawn)
                   / ((drawn + 1) * (total - count - depth + drawn + 1)))
        weights[drawn + 1 - lower] = weight
    normaliser = math.fsum(weights)
    positive = range(max(1, lower), upper + 1)
    presence = math.fsum(weights[x - lower] for x in positive) / normaliser
    entropy = math.fsum(-weights[x - lower] * (x / depth) * math.log(x / depth)
                        for x in positive) / normaliser
    return presence, entropy


@lru_cache(maxsize=4096)
def interpolate(counts: tuple[int, ...], depth: int) -> dict[str, float | None]:
    """Expected q0, exp(expected Shannon), and reciprocal expected concentration."""
    total = sum(counts)
    if not 0 <= depth <= total:
        raise ValueError("depth must be between zero and total support")
    if depth == 0:
        return {"richness": 0, "hill_shannon": None, "inverse_simpson": None}
    if depth == total:
        p = [count / total for count in counts]
        return {"richness": len(counts),
                "hill_shannon": math.exp(-math.fsum(x * math.log(x) for x in p)),
                "inverse_simpson": 1 / math.fsum(x * x for x in p)}
    terms = [(multiplicity, group_expectations(count, total, depth))
             for count, multiplicity in Counter(counts).items()]
    concentration = sum(count * (count - 1) for count in counts) / (total * (total - 1))
    return {"richness": math.fsum(n * value[0] for n, value in terms),
            "hill_shannon": math.exp(math.fsum(n * value[1] for n, value in terms)),
            "inverse_simpson": depth / (1 + (depth - 1) * concentration)}


def expected_detail(counts: tuple[int, ...], minimum: int, target: int | None) -> dict:
    total, richness = sum(counts), len(counts)
    status = ("empty_repertoire" if total == 0 else
              "below_minimum_support" if total < minimum else "available")
    rarefaction_status = ("not_configured" if target is None else
                         status if status != "available" else
                         "below_target_support" if total < target else "available")
    result = {"grouped_consensus_count": total, "richness": richness,
              "diversity_status": status, "rarefaction_status": rarefaction_status,
              **dict.fromkeys(OBSERVED_METRICS)}
    if status == "available":
        p = sorted((count / total for count in counts), reverse=True)
        entropy = -math.fsum(x * math.log(x) for x in p)
        result.update(shannon_entropy=entropy, hill_shannon=math.exp(entropy),
                      inverse_simpson=1 / math.fsum(x * x for x in p),
                      clonality=min(1, max(0, 1 - entropy / math.log(richness)))
                      if richness >= 2 else None)
        result.update({f"top_{n}_support_fraction": math.fsum(p[:n]) for n in (1, 5, 10)})
    target_values = (interpolate(counts, target) if rarefaction_status == "available"
                     else dict.fromkeys(("richness", "hill_shannon", "inverse_simpson")))
    result.update({f"rarefied_{key}": value for key, value in target_values.items()})
    concentration = (sum(count * (count - 1) for count in counts) / (total * (total - 1))
                     if total >= 2 else None)
    asymptote_status = ("insufficient_support" if total < 2 else
                        "no_repeated_support" if concentration == 0 else "available")
    result.update(asymptote_status=asymptote_status,
                  inverse_simpson_asymptote=1 / concentration if asymptote_status == "available" else None)
    return result


def verify_ranks(vector: list[tuple[int, dict]], junctions: dict, context: object) -> None:
    """Check group_rank and the row order of one library and locus.

    A group's rank is one plus the number of groups with strictly more support,
    so tied groups share the smaller rank. Rows run in rank order; tied
    clonotypes follow clonotype_key, and tied IGH clones follow their smallest
    member junction in the library, then clone_id.
    """
    counts = sorted(count for count, _ in vector)
    for count, row in vector:
        check_fields(row, {"group_rank": 1 + len(counts) - bisect_right(counts, count)}, context)
    order = [(int(row["group_rank"]), junctions.get(row["group_key"], ""), row["group_key"])
             for _, row in vector]
    if order != sorted(order):
        raise AssertionError(f"support-vector rows are not in rank order: {context}")


def verify_scope(directory: Path) -> tuple[int, int]:
    """Check every released curve point, endpoint, target, total, status and dominance.

    Returns the number of diversity rows and the number of curve points checked.
    Diagnostic curves and the q2 asymptote remain available below a
    configured diversity minimum, so this checks them there as well. Each
    support vector's n_members must equal the member rows of its clone or
    clonotype in that library, and its group_rank and row order must follow
    its support.
    """
    metrics = Path(directory) / "metrics"
    clone_members = read_table(metrics / "clone_members.tsv")
    members = Counter((row["library_id"], "IGH", row["clone_id"]) for row in clone_members)
    clonotype_keys = {row["clonotype_id"]: row["clonotype_key"]
                      for row in read_table(metrics / "clonotypes.tsv")}
    members.update((row["library_id"], row["locus"], clonotype_keys[row["clonotype_id"]])
                   for row in read_table(metrics / "clonotype_members.tsv"))
    junctions = defaultdict(dict)
    for row in clone_members:
        smallest = junctions[row["library_id"]].get(row["clone_id"])
        if smallest is None or row["junction"] < smallest:
            junctions[row["library_id"]][row["clone_id"]] = row["junction"]
    vectors = defaultdict(list)
    identities = set()
    for row in read_table(metrics / "group_abundance.tsv"):
        key = (row["library_id"], row["locus"], row["group_unit"])
        identity = (*key, row["group_key"])
        if identity in identities:
            raise AssertionError(f"duplicate support-vector key: {identity}")
        identities.add(identity)
        check_fields(row, {"n_members": members[(row["library_id"], row["locus"], row["group_key"])]},
                     identity)
        count = number(row["consensus_count"])
        if count is None or count <= 0 or not count.is_integer():
            raise AssertionError(f"invalid support count at {identity}")
        vectors[key].append((int(count), row))
    curves = defaultdict(list)
    for row in read_table(metrics / "rarefaction_curves.tsv"):
        curves[(row["library_id"], row["locus"], row["group_unit"])].append(row)
    details = read_table(metrics / "diversity.tsv")
    seen = set()
    points = 0
    for row in details:
        key = (row["library_id"], row["locus"], row["group_unit"])
        if key in seen:
            raise AssertionError(f"duplicate diversity row: {key}")
        seen.add(key)
        counts = tuple(sorted(count for count, _ in vectors[key]))
        total = sum(counts)
        for count, group in vectors[key]:
            check_fields(group, {"support_fraction": count / total}, key)
        verify_ranks(vectors[key], junctions[key[0]] if key[1] == "IGH" else {}, key)
        minimum = int(row["min_consensus_count_for_diversity"])
        target_value = number(row["rarefaction_target"])
        target = int(target_value) if target_value is not None else None
        check_fields(row, expected_detail(counts, minimum, target), key)
        curve = curves[key]
        depths = [int(point["support_depth"]) for point in curve]
        if target is None:
            if curve:
                raise AssertionError(f"unconfigured curve: {key}")
            continue
        if depths != sorted(set(depths)) or total not in depths or len(depths) > MAX_CURVE_POINTS:
            raise AssertionError(f"curve grid or endpoint is invalid: {key}")
        if target <= total and target not in depths:
            raise AssertionError(f"target is absent from curve: {key}")
        for point, depth in zip(curve, depths):
            check_fields(point, {**interpolate(counts, depth), "is_target": depth == target,
                                 "is_endpoint": depth == total}, (*key, depth))
        points += len(curve)
    if set(vectors) - seen or set(curves) - seen:
        raise AssertionError("support vectors or curves have no diversity row")
    return len(seen), points
