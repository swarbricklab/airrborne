"""Reconstruct shm results independently from retained scientific inputs."""

from pathlib import Path
from test.support.reference_tables import check_fields, number, read_table, index_rows
from collections import defaultdict
import math
import statistics
from test.components.isotype.integration.reference_calculations import resolve_call
REGIONS = ("fwr1", "cdr1", "fwr2", "cdr2", "fwr3")
SHM_GROUPS = ("all_IGH", "IgM", "IgG", "IgA", "IgG1", "IgG2", "IgG3", "IgG4", "IgA1", "IgA2")


def quantile(values: list[float], probability: float) -> float | None:
    if not values:
        return None
    values = sorted(values)
    position = (len(values) - 1) * probability
    lower = math.floor(position)
    upper = math.ceil(position)
    return values[lower] + (position - lower) * (values[upper] - values[lower])


def region_values(features: dict[str, str], region: str, settings: dict) -> tuple:
    regions = REGIONS if region == "whole_v" else (region,)
    rate = number(features["shm_frequency" if region == "whole_v" else f"shm_frequency_{region}"])
    occupied = [number(features[f"germline_occupied_{name}"]) for name in regions]
    supported = [number(features[f"informative_supported_{name}"]) for name in regions]
    covered = all(o is not None and s is not None and o > 0 and
                  s / o >= settings["min_supported_subregion_fraction"]
                  for s, o in zip(supported, occupied))
    eligible = rate is not None and covered
    if region == "whole_v":
        eligible = eligible and sum(occupied) >= settings["min_germline_occupied_v_positions"]
    mutations = [number(features[f"mutations_{name}"]) for name in regions]
    votes = [number(features[f"median_majority_votes_{name}"]) for name in regions]
    return (rate, eligible, sum(mutations) if all(x is not None for x in mutations) else None,
            statistics.mean(votes) if all(x is not None for x in votes) else None)


def supported_positions(features: dict[str, str], regions: tuple[str, ...]) -> float | None:
    """The source's vote-supported positions over `regions`, or None where a count is missing."""
    counts = [number(features[f"informative_supported_{name}"]) for name in regions]
    return sum(counts) if all(count is not None for count in counts) else None


def sequence_shm_evidence(identity: str, records: dict, mapping: dict,
                          features: dict, settings: dict) -> dict:
    """Reconstruct each region from one independently qualifying source contig.

    Each region gives the selected source's ID, its measurements and the supported
    positions of the source ranked first, which ineligible members also report.
    """
    canonical = features[(identity,)]
    candidates = [key for key, row in mapping.items()
                  if row["is_backbone_representative"] == "TRUE" and
                  row["inference_representative_sequence_id"] == identity]
    assert candidates, (identity, "missing retained backbone sources")
    sequence_fields = ["qualified_c_call", "shm_frequency", "informative_v",
                       *(f"{prefix}_{region}" for region in REGIONS
                         for prefix in ("shm_frequency", "mutations", "informative", "germline_occupied"))]
    for key in candidates:
        assert records[key]["sequence"] == records[(identity,)]["sequence"], (identity, key, "evidence DNA")
        assert all(features[key][column] == canonical[column] for column in sequence_fields), \
            (identity, key, "sequence-derived features disagree")
    result = {}
    for region in ("whole_v", *REGIONS):
        covered_regions = REGIONS if region == "whole_v" else (region,)
        def coverage(key: tuple) -> float:
            count = supported_positions(features[key], covered_regions)
            return math.inf if count is None else -count
        # A qualifying source ranks first. Coverage then wins before native record
        # support. The row's own support is independent of the support later
        # transferred from its siblings.
        chosen = sorted(candidates, key=lambda key: (
            not region_values(features[key], region, settings)[1], coverage(key),
            -int(records[key]["consensus_count"]), key[0]))[0]
        supported = supported_positions(features[chosen], covered_regions)
        if region_values(features[chosen], region, settings)[1]:
            result[region] = (chosen[0], region_values(features[chosen], region, settings), supported)
        else:
            rate, _, mutations, _ = region_values(canonical, region, settings)
            result[region] = ("", (rate, False, mutations, None), supported)
    return result


# Released clone_members.tsv names for three donor shard columns.
MEMBER_RELEASE_NAMES = {"v_gene_key": "v_gene", "j_gene_key": "j_gene",
                        "member_consensus_count": "consensus_count"}


def verify_shm(root: Path, scope: Path, settings: dict) -> tuple[int, int]:
    """Rebuild member SHM evidence and both SHM tables from the retained features.

    The donor member shards keep every measured rate and the evidence IDs, so
    the evidence check and the clone SHM rebuild read them. Each released member
    column must equal its shard column, under the release name for three. The
    released member SHM table must hold one row per member and region with a
    finite rate, its count, rate and votes empty where the region is not
    SHM-eligible. The released clone SHM table must hold exactly the rebuilt rows
    with an eligible member.
    """
    released = read_table(scope / "metrics/clone_members.tsv")
    members = [row for path in sorted((scope / "intermediate").glob("*/*_clone_members.tsv"))
               for row in read_table(path)]
    shards = index_rows(members, ("library_id", "sequence_id"))
    assert set(shards) == {(row["library_id"], row["sequence_id"]) for row in released}
    for row in released:
        shard = shards[(row["library_id"], row["sequence_id"])]
        for column, value in row.items():
            expected = shard[MEMBER_RELEASE_NAMES.get(column, column)]
            assert value == expected, (row["sequence_id"], column, value, expected)
    member_shm = index_rows(read_table(scope / "metrics/clone_member_shm.tsv"),
                            ("library_id", "sequence_id", "v_region"))
    expected_member_shm = set()
    features = {}
    records = {}
    mappings = {}
    for library in {row["library_id"] for row in members}:
        directory = root / "libraries" / library
        table = directory / "metrics" / f"{library}_igh_features.tsv"
        features[library] = index_rows(read_table(table), ("sequence_id",))
        records[library] = index_rows(read_table(directory / "airr" / f"{library}_airr_qc_pass.tsv"), ("sequence_id",))
        mappings[library] = index_rows(read_table(directory / "metrics" / f"{library}_igh_grouping.tsv"), ("sequence_id",))
    units = defaultdict(list)
    for row in members:
        library = row["library_id"]
        assert row["sequence_id"].startswith(library + "~")
        identity = row["sequence_id"][len(library) + 1:]
        feature = features[library][(identity,)]
        evidence = sequence_shm_evidence(identity, records[library], mappings[library], features[library], settings)
        assert row["qualified_c_call"] == feature["qualified_c_call"]
        major, minor = resolve_call(row["qualified_c_call"])
        groups = {"all_IGH"} | ({major, minor} & set(SHM_GROUPS))
        for group in groups:
            units[(library, row["clone_id"], group)].append((row, evidence))
        for region, (source_id, (rate, eligible, mutations, votes), supported) in evidence.items():
            wanted = {f"shm_eligible_{region}": eligible,
                      f"shm_evidence_sequence_id_{region}": library + "~" + source_id if source_id else ""}
            if region in ("whole_v", "cdr2", "fwr3"):
                wanted.update({f"shm_frequency_{region}": rate, f"mutations_{region}": mutations})
            check_fields(row, wanted, row["sequence_id"])
            if rate is None:
                continue
            key = (library, row["sequence_id"], region)
            expected_member_shm.add(key)
            assert key in member_shm, (key, "member SHM row with a measured rate is absent")
            whole_v = region == "whole_v"
            check_fields(member_shm[key], {
                "clone_id": row["clone_id"], "shm_eligible": eligible,
                "n_informative_positions": number(feature["informative_v" if whole_v else f"informative_{region}"]),
                "n_germline_occupied_positions": sum(number(feature[f"germline_occupied_{name}"])
                                                     for name in (REGIONS if whole_v else (region,))),
                "n_vote_supported_positions": supported,
                "n_mutations": mutations if eligible else None,
                "shm_frequency": rate if eligible else None,
                "median_trust4_consensus_base_votes": votes if eligible else None}, key)
    assert set(member_shm) == expected_member_shm
    detail = index_rows(read_table(scope / "metrics/clone_shm.tsv"),
                        ("library_id", "clone_id", "isotype_group", "v_region"))
    expected_rows = set()
    cells = defaultdict(list)
    for key, records in units.items():
        for region in ("whole_v", *REGIONS):
            identity = (*key, region)
            values = [evidence[region][1] for _, evidence in records]
            evaluated = [value for value in values if value[1]]
            rates = [value[0] for value in evaluated]
            mutations = [value[2] for value in evaluated]
            zero = sum(value == 0 for value in mutations)
            positive = sum(value > 0 for value in mutations)
            median = statistics.median(rates) if rates else None
            cells[(key[0], region, key[2])].append({"median": median, "zero": zero, "n": len(rates),
                                                   "zero_median": statistics.median(mutations) == 0 if mutations else False,
                                                   "measured": sum(value[0] is not None for value in values),
                                                   "votes": [value[3] for value in evaluated if value[3] is not None]})
            if not rates:
                continue
            expected_rows.add(identity)
            assert identity in detail, (identity, "clone SHM row with eligible members is absent")
            state = ("all_zero" if zero == len(rates) else "all_positive" if positive == len(rates)
                     else "mixed")
            check_fields(detail[identity], {
                "consensus_count": sum(int(member["consensus_count"]) for member, _ in records),
                "eligible_consensus_count": sum(int(member["consensus_count"]) for member, evidence in records
                                                if evidence[region][1][1]),
                "n_members": len(records), "n_members_eligible": len(rates),
                "n_members_zero_shm": zero, "n_members_positive_shm": positive,
                "median_shm_frequency": median, "mutation_state": state}, identity)
    assert set(detail) == expected_rows
    summaries = read_table(scope / "metrics/shm.tsv")
    for row in summaries:
        key = (row["library_id"], row["v_region"], row["isotype_group"])
        clones = cells[key]
        medians = [clone["median"] for clone in clones if clone["median"] is not None]
        zero = sum(clone["zero_median"] for clone in clones)
        votes = [vote for clone in clones for vote in clone["votes"]]
        q25, q75 = quantile(medians, .25), quantile(medians, .75)
        check_fields(row, {"n_members_measured": sum(clone["measured"] for clone in clones),
                           "n_members_eligible": sum(clone["n"] for clone in clones),
                           "n_members_zero_shm": sum(clone["zero"] for clone in clones),
                           "n_clones_eligible": len(medians), "n_clones_zero_median_shm": zero,
                           "fraction_clones_zero_median_shm": zero / len(medians) if medians else None,
                           "mean_clone_median_shm": statistics.mean(medians) if medians else None,
                           "median_clone_median_shm": statistics.median(medians) if medians else None,
                           "q25_clone_median_shm": q25, "q75_clone_median_shm": q75,
                           "iqr_clone_median_shm": q75 - q25 if medians else None,
                           "median_trust4_consensus_base_votes": statistics.median(votes) if votes else None}, key)
    return len(detail), len(summaries)
