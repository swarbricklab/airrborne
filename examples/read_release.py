"""Read a release without changing identifiers or multiplying member support.

Run from the pipeline checkout:
python examples/read_release.py RELEASE --library SYN_D1_T1
Add --metadata FILE for a CSV with a unique tissue_id column.
"""

import argparse
import csv
import gzip
from pathlib import Path


def read_rows(path, delimiter="\t"):
    """Yield string cells, preserving literal NA and leading zeroes; blank stays blank."""
    path = Path(path)
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", newline="", encoding="utf-8-sig") as handle:
        yield from csv.DictReader(handle, delimiter=delimiter)


def number(value):
    """Convert a numeric blank to None after identifier columns have been preserved."""
    return None if value == "" else float(value)


def unique_index(rows, keys):
    result = {}
    for row in rows:
        key = tuple(row[name] for name in keys)
        if any(value == "" for value in key) or key in result:
            raise ValueError(f"Empty or duplicate lookup key: {key}")
        result[key] = row
    return result


def add_metadata(headline, metadata, key="tissue_id"):
    """Retain every headline row and reject metadata that would multiply it."""
    lookup = unique_index(metadata, (key,))
    columns = {name for row in metadata for name in row if name != key}
    overlap = columns.intersection(name for row in headline for name in row)
    for row in headline:
        matched = lookup.get((row[key],))
        if matched is not None and any(matched[name] != row[name] for name in overlap):
            raise ValueError("Metadata conflicts with scientific columns")
    columns -= overlap
    joined = [{**row, **{name: lookup.get((row[key],), {}).get(name, "")
                        for name in columns}} for row in headline]
    assert len(joined) == len(headline)
    assert [row["library_id"] for row in joined] == [row["library_id"] for row in headline]
    return joined


def annotate_vj_pairs(pairs, gene_usage):
    """Map V and J entities separately with a unique lookup in each library and locus."""
    context = ("scope_id", "library_id", "locus")
    lookup = unique_index(gene_usage, (*context, "gene_type", "gene_key"))
    result = []
    for row in pairs:
        out = dict(row)
        for segment in ("V", "J"):
            key = (*[row[name] for name in context], segment, row[segment.lower() + "_gene_key"])
            out[segment.lower() + "_gene_equivalence_group"] = lookup[key]["gene_equivalence_group"]
        result.append(out)
    assert len(result) == len(pairs)
    assert sum(int(row["consensus_count"]) for row in result) == sum(int(row["consensus_count"]) for row in pairs)
    return result


def join_member_airr(members, airr):
    """Join one member to its representative, retaining transferred support once."""
    lookup = unique_index((r for r in airr if r["is_member_representative"] == "T"),
                          ("library_id", "sequence_id"))
    unique_index(members, ("scope_id", "sequence_id"))
    joined = []
    for row in members:
        record = lookup[(row["library_id"], row["sequence_id"])]
        assert int(record["member_consensus_count"]) == int(row["member_consensus_count"])
        joined.append({**row, "sequence": record["sequence"],
                       "source_consensus_count": record["consensus_count"]})
    assert len(joined) == len(members)
    assert sum(int(r["member_consensus_count"]) for r in joined) == sum(int(r["member_consensus_count"]) for r in members)
    return joined


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("release", type=Path)
    parser.add_argument("--library", default="SYN_D1_T1")
    parser.add_argument("--metadata", type=Path)
    args = parser.parse_args()
    headline = list(read_rows(args.release / "metrics/headline.tsv"))
    unique_index(headline, ("scope_id", "library_id"))
    if args.metadata:
        headline = add_metadata(headline, list(read_rows(args.metadata, delimiter=",")))
    selected = [row for row in headline if row["library_id"] == args.library]
    if len(selected) != 1:
        raise ValueError("Select exactly one library in this release")
    row = selected[0]
    print("Library:", row["library_id"], "input pairs:", number(row["input_read_pairs"]))
    print("IGH yield per million pairs:", number(row["igh_yield_per_million_read_pairs"]))
    print("Grouped IGH support:", number(row["igh_grouped_consensus_count"]))
    shm = [r for r in read_rows(args.release / "metrics/shm.tsv")
           if r["library_id"] == args.library and r["isotype_group"] == "all_IGH" and r["v_region"] == "whole_v"]
    assert len(shm) == 1
    print("Whole-V eligible/measured members:", shm[0]["n_members_eligible"], shm[0]["n_members_measured"])
    members = [r for r in read_rows(args.release / "metrics/detail/igh/clone_members.tsv")
               if r["library_id"] == args.library]
    airr = read_rows(args.release / "airr" / (args.library + "_rearrangements.tsv.gz"))
    joined = join_member_airr(members, airr)
    print("IGH members:", len(joined), "transferred support:", sum(int(r["member_consensus_count"]) for r in joined))
    pairs = [r for r in read_rows(args.release / "metrics/vj_pair_usage.tsv")
             if r["library_id"] == args.library and r["locus"] == "IGH"]
    genes = [r for r in read_rows(args.release / "metrics/gene_usage.tsv")
             if r["library_id"] == args.library and r["locus"] == "IGH"]
    annotate_vj_pairs(pairs, genes)


if __name__ == "__main__":
    main()
