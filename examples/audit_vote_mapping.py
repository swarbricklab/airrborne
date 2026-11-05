"""Count unavailable vote attribution in an originating results tree.

The portable release omits vote vectors. Supply outs.results_dir from the run.
Counts describe common-QC source records, not distinct grouped members or cells.
"""
import argparse
from collections import Counter
import json
from pathlib import Path


def audit(results):
    counts = Counter()
    libraries = []
    for path in sorted((Path(results) / "libraries").glob("*/airr/*_igh_consensus_votes.tsv")):
        local = Counter()
        with path.open("rb") as handle:
            header = handle.readline().rstrip(b"\r\n").split(b"\t")
            if header[-1] != b"junction_vote_status":
                raise ValueError(f"Unexpected vote-table column order: {path}")
            mapped_column = header.index(b"imgt_v_positions_mapped")
            for line in handle:
                local["igh_source_records"] += 1
                if line.rstrip(b"\r\n").endswith(b"\topposite_contig_orientation"):
                    local["opposite_contig_orientation"] += 1
                    if line.split(b"\t", mapped_column + 1)[mapped_column] != b"0":
                        raise ValueError(f"Opposite orientation retained mapped positions: {path}")
        counts.update(local)
        libraries.append({"library_id": path.parent.parent.name, **local})
    if not libraries:
        raise ValueError("No IGH consensus-vote tables found in the results tree")
    counts.setdefault("opposite_contig_orientation", 0)
    return {"unit": "common-QC IGH source records", "n_libraries": len(libraries),
            "counts": dict(counts), "libraries": libraries}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("results", type=Path)
    args = parser.parse_args()
    print(json.dumps(audit(args.results), indent=2))
