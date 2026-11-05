"""Reconstruct inputs results independently from retained scientific inputs."""

from pathlib import Path
from test.support.reference_tables import check_fields, number, read_table, index_rows
import csv
import re


def verify_input_depth(root: Path, libraries: set[str]) -> int:
    """Require one native and analytical depth row per selected library."""
    native_columns = ("gussid", "n_fastq_pairs", "total_read_pairs", "total_reads", "total_bases")
    analytical_columns = ("library_id", "n_fastq_file_pairs", "input_read_pairs", "total_bases")
    for library in sorted(libraries):
        directory = root / "libraries" / library / "metrics"
        tables = []
        for suffix, columns in (("seq_depth", native_columns), ("input_depth", analytical_columns)):
            with (directory / f"{library}_{suffix}.tsv").open(newline="") as handle:
                reader = csv.DictReader(handle, delimiter="\t")
                assert tuple(reader.fieldnames or ()) == columns, f"{library}: invalid depth schema"
                rows = list(reader)
            assert len(rows) == 1, f"{library}: expected one depth row"
            row = rows[0]
            assert tuple(row) == columns, f"{library}: invalid depth row"
            assert row[columns[0]] == library, f"{library}: incorrect depth library"
            assert all(re.fullmatch(r"[0-9]+", row[column] or "") for column in columns[1:]), \
                f"{library}: depth counts must be non-negative integers"
            tables.append(row)
        native, adapted = tables
        assert int(native["total_reads"]) == 2 * int(native["total_read_pairs"]), \
            f"{library}: total_reads must equal twice total_read_pairs"
        for column, source in (("n_fastq_file_pairs", "n_fastq_pairs"),
                               ("input_read_pairs", "total_read_pairs"),
                               ("total_bases", "total_bases")):
            assert int(adapted[column]) == int(native[source]), f"{library}: incorrect {column}"
    return len(libraries)


def verify_yield(root: Path, scope: Path) -> int:
    """Rebuild the yield table's depth columns and yields from each library's input depth."""
    rows = read_table(scope / "metrics/yield.tsv")
    for row in rows:
        library = row["library_id"]
        depth, = read_table(root / "libraries" / library / "metrics" / f"{library}_input_depth.tsv")
        pairs = int(depth["input_read_pairs"])
        expected = {column: depth[column] for column in ("n_fastq_file_pairs", "input_read_pairs")}
        expected["input_bases"] = depth["total_bases"]
        expected["yield_status"] = "available" if pairs else "no_input_read_pairs"
        for locus in ("igh", "igk", "igl", "tra", "trb", "trd", "trg"):
            expected[f"{locus}_yield_per_million_read_pairs"] = int(row[f"{locus}_trust4_consensus_count"]) / pairs * 1e6 if pairs else None
        check_fields(row, expected, library)
    return len(rows)
