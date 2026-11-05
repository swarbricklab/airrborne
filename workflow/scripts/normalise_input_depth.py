"""Convert a library's native depth record to the input-depth schema.

Reads the <library_id>_seq_depth.tsv that count_seq_depth.sh wrote beside
assembly and writes <library_id>_input_depth.tsv with library_id,
n_fastq_file_pairs, input_read_pairs and total_bases. It reads no FASTQ, so a
--post-trust4 run can reuse a retained depth record.
"""

import csv
from pathlib import Path

SOURCE_COLUMNS = ["gussid", "n_fastq_pairs", "total_read_pairs", "total_reads", "total_bases"]
DEPTH_COLUMNS = ["library_id", "n_fastq_file_pairs", "input_read_pairs", "total_bases"]


def normalise_input_depth(source: Path, destination: Path, library_id: str) -> None:
    """Rename the one retained row and drop its redundant total_reads cell.

    Retained rows must preserve the counter's paired-read invariant even when
    the original FASTQs are unavailable. The analytical schema keeps the
    read-pair denominator, in which both mates of a fragment are one unit.
    """
    with source.open(newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader, [])
        rows = list(reader)
    if header != SOURCE_COLUMNS or len(rows) != 1 or len(rows[0]) != len(header):
        raise ValueError(f"Invalid retained sequencing-depth table: {source}")
    if rows[0][0] != library_id:
        raise ValueError(f"Sequencing-depth library differs from {library_id}: {source}")
    if any(not value.isascii() or not value.isdigit() for value in rows[0][1:]):
        raise ValueError(f"Sequencing-depth counts must be non-negative integers: {source}")
    if int(rows[0][3]) != 2 * int(rows[0][2]):
        raise ValueError(f"Sequencing-depth total_reads must equal twice total_read_pairs "
                         f"for {library_id}: {source}")
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = destination.with_suffix(".tmp")
    try:
        with temporary.open("w", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(DEPTH_COLUMNS)
            writer.writerow(rows[0][:3] + rows[0][4:])
        temporary.replace(destination)
    finally:
        temporary.unlink(missing_ok=True)


if "snakemake" in globals():
    normalise_input_depth(Path(snakemake.input.depth), Path(snakemake.output.depth),
                           str(snakemake.wildcards.library_id))
    Path(snakemake.log[0]).write_text(f"Normalised retained depth for {snakemake.wildcards.library_id}\n")
