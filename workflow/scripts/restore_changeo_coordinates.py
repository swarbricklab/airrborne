"""Restore native coordinates for explicit Change-O reuse of a released AIRR file.

Writes supported IG rows to --output and all other rows to --unavailable.
The two outputs retain original extension evidence. No sequence or germline is
modified. Use the pinned IMGT-gapped reference with CreateGermlines.
"""
import argparse
from contextlib import ExitStack
import csv
import gzip
from pathlib import Path
import tempfile

COORDINATES = ("v_germline_start", "v_germline_end", "d_germline_start", "d_germline_end", "d_cigar")


def publish_outputs(temporary: list[Path], destinations: tuple[Path, Path]) -> None:
    """Restore the previous output pair if either replacement fails."""
    backups: dict[Path, Path] = {}
    published: list[Path] = []
    try:
        for destination in destinations:
            if destination.is_dir():
                raise IsADirectoryError(f"Output must be a file: {destination}")
        for destination in destinations:
            if destination.exists() or destination.is_symlink():
                with tempfile.NamedTemporaryFile(
                        dir=destination.parent, prefix=destination.name + ".backup.", delete=False) as handle:
                    backup = Path(handle.name)
                try:
                    destination.replace(backup)
                except BaseException:
                    backup.unlink()
                    raise
                backups[destination] = backup
        for path, destination in zip(temporary, destinations):
            path.replace(destination)
            published.append(destination)
    except BaseException:
        for destination in published:
            destination.unlink()
        for destination, backup in backups.items():
            backup.replace(destination)
        raise
    else:
        for backup in backups.values():
            backup.unlink()


def restore(source: Path, output: Path, unavailable: Path) -> tuple[int, int]:
    """Partition IG reconstruction inputs, preserving unsupported rows separately."""
    if len({path.resolve() for path in (source, output, unavailable)}) != 3:
        raise ValueError("Input, output and unavailable paths must differ")
    opener = gzip.open if source.suffix == ".gz" else open
    counts, temporary = [0, 0], []
    try:
        with ExitStack() as stack:
            handle = stack.enter_context(opener(source, "rt", newline=""))
            reader = csv.DictReader(handle, delimiter="\t")
            required = {"germline_reconstruction_status", *("changeo_" + field for field in COORDINATES)}
            if not required <= set(reader.fieldnames or []):
                raise ValueError("AIRR file lacks native reconstruction extensions")
            writers = []
            for destination in (output, unavailable):
                stream = stack.enter_context(tempfile.NamedTemporaryFile(
                    mode="w", newline="", dir=destination.parent, prefix=destination.name + ".", delete=False))
                temporary.append(Path(stream.name))
                writer = csv.DictWriter(stream, fieldnames=reader.fieldnames, delimiter="\t", lineterminator="\n")
                writer.writeheader()
                writers.append(writer)
            for row in reader:
                if None in row or None in row.values():
                    raise ValueError(f"Malformed AIRR row {reader.line_num}")
                group = 0 if row["germline_reconstruction_status"] == "available" else 1
                if group == 0:
                    for field in COORDINATES:
                        row[field] = row["changeo_" + field]
                writers[group].writerow(row)
                counts[group] += 1
        publish_outputs(temporary, (output, unavailable))
    finally:
        for path in temporary:
            path.unlink(missing_ok=True)
    return tuple(counts)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--unavailable", required=True, type=Path)
    args = parser.parse_args()
    supported, unsupported = restore(args.input, args.output, args.unavailable)
    print(f"Wrote {supported} supported IG inputs and {unsupported} unavailable rows")
