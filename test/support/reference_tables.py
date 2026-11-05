"""Read and compare reference tables without production analysis imports."""

import csv
from pathlib import Path

MISSING = {"", "NA", "NaN", "nan"}
TOLERANCE = 1e-9

def read_table(path: Path) -> list[dict[str, str]]:
    with Path(path).open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))

def number(value: str) -> float | None:
    return None if value in MISSING else float(value)

def check_fields(row: dict[str, str], expected: dict, context: object) -> None:
    for column, value in expected.items():
        observed = row[column]
        if value is None:
            matched = observed in MISSING
        elif isinstance(value, str):
            matched = observed == value
        elif isinstance(value, bool):
            matched = observed == str(value).upper()
        else:
            parsed = number(observed)
            matched = parsed is not None and abs(parsed - value) <= TOLERANCE
        if not matched:
            raise AssertionError(f"{context}, {column}: observed {observed!r}, expected {value!r}")

def index_rows(rows: list[dict[str, str]], columns: list[str] | tuple[str, ...]) -> dict:
    result = {}
    for row in rows:
        key = tuple(row[column] for column in columns)
        if key in result:
            raise AssertionError(f"duplicate key: {key}")
        result[key] = row
    return result
