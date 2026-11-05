"""Compare reviewed values with scientific output tables using explicit tolerances."""

from pathlib import Path
from test.support.reference_tables import check_fields, number, read_table, index_rows
from decimal import Decimal, InvalidOperation
import math
import json
MISSING_VALUES = {"", "NA", "NaN", "nan", "N/A", "NULL", "null"}

def matches_expected(observed: str, expected: str, kind: str) -> bool:
    """Compare one published cell with its reviewed value under a declared type.

    A `string` expectation compares exactly, so leading zeros in an identifier
    stay significant and the missing-value tokens stay distinct from each other
    and from zero. An `integer` compares as a decimal, which remains exact above
    float precision. A `number` compares with relative tolerance 1e-8 and
    absolute tolerance 1e-10. A missing token never compares as a quantity.
    """
    if observed == expected:
        return True
    if kind == "string" or observed in MISSING_VALUES or expected in MISSING_VALUES:
        return False
    try:
        if kind == "integer":
            return Decimal(observed) == Decimal(expected)
        return math.isclose(float(observed), float(expected), rel_tol=1e-8, abs_tol=1e-10)
    except (ValueError, InvalidOperation):
        return False

def verify_expected_values(results: Path, examples: Path, scope: str | None = None,
                           libraries: set[str] | None = None) -> int:
    """Check the committed expectations against the tables this run produced.

    Each example names a table relative to the results root, a JSON selector
    that must match exactly one row, the column, the reviewed value and the type
    that decides how the cell is compared. Selector cells always compare as
    strings. `review_basis` records how a person reviewed the value and is not
    itself checked. Returns the number of expectations checked.
    """
    failures = []
    examples_rows = read_table(examples)
    if scope:
        examples_rows = [row for row in examples_rows
                         if Path(row["table"]).parts[:2] == ("scopes", scope)
                         or (Path(row["table"]).parts[0] == "libraries"
                             and Path(row["table"]).parts[1] in (libraries or set()))]
    assert examples_rows, "No committed expectations apply to this run"
    for example in examples_rows:
        context = f"{example['table']} {example['selector']} {example['column']}"
        try:
            published = read_table(results / example["table"])
        except OSError as error:
            failures.append(f"{context}: {error}")
            continue
        selector = json.loads(example["selector"])
        selected = [row for row in published
                    if all(row.get(column) == value for column, value in selector.items())]
        if len(selected) != 1:
            failures.append(f"{context}: selector matches {len(selected)} rows, not one")
            continue
        if example["column"] not in selected[0]:
            failures.append(f"{context}: the published table has no such column")
            continue
        observed = selected[0][example["column"]]
        if not matches_expected(observed, example["expected"], example["type"]):
            failures.append(f"{context}: observed={observed!r} "
                            f"expected={example['expected']!r}")
    if failures:
        raise AssertionError("Committed expectations failed:\n  " + "\n  ".join(failures))
    return len(examples_rows)
