#!/usr/bin/env python3
"""Render prepared scope data as a portable folder of local browser assets."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import shutil
import tempfile
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from workflow.lib.release_bundle import release_provenance
from workflow.lib.release_methods import read_provenance, read_methods_inputs, methods_paragraphs

REPORT = Path(__file__).resolve().parents[1] / "report"


def render(data_dir: Path, output_dir: Path, provenance: Path, effective_config: Path,
           manifest: Path, threshold: Path) -> None:
    """Build the report folder in a temporary directory, then replace output_dir.

    Adds the compact provenance, with its files named at their release paths
    as the released provenance.json names them, and the methods paragraphs to
    the prepared payload, and writes it as data.js beside index.html, the
    assets, curves.js and one script per library detail and per locus gene-use
    payload.
    """
    data = json.loads((data_dir / "report.json").read_text())
    id_column = data["libraries"]["columns"].index("library_id")
    library_ids = {row[id_column] for row in data["libraries"]["rows"]}
    record = release_provenance(read_provenance(provenance, data["scope_id"], library_ids))
    config, libraries, calibration = read_methods_inputs(effective_config, manifest, threshold)
    data["provenance"] = record
    data["methods"] = methods_paragraphs(record, config, libraries, calibration)
    output_dir.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=".report-", dir=output_dir.parent) as directory:
        staging = Path(directory) / "report"
        assets = staging / "assets"
        assets.mkdir(parents=True)
        # plotly.min.js is plotly.js 3.7.0 (MIT licence), the bundle plotly 6.9.0's
        # get_plotlyjs() returns. The browser draws every figure from it.
        for name in ("plotly.min.js", "LICENSE.plotly.txt", "selection_controller.js",
                     "report_view.js", "report_sections.js", "report.css"):
            shutil.copyfile(REPORT / "assets" / name, assets / name)
        # Cohort rarefaction loads curve rows without every library's clone detail.
        cohort_curves = {"columns": [], "rows": []}
        for library_id, filename in list(data["detail_files"].items()):
            payload = json.loads((data_dir / filename).read_text())
            if payload["curves"]["rows"]:
                if not cohort_curves["columns"]:
                    cohort_curves["columns"] = payload["curves"]["columns"]
                if cohort_curves["columns"] != payload["curves"]["columns"]:
                    raise ValueError("Library rarefaction curves must have the same columns.")
                cohort_curves["rows"].extend(payload["curves"]["rows"])
            target = Path(filename).with_suffix(".js")
            (staging / target).parent.mkdir(parents=True, exist_ok=True)
            (staging / target).write_text(f"window.REPORT_DETAIL[{json.dumps(library_id)}] = {json.dumps(payload, ensure_ascii=True, allow_nan=False, separators=(',', ':'))};\n", encoding="utf-8")
            data["detail_files"][library_id] = str(target)
        data["curve_file"] = "curves.js"
        (staging / data["curve_file"]).write_text(
            'window.REPORT_CURVES["cohort"] = ' +
            json.dumps(cohort_curves, ensure_ascii=True, allow_nan=False, separators=(',', ':')) + ";\n",
            encoding="utf-8",
        )
        for locus, filename in list(data["gene_usage_files"].items()):
            payload = json.loads((data_dir / filename).read_text())
            target = Path(filename).with_suffix(".js")
            (staging / target).parent.mkdir(parents=True, exist_ok=True)
            (staging / target).write_text(f"window.REPORT_GENES[{json.dumps(locus)}] = {json.dumps(payload, ensure_ascii=True, allow_nan=False, separators=(',', ':'))};\n", encoding="utf-8")
            data["gene_usage_files"][locus] = str(target)
        (staging / "data.js").write_text(f"window.REPORT_DATA = {json.dumps(data, ensure_ascii=True, allow_nan=False, separators=(',', ':'))};\n", encoding="utf-8")
        shutil.copyfile(REPORT / "index.html", staging / "index.html")
        if output_dir.exists():
            shutil.rmtree(output_dir)
        staging.rename(output_dir)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", type=Path, required=True,
                        help="the payload directory prepare_scope_report.py wrote")
    parser.add_argument("--output-dir", type=Path, required=True,
                        help="report directory to replace")
    parser.add_argument("--provenance", type=Path, required=True, help="the scope's provenance.json")
    parser.add_argument("--effective-config", type=Path, required=True,
                        help="the scope's effective_config.yaml")
    parser.add_argument("--manifest", type=Path, required=True, help="the scope manifest")
    parser.add_argument("--threshold", type=Path, required=True,
                        help="the scope's one-row operative clone threshold table")
    args = parser.parse_args()
    render(args.data_dir, args.output_dir, args.provenance, args.effective_config, args.manifest, args.threshold)
