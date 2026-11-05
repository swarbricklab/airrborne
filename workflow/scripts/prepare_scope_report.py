#!/usr/bin/env python3
"""Prepare one scope's supplied measurements for the offline report."""

import argparse
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from workflow.lib.report_data import prepare


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--scope-dir", type=Path, required=True,
                        help="scope directory holding manifest/ and metrics/")
    parser.add_argument("--output-dir", type=Path, required=True,
                        help="directory for report.json and the details/ and genes/ payloads")
    parser.add_argument("--metadata-sources", type=json.loads, default={},
                        help="JSON declarations for this scope's tissue and library report metadata.")
    parser.add_argument("--calibration", type=Path, nargs="*", default=[],
                        help="the scope's clone threshold, donor threshold and distance diagnostic tables")
    args = parser.parse_args()
    prepare(args.scope_dir, args.output_dir, args.metadata_sources, args.calibration)
