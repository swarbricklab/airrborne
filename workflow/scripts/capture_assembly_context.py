#!/usr/bin/env python3
"""Capture source state only after a new assembly has actually executed.

The assembly image receipt, when given, supplies the identity of the image the
job ran, so the producer record names the actual SIF or Docker image.
"""

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from workflow.lib.provenance import capture_source

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--repo", type=Path, required=True,
                    help="the pipeline checkout whose source state is recorded")
parser.add_argument("--context", required=True,
                    help="the run_trust4 producer context as JSON; the source state is added "
                         "and the result is printed")
parser.add_argument("--image-receipt", type=Path,
                    help="the assembly image receipt obtain_image.py wrote; its identity replaces "
                         "any image in the context")
args = parser.parse_args()
context = json.loads(args.context)
context["source"] = capture_source(args.repo)
if args.image_receipt:
    context["image"] = json.loads(args.image_receipt.read_text())
print(json.dumps(context, sort_keys=True))
