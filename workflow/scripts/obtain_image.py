#!/usr/bin/env python3
"""Obtain a pinned image for the selected runtime and record its identity.

The obtain_analysis_image and obtain_assembly_image rules run this script.
Apptainer and Singularity run a SIF pulled from the reference into --sif when
that file is absent. Docker pulls the reference into its daemon and ignores
--sif. The receipt at --receipt holds the identity every rule that runs the
image declares as its image input. It is rewritten only when its content
changes.
"""

import argparse
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from workflow.lib.provenance import obtain_image


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--reference", required=True,
                        help="digest-pinned image reference, such as immcantation/suite@sha256:...")
    parser.add_argument("--sif", type=Path, required=True,
                        help="local SIF path for Apptainer and Singularity, pulled when absent")
    parser.add_argument("--receipt", type=Path, required=True,
                        help="identity record to write, the image input that rules running it declare")
    parser.add_argument("--runtime", choices=("apptainer", "singularity", "docker"),
                        help="runtime to obtain the image for; the default follows container_exec.sh")
    args = parser.parse_args()
    try:
        identity = obtain_image(args.reference, args.sif, args.receipt, runtime=args.runtime)
    except (OSError, ValueError, subprocess.CalledProcessError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        sys.exit(1)
    print(f"{identity['runtime']} runs {identity['configured_reference']} as {identity['format']}; "
          f"receipt {args.receipt}")
