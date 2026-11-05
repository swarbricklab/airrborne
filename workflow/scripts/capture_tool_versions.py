#!/usr/bin/env python3
"""Measure installed tools once per immutable image, outside scientific inputs."""

import argparse
import json
import sys
import re
import subprocess
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from workflow.lib.provenance import capture_tool_versions, identify_image, version_queries


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--image", required=True,
                        help="container_exec.sh's first argument: a local SIF path or a docker:// reference")
    parser.add_argument("--docker-image", required=True, help="the digest-pinned reference the runtime obtains")
    parser.add_argument("--image-record", type=Path,
                        help="identity receipt already written for this image; without it the image is identified here")
    parser.add_argument("--role", choices=("analysis", "assembly"), required=True,
                        help="which image's version queries to run")
    parser.add_argument("--cache-dir", type=Path, required=True,
                        help="directory holding one measured record per image digest")
    parser.add_argument("--output", type=Path, required=True, help="the version record to write")
    parser.add_argument("--container-exec", type=Path,
                        default=Path(__file__).with_name("container_exec.sh"),
                        help="the container wrapper; defaults to container_exec.sh beside this script")
    parser.add_argument("--refresh", action="store_true",
                        help="measure again and replace the cached record for this image")
    args = parser.parse_args()
    try:
        if args.image_record:
            image = json.loads(args.image_record.read_text())
        else:
            image = identify_image(args.image, args.docker_image)
        versions = capture_tool_versions(args.cache_dir, image,
                                         [str(args.container_exec.resolve()), args.image, args.docker_image],
                                         args.role, refresh=args.refresh)
    except (OSError, ValueError, subprocess.CalledProcessError) as error:
        if args.role != "assembly":
            raise
        # Analysis-only handover must remain possible without the assembly image.
        # Its configured digest identifies a requested tool, not a reused producer.
        digest = re.search(r"@sha256:([0-9a-f]{64})$", args.docker_image)
        if not digest:
            raise
        image = {"format": "oci", "digest": "sha256:" + digest.group(1),
                 "configured_reference": args.docker_image, "status": "unavailable"}
        versions = {"query_schema_version": 1, "image": image, "software": [
            {"name": name, "version": None, "status": "unknown", **query,
             "error": {"message": "Configured assembly image is unavailable: " + str(error)}}
            for name, query in version_queries("assembly").items()]}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(versions, indent=2, sort_keys=True) + "\n")
