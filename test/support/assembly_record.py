"""Expected assembly records for producer and publication tests."""

import copy
import hashlib
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any


def file_identity(path: Path) -> dict[str, Any]:
    """Describe the constructed file independently of the production recorder."""
    content = path.read_bytes()
    return {"name": path.name, "sha256": hashlib.sha256(content).hexdigest(),
            "size_bytes": len(content)}


def build_assembly_provenance(library_id: str, *, source: Mapping[str, Any],
                              image: Mapping[str, Any], settings: Mapping[str, Any],
                              reference: Path, read_pairs: Sequence[Mapping[str, str]],
                              outputs: Mapping[str, Path],
                              software: Sequence[Mapping[str, Any]] = ()) -> dict[str, Any]:
    """Construct the expected record for comparison with the Perl recorder."""
    return {"schema_version": "2.0", "stage": "assembly", "library_id": library_id,
            "producer": {"source": copy.deepcopy(dict(source)), "image": copy.deepcopy(dict(image)),
                         "settings": copy.deepcopy(dict(settings)), "software": copy.deepcopy(list(software))},
            "inputs": {"read_pairs": copy.deepcopy(list(read_pairs)), "reference": file_identity(reference)},
            "outputs": {name: file_identity(path) for name, path in sorted(outputs.items())}}
