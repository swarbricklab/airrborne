"""Describe a release from recorded inputs, without changing analysis dependencies.

The module identifies the pipeline source, the pinned images and the reference
files, obtains the analysis image, measures tool versions and builds the
per-library processing identity and the scope provenance record. Run image
queries once outside scientific rules, then pass their cached records to each
scope. An assembly producer record is written when assembly finishes.
Packaging cannot recover missing history from the current configuration.
"""

import copy
import hashlib
import json
import os
import re
import shutil
import subprocess
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any


SCOPE_SCHEMA_VERSION = "3.2"
ASSEMBLY_SCHEMA_VERSION = "2.0"
SEED = 1337
SOURCE_PATHS = ("Snakefile", "workflow", "envs",
                "profiles", "run_module.sh")
IMAGE_SETTINGS = ("immcantation_docker", "trust4_container")
RUNTIMES = ("apptainer", "singularity", "docker")


def _json(value: Any) -> str:
    return json.dumps(value, sort_keys=True, ensure_ascii=True, allow_nan=False, separators=(",", ":"))


def _hash(value: bytes) -> str:
    return hashlib.sha256(value).hexdigest()


def file_identity(path: Path) -> dict[str, Any]:
    """Hash an actual reference, image or output, not a configured name."""
    path = Path(path)
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return {"name": path.name, "sha256": digest.hexdigest(), "size_bytes": path.stat().st_size}


def _git(repo: Path, *arguments: str) -> bytes:
    return subprocess.run(["git", "-C", str(repo), *arguments], check=True,
                          stdout=subprocess.PIPE, stderr=subprocess.PIPE).stdout


def capture_source(repo: Path, *, relevant_paths: Sequence[str] = SOURCE_PATHS) -> dict[str, Any]:
    """Identify pipeline source, including tracked edits and untracked code.

    Configuration and metadata are data inputs. Effective scientific settings
    are recorded separately, while optional metadata belongs to the report.
    The dirty flag describes the source paths below, not unrelated repository
    files. Source edits remain identifiable beyond their base revision.
    """
    repo = Path(repo).resolve()
    revision = _git(repo, "rev-parse", "HEAD").decode().strip()
    tracked = _git(repo, "status", "--porcelain=v1", "-z", "--untracked-files=no", "--", *relevant_paths)
    working_diff = _git(repo, "diff", "--binary", "HEAD", "--", *relevant_paths)
    index_diff = _git(repo, "diff", "--cached", "--binary", "--", *relevant_paths)
    changed = _git(repo, "diff", "--name-only", "-z", "HEAD", "--", *relevant_paths)
    staged = _git(repo, "diff", "--cached", "--name-only", "-z", "--", *relevant_paths)
    untracked = _git(repo, "ls-files", "--others", "--exclude-standard", "-z", "--", *relevant_paths)
    changed_names = {name.decode() for name in (changed + staged).split(b"\0") if name}
    untracked_names = {name.decode() for name in untracked.split(b"\0") if name}
    changes = []
    for name in sorted(changed_names | untracked_names):
        path = repo / name
        entry = {"path": name, "status": "untracked" if name in untracked_names else "tracked"}
        if path.is_symlink():
            entry["symlink_target_sha256"] = _hash(str(path.readlink()).encode())
        elif path.is_dir():
            # A changed gitlink refers to the submodule's actual checkout.
            entry.update(revision=_git(path, "rev-parse", "HEAD").decode().strip(),
                         diff_sha256=_hash(_git(path, "diff", "--binary", "HEAD")))
        elif path.exists():
            entry.update(file_identity(path))
        else:
            entry["status"] = "deleted"
        changes.append(entry)
    state = {"revision": revision, "dirty": bool(tracked or untracked_names),
             "tracked_dirty": bool(tracked), "working_diff_sha256": _hash(working_diff),
             "index_diff_sha256": _hash(index_diff), "relevant_changes": changes}
    state["state_sha256"] = _hash(_json(state).encode())
    return state


def select_runtime(runtime: str | None = None) -> str:
    """Name the runtime container_exec.sh selects: an explicit choice, else the first installed."""
    runtime = runtime or os.environ.get("AIRRBORNE_CONTAINER_RUNTIME") or next(
        (name for name in RUNTIMES if shutil.which(name)), None)
    if runtime not in RUNTIMES:
        raise ValueError("Select an installed apptainer, singularity or docker runtime")
    return runtime


def configured_digest(docker_reference: str) -> str:
    """Return the hex digest a pinned reference names; a tag is not an identity."""
    match = re.search(r"@sha256:([0-9a-f]{64})$", docker_reference)
    if not match:
        raise ValueError("Image identity requires a digest-pinned Docker reference")
    return match.group(1)


def image_sif_path(resources_dir: Path, docker_reference: str) -> Path:
    """Name the SIF the workflow pulls from a pinned reference, beside its other obtained resources."""
    repository = re.sub(r"[^A-Za-z0-9_.-]", "_", docker_reference.split("@", 1)[0])
    return Path(resources_dir) / f"{repository}-{configured_digest(docker_reference)[:12]}.sif"


def image_receipt_path(resources_dir: Path, docker_reference: str) -> Path:
    """Name the identity receipt written beside the SIF, which the rules running the image declare.

    The receipt stays with the obtained resources rather than in a results
    directory. A fresh results directory therefore reuses it, and the shared
    reference files that declare it are not rebuilt for each new run.
    """
    return image_sif_path(resources_dir, docker_reference).with_suffix(".json")


def _sif_identity(path: Path, docker_reference: str) -> dict[str, Any]:
    """Hash a SIF once and require the pinned reference it was pulled from in its metadata."""
    needle = docker_reference.encode()
    digest = hashlib.sha256()
    found = False
    carry = b""
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
            if not found:
                found = needle in carry + block
                carry = block[-(len(needle) - 1):]
    if not found:
        raise ValueError(f"{path} does not record the configured reference {docker_reference}; "
                         "remove the file so the workflow obtains the image again")
    return {"name": path.name, "sha256": digest.hexdigest(), "size_bytes": path.stat().st_size}


def identify_image(image: str, docker_reference: str, *, runtime: str | None = None) -> dict[str, Any]:
    """Identify a pinned image by its OCI digest and record how the selected runtime materialises it.

    `image` is container_exec.sh's first argument: a local SIF path or the
    docker:// form of the reference. The digest is the identity in every
    runtime. The format and its evidence (Docker image ID, SIF hash, or the
    registry manifest) describe this run and stay out of the identity key.
    """
    identity = {"runtime": select_runtime(runtime), "configured_reference": docker_reference,
                "digest": "sha256:" + configured_digest(docker_reference)}
    if identity["runtime"] == "docker":
        result = subprocess.run(["docker", "image", "inspect", "--format", "{{json .}}", docker_reference],
                                check=True, capture_output=True, text=True)
        inspected = json.loads(result.stdout)
        if docker_reference not in inspected.get("RepoDigests", []):
            raise ValueError("The inspected image does not carry the configured immutable digest")
        identity.update(format="docker", image_id=inspected["Id"], architecture=inspected["Architecture"],
                        os=inspected["Os"])
    elif image == "docker://" + docker_reference:
        # The runtime resolves this immutable registry manifest itself.
        identity["format"] = "oci"
    elif Path(image).is_file():
        identity.update(format="sif", file=_sif_identity(Path(image), docker_reference))
    else:
        raise ValueError(f"The configured SIF is absent: {image}")
    return identity


def obtain_image(docker_reference: str, sif: Path, receipt: Path, *,
                 runtime: str | None = None) -> dict[str, Any]:
    """Obtain a pinned image for the selected runtime and record its identity once.

    The workflow obtains the analysis and assembly images this way, from the
    Snakemake main process.

    Apptainer and Singularity run `sif`, pulled from the reference when the file
    is absent. The pull disables the runtime's own cache so image layers land in
    the temporary directory rather than the home directory. Docker pulls the
    reference into its daemon when it is absent and ignores `sif`. The receipt
    is rewritten only when its content changes, so a repeated run leaves the
    outputs that declare it valid.
    """
    runtime = select_runtime(runtime)
    sif, receipt = Path(sif), Path(receipt)
    if runtime == "docker":
        present = subprocess.run(["docker", "image", "inspect", docker_reference],
                                 capture_output=True).returncode == 0
        if not present:
            subprocess.run(["docker", "pull", "--platform", "linux/amd64", docker_reference], check=True)
    elif not sif.is_file():
        sif.parent.mkdir(parents=True, exist_ok=True)
        temporary = sif.with_name(sif.name + ".tmp")
        temporary.unlink(missing_ok=True)
        subprocess.run([runtime, "pull", "--disable-cache", "--arch", "amd64", str(temporary),
                        "docker://" + docker_reference], check=True)
        temporary.replace(sif)
    identity = identify_image(str(sif), docker_reference, runtime=runtime)
    rendered = json.dumps(identity, indent=2, sort_keys=True) + "\n"
    if not (receipt.exists() and receipt.read_text() == rendered):
        receipt.parent.mkdir(parents=True, exist_ok=True)
        temporary = receipt.with_name(receipt.name + ".tmp")
        temporary.write_text(rendered)
        temporary.replace(receipt)
    return identity


def image_cache_key(image: Mapping[str, Any]) -> str:
    """Key an image by its immutable digest alone, not by runtime, format or path."""
    if not re.fullmatch(r"sha256:[0-9a-f]{64}", str(image.get("digest", ""))):
        raise ValueError("An actual immutable image identity is required")
    return _hash(_json({"digest": image["digest"]}).encode())


def library_processing_identity(code_files: Mapping[str, Path], settings: Mapping[str, Any],
                                 image: Mapping[str, Any], references: Mapping[str, Path]) -> str:
    """Identify only the code, settings, image and resources consumed by library analysis.

    The caller supplies the actual library modules and scientific parameter
    subset. Do not pass the whole configuration, scope thresholds, metadata or
    the descriptive repository state: those would change IDs for unrelated work.
    """
    components = {"code": {name: file_identity(path)["sha256"] for name, path in sorted(code_files.items())},
                  "settings": settings, "image": image_cache_key(image),
                  "references": {name: file_identity(path)["sha256"] for name, path in sorted(references.items())}}
    return _hash(_json(components).encode())


def version_queries(role: str) -> dict[str, dict[str, Any]]:
    """Commands interrogate installed programs inside the selected image."""
    if role == "assembly":
        # run-trust4 has no --version flag: that command prints a banner and
        # fails. Read the installed launcher's own declaration with a command
        # that fails unless the declaration was actually found.
        command = r'''perl -ne 'if (/^my \$version = "([^"]+)"/) { print "$1\n"; $found=1; last } END { exit($found ? 0 : 1) }' "$(command -v run-trust4)"'''
        return {"TRUST4": {"command": ["bash", "-c", command],
                           "method": "installed run-trust4 version declaration"}}
    if role != "analysis":
        raise ValueError("Tool version role must be analysis or assembly")
    queries = {
        "IgBLAST": {"command": ["igblastn", "-version"], "pattern": r"igblastn:\s+(\S+)",
                    "method": "program version command"},
        "Change-O": {"command": ["python3", "-c", "from importlib.metadata import version; print(version('changeo'))"],
                     "method": "installed Python package version"},
        "airr-python": {"command": ["python3", "-c", "from importlib.metadata import version; print(version('airr'))"],
                        "method": "installed Python package version"},
        "R": {"command": ["Rscript", "--vanilla", "-e", "cat(paste(R.version$major, R.version$minor, sep='.'))"],
              "method": "R version object"},
    }
    for package in ("alakazam", "shazam", "scoper", "airr"):
        name = "airr-R" if package == "airr" else package
        queries[name] = {"command": ["Rscript", "--vanilla", "-e",
                                      f"cat(as.character(packageVersion('{package}')))"],
                         "method": "installed R package version"}
    return queries


def capture_tool_versions(cache_dir: Path, image: Mapping[str, Any],
                           command_prefix: Sequence[str], role: str, *,
                           refresh: bool = False) -> dict[str, Any]:
    """Cache measured versions once per image; failed queries remain unknown.

    command_prefix is container_exec.sh followed by the SIF and Docker names.
    Its runtime must match identify_image(). A dedicated workflow job owns this
    cache capture; scope packaging only reads the resulting record.
    """
    key = image_cache_key(image)
    path = Path(cache_dir) / (key + ".json")
    queries = version_queries(role)
    if path.exists() and not refresh:
        record = json.loads(path.read_text())
        if image_cache_key(record["image"]) != key or record["query_schema_version"] != 1:
            raise ValueError("Tool version cache identity or format differs")
        if {entry["name"] for entry in record["software"]} != set(queries):
            raise ValueError("Tool version cache has a different query set; refresh it explicitly")
        if any(any(entry.get(field) != value for field, value in queries[entry["name"]].items())
               for entry in record["software"]):
            raise ValueError("Tool version cache has a different query definition; refresh it explicitly")
        return record
    software = []
    for name, query in queries.items():
        entry = {"name": name, "version": None, "status": "unknown", **query}
        try:
            result = subprocess.run([*command_prefix, *query["command"]], capture_output=True,
                                    text=True, timeout=60, check=False)
            output = result.stdout.strip()
            match = re.search(query["pattern"], output) if "pattern" in query else None
            version = match.group(1) if match else output if "pattern" not in query else ""
            if result.returncode != 0:
                entry["error"] = {"exit_code": result.returncode, "stderr": result.stderr.strip()}
            elif (not version or version.lower() in {"unknown", "0+unknown", "na"}
                  or any(character.isspace() for character in version)):
                entry["error"] = {"message": "The version query did not return one recognised version value"}
            else:
                entry.update(version=version, status="measured")
        except (OSError, subprocess.TimeoutExpired) as error:
            entry["error"] = {"message": str(error)}
        software.append(entry)
    record = {"query_schema_version": 1, "image": copy.deepcopy(dict(image)), "software": software}
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(".json.tmp")
    temporary.write_text(json.dumps(record, indent=2, sort_keys=True, allow_nan=False) + "\n")
    temporary.replace(path)
    return record


def effective_scope_config(config: Mapping[str, Any], scope_id: str) -> dict[str, Any]:
    """Copy the effective scientific settings; report annotations are separate."""
    effective = {"schema_version": SCOPE_SCHEMA_VERSION, "scope_id": scope_id, "seed": SEED,
                 "params": copy.deepcopy(config["params"]), "refs": copy.deepcopy(config["refs"]),
                 "images": {key: config["deps"][key] for key in IMAGE_SETTINGS}}
    return effective


def load_assembly_provenance(library_id: str, sidecar: Path | None,
                             outputs: Mapping[str, Path]) -> dict[str, Any]:
    """Preserve recorded producer history and reject a sidecar for other bytes."""
    identities = {name: file_identity(path) for name, path in sorted(outputs.items())}
    if sidecar is None or not Path(sidecar).is_file():
        return {"library_id": library_id, "status": "unknown",
                "reason": "No assembly producer record was supplied for these reused outputs.",
                "outputs": identities}
    record = json.loads(Path(sidecar).read_text())
    if (record.get("schema_version") != ASSEMBLY_SCHEMA_VERSION or record.get("stage") != "assembly"
            or record.get("library_id") != library_id or not record.get("producer")):
        raise ValueError(f"Invalid assembly producer record for {library_id!r}")
    if set(record.get("outputs", {})) != set(outputs):
        raise ValueError(f"Assembly producer record names different outputs for {library_id!r}")
    for name, identity in identities.items():
        if record["outputs"][name] != identity:
            raise ValueError(f"Assembly producer record does not match {library_id!r} output {name!r}")
    return {"library_id": library_id, "status": "recorded", "record": record}


def compact_source(source: Mapping[str, Any]) -> dict[str, Any]:
    """Keep the revision and dirty flag, and the source-state hash of a dirty source.

    For the scope source, the hash names the directory write_source_evidence()
    fills with the patches and untracked files, so a dirty run stays
    identifiable beyond its revision.
    """
    result = {key: source[key] for key in ("revision", "dirty")}
    if source["dirty"]:
        result["state_sha256"] = source["state_sha256"]
    return result


def write_source_evidence(repo: Path, directory: Path, source: Mapping[str, Any]) -> None:
    """Retain the patches and untracked code needed to identify a dirty run."""
    if capture_source(repo) != source:
        raise ValueError("Source changed after workflow construction; rebuild the DAG")
    destination = directory / source["state_sha256"]
    destination.mkdir(parents=True, exist_ok=True)
    (destination / "source.json").write_text(json.dumps(source, indent=2, sort_keys=True) + "\n")
    (destination / "working.patch").write_bytes(_git(repo, "diff", "--binary", "HEAD", "--", *SOURCE_PATHS))
    (destination / "index.patch").write_bytes(_git(repo, "diff", "--cached", "--binary", "--", *SOURCE_PATHS))
    for entry in source["relevant_changes"]:
        if entry["status"] == "untracked":
            original = repo / entry["path"]
            target = destination / "untracked" / entry["path"]
            target.parent.mkdir(parents=True, exist_ok=True)
            if original.is_symlink():
                if not target.is_symlink():
                    target.symlink_to(original.readlink())
            else:
                shutil.copy2(original, target)


def build_scope_provenance(scope_id: str, libraries: Sequence[Mapping[str, str]],
                           config: Mapping[str, Any], *, source: Mapping[str, Any],
                           images: Mapping[str, Mapping[str, Any]],
                           tool_versions: Mapping[str, Mapping[str, Any]],
                           references: Mapping[str, Path], reference_receipt: Mapping[str, Any],
                           assembly_records: Mapping[str, Mapping[str, Any]] | None = None) -> dict[str, Any]:
    """Identify a scope while its manifest, configuration and calibration own values.

    The record holds:

    - `schema_version` (SCOPE_SCHEMA_VERSION) and `scope_id`;
    - `source`, the compacted pipeline source identity (see compact_source());
    - `images`, each image's format and digest, keyed by image_cache_key();
    - `software`, the analysis and configured assembly images with each tool's
      name, version and status, where a failed query stays "unknown" with no
      version;
    - `references`, the name and SHA-256 of each reference file and the image
      paths of the IgBLAST databases and germline directory;
    - `files`, the release-relative paths of the manifest, settings and
      threshold tables, which hold the values this record does not repeat;
    - `assembly_producers`, one producer identity per library with recorded
      history, and the libraries whose reused outputs have no producer record,
      which stay explicitly unknown.

    reference_receipt is the record audit_imgt_reference.py wrote for one of
    the references. Its IMGT download date, agreement flags, auxiliary-data
    pair checks and added entries, and the J alleles without auxiliary data are
    copied onto that reference's entry, and each assembly producer's reference
    is compared with it. A producer keeps its compacted source, image, measured
    tools, reference identity and actual assembly settings. Detailed reference
    sequence classes remain in the shared reference receipt.
    """
    crosswalk = [{key: row[key] for key in ("library_id", "tissue_id", "donor_id")} for row in libraries]
    if (not crosswalk or any(not isinstance(value, str) or not value for row in crosswalk for value in row.values())
            or len({row["library_id"] for row in crosswalk}) != len(crosswalk)):
        raise ValueError("Scope provenance requires one row per string library ID")
    if set(images) != {"analysis", "assembly"} or set(tool_versions) != set(images):
        raise ValueError("Both analysis and configured assembly image/version records are required")
    image_catalog = {}

    def image_id(image):
        key = image_cache_key(image)
        image_catalog[key] = {field: image[field] for field in ("format", "digest")}
        return key

    def measured_tools(entries):
        result = []
        for entry in entries:
            if (entry.get("status") not in {"measured", "unknown"}
                    or (entry["status"] == "measured" and not entry.get("version"))
                    or (entry["status"] == "unknown" and entry.get("version") is not None)):
                raise ValueError("Software status and measured version disagree")
            result.append({field: entry[field] for field in ("name", "version", "status")})
        return result

    software = {}
    for role, image in images.items():
        measured = tool_versions[role]
        if image_cache_key(measured["image"]) != image_cache_key(image):
            raise ValueError(f"Measured tool versions do not belong to the {role} image")
        software["configured_assembly" if role == "assembly" else role] = {
            "image": image_id(image), "tools": measured_tools(measured["software"])}
    resource_ids = {name: {key: value for key, value in file_identity(path).items()
                           if key in ("name", "sha256")} for name, path in sorted(references.items())}
    audited = [name for name, identity in resource_ids.items()
               if identity["sha256"] == reference_receipt["reference"]["sha256"]]
    if len(audited) != 1 or image_cache_key(reference_receipt["image"]) != software["analysis"]["image"]:
        raise ValueError("The reference receipt must describe one supplied reference and the analysis image")
    if not all(reference_receipt["agreement"].values()):
        raise ValueError("The reference receipt must record agreement; the audit writes no receipt otherwise")
    annotation_sha256 = resource_ids[audited[0]]["sha256"]
    auxiliary = reference_receipt["auxiliary_data"]
    resource_ids[audited[0]].update(
        imgt_snapshot={key: reference_receipt["imgt_snapshot"][key] for key in ("source", "date")},
        audit={"agreement": dict(reference_receipt["agreement"]),
               "auxiliary_data": {
                   "identical_pairs_checked": auxiliary["identical_pairs_checked"],
                   "sibling_pairs_checked": auxiliary["sibling_pairs_checked"],
                   "supplemented": [{"allele": entry["allele"], "rule": entry["rule"],
                                     "copied_from": list(entry["copied_from"])}
                                    for entry in auxiliary["supplemented"]],
                   "j_alleles_without_entry": [{"allele": entry["allele"], "functionality": entry["functionality"]}
                                               for entry in auxiliary["j_alleles_without_entry"]]}})
    if "reference_metadata" in references:
        metadata = reference_receipt["reference_metadata"]
        if resource_ids["reference_metadata"]["sha256"] != metadata["sha256"]:
            raise ValueError("Reference metadata differs from the audited reference snapshot")
    for name in ("igblast_db", "germline_dir"):
        if name in config["params"].get("igblast", {}):
            resource_ids[name] = {"image": software["analysis"]["image"],
                                  "path": config["params"]["igblast"][name]}
    assembly_records = assembly_records or {}
    known, producers, unknown = {}, {}, []
    for row in crosswalk:
        library = row["library_id"]
        entry = assembly_records.get(library) or load_assembly_provenance(library, None, {})
        if entry.get("library_id") != library or entry.get("status") not in {"recorded", "unknown"}:
            raise ValueError(f"Assembly producer record belongs to another library: {library!r}")
        if entry["status"] == "unknown":
            unknown.append(library)
            continue
        original = entry.get("record", {})
        if original.get("library_id") != library or original.get("schema_version") != ASSEMBLY_SCHEMA_VERSION:
            raise ValueError(f"Nested assembly producer record belongs to another library: {library!r}")
        producer = original["producer"]
        reference = original["inputs"]["reference"]
        identity = _hash(_json({"producer": producer, "reference": reference}).encode())
        known[library] = identity
        producers[identity] = {
            "source": compact_source(producer["source"]), "image": image_id(producer["image"]),
            "software": measured_tools(producer.get("software", [])),
            "settings": copy.deepcopy(producer["settings"]),
            # A reused assembly can predate the annotation reference. The flag
            # says whether TRUST4 read the same bytes IgBLAST annotation did.
            "reference": {key: reference[key] for key in ("name", "sha256")}
            | {"same_as_annotation_reference": reference["sha256"] == annotation_sha256}}
    return {
        "schema_version": SCOPE_SCHEMA_VERSION, "scope_id": scope_id,
        "source": compact_source(source), "images": image_catalog, "software": software,
        "references": resource_ids,
        "files": {"manifest": "manifest/scope_manifest.tsv", "effective_config": "effective_config.yaml",
                  "threshold": f"clones/{scope_id}_clone_threshold.tsv",
                  "donor_thresholds": f"clones/{scope_id}_donor_clone_thresholds.tsv",
                  "distance_diagnostics": f"clones/{scope_id}_clone_distance_diagnostics.tsv"},
        "assembly_producers": {"libraries": known, "producers": producers,
                               "unknown": {"libraries": unknown,
                                           "reason": "No assembly producer record was supplied for these reused outputs."}}}


def write_scope_provenance(record: Mapping[str, Any], effective_config: Mapping[str, Any],
                           provenance_path: Path, effective_config_path: Path) -> None:
    """Write compact identity and authoritative settings as separate atomic files."""
    import yaml

    for path, text in (
        (Path(provenance_path), json.dumps(record, indent=2, sort_keys=True, allow_nan=False) + "\n"),
        (Path(effective_config_path), yaml.safe_dump(dict(effective_config), sort_keys=False, allow_unicode=False)),
    ):
        path.parent.mkdir(parents=True, exist_ok=True)
        temporary = path.with_suffix(".tmp")
        temporary.write_text(text)
        temporary.replace(path)
