"""Package one complete scope handover with an explicit file allow-list.

Stage the whole folder in a temporary directory beside its destination, validate
those staged files, list them in manifest/release_files.tsv, and only then
replace the published folder. The replacement is two atomic renames, and a small
transaction record restores the previous folder if publication stops between
them. One workflow rule owns the released folder and declares only the receipt
.complete.json, which its script writes after publication. Source tables must
also be replaced atomically: rewriting a hard-linked source in place would alter
an already released table.
"""

import copy
import csv
import gzip
import hashlib
import json
import os
from pathlib import Path, PurePosixPath
import re
import shutil
import stat
import tempfile
from typing import Any, Callable, Mapping, Sequence


# The scope holding every library, as the Snakefile names it. Only that bundle
# carries the run-wide table below; the Snakefile asserts the two names agree.
RESERVED_COHORT_SCOPE = "cohort"

# Pipeline output files the data-dictionary contract deliberately skips, with the
# reason for each. They are TRUST4's and Change-O's own output, kept for diagnosis
# rather than tables this pipeline defines, so their schemas belong to those tools.
# None of them is part of a release bundle.
EXEMPT_SUFFIXES = {
    "_makedb_pass.tsv.gz": "raw MakeDb rows before CIGAR normalisation, kept for diagnosis",
    "_makedb_fail.tsv": "records MakeDb could not convert, kept for diagnosis",
    "_report.tsv": "TRUST4's own report, retained because cleanup is disabled",
    "_airr_align.tsv": "TRUST4's own alignment table, retained with the report",
}

# Files whose columns are under contract but whose values are not yet. The
# reannotated table is written before the common QC gate and again after it, in
# the same schema, but several of the declared vocabularies are what that gate
# establishes: productive is documented as T because every retained row is
# productive, which is true of the filtered table and not of its input. Checking
# the pre-filter copy against them would mean widening the vocabulary and losing
# the assertion where it means something.
HEADER_ONLY_SUFFIXES = {
    "_igblast.tsv": "the reannotated table before the common QC gate, so the "
                    "vocabularies that gate establishes do not hold in it yet",
}

# What a project receives, relative to its scope directory. `{scope}` is
# substituted with the scope name. Everything else a scope directory holds stays
# behind: the donor work, including the complete member shards; the per-donor
# clone membership files, whose columns clone_members.tsv carries, with v_gene,
# j_gene and consensus_count as v_gene_key, j_gene_key and
# member_consensus_count; and the tables SCOPE_WORKING_TABLES names.
SCOPE_RELEASE_FILES = [
    "manifest/scope_manifest.tsv",
    # Collected library tables and scope audits describe this scope's libraries.
    "metrics/yield.tsv",
    "metrics/qc.tsv",
    "metrics/qc_trail.tsv",
    "metrics/diversity.tsv",
    "metrics/gene_usage.tsv",
    "metrics/vj_pair_usage.tsv",
    "metrics/clonotypes.tsv",
    "metrics/clonotype_members.tsv",
    "metrics/isotype_usage.tsv",
    "metrics/rarefaction_curves.tsv",
    "metrics/clone_shm.tsv",
    # Scope-native IGH tables, computed over this scope's libraries.
    "metrics/shm.tsv",
    "metrics/selection_pressure.tsv",
    "metrics/clones.tsv",
    "metrics/clone_isotype.tsv",
    "metrics/clone_members.tsv",
    "metrics/clone_member_shm.tsv",
    "metrics/clone_member_shm_coverage.tsv",
    # Both arms: the clone-size vector, and the clones a donor's libraries share.
    "metrics/group_abundance.tsv",
    "metrics/donor_overlap_pairs.tsv",
    "metrics/donor_overlap_summary.tsv",
    "metrics/shared_groups.tsv",
    # How the IGH lineage cutoff was chosen, so a reader can judge the unit the
    # IGH tables count. release_path() places these in igh_clone_threshold/.
    "clones/{scope}_clone_threshold.tsv",
    "clones/{scope}_donor_clone_thresholds.tsv",
    "clones/{scope}_clone_distance_diagnostics.tsv",
    "metrics/headline.tsv",
    "provenance.json",
    "effective_config.yaml",
    "airr/Repertoire.json",
    "airr/validation.json",
    "report",
]


def scope_release_files(scope: str, membership: bool) -> list[str]:
    """The bundle's contents for one scope, as scope-relative paths.

    One file is not a property of a scope at all: scope_assignment.tsv records what
    the run did with every library, including the ones no project claimed and the
    ones QC rejected. It is written once, beside the reserved cohort scope's
    manifest, and ships in that bundle only. A project receives the libraries it
    was given, not the register of everybody else's. Without a membership file
    the register has no projects to record, so `membership` false leaves it out.
    """
    files = [path.format(scope=scope) for path in SCOPE_RELEASE_FILES]
    if scope == RESERVED_COHORT_SCOPE and membership:
        files.append("manifest/scope_assignment.tsv")
    return files


def scope_release_source_files(scope: str, membership: bool) -> list[str]:
    """Producer inputs exclude validation generated by the publisher."""
    return [path for path in scope_release_files(scope, membership)
            if path != "airr/validation.json"]


# The dictionary table of a scope-level metrics file is its stem. The allow-list
# is the one owner of that set, so a table added there is recognised here too.
SCOPE_METRIC_TABLES = frozenset(
    PurePosixPath(entry).stem for entry in SCOPE_RELEASE_FILES
    if entry.startswith("metrics/") and entry.endswith(".tsv"))

# Scope metric tables the release leaves out. The headline reads the callable
# totals and IgG:IgA ratio of isotype_summary.tsv, and the headline and
# isotype_usage.tsv carry its other values. The collected qc_audit.tsv stays in
# the working tree, and qc_trail.tsv is its released view.
SCOPE_WORKING_TABLES = frozenset({"isotype_summary", "qc_audit"})

# Tables with one row per clone, member, clonotype or shared observation, and
# the per-stage QC trail. The release places them in metrics/detail/, apart
# from the per-library summaries, while the scope working tree keeps every table
# in metrics/.
RELEASE_DETAIL_TABLES = frozenset({
    "clone_members", "clone_member_shm", "clonotype_members", "clonotypes", "clones",
    "clone_shm", "clone_isotype", "group_abundance", "shared_groups", "qc_trail", "clone_member_shm_coverage"})
assert RELEASE_DETAIL_TABLES <= SCOPE_METRIC_TABLES


# The IGH threshold files. The scope working tree keeps them as
# clones/<scope><suffix>, the paths provenance.py records, and the release places
# each one in igh_clone_threshold/ under its dictionary table's name.
OUTPUT_SCHEMA_VERSION = "2.0.0"
IGH_DETAIL_TABLES = frozenset({"clones", "clone_members", "clone_member_shm",
                               "clone_shm", "clone_isotype", "clone_member_shm_coverage"})
NONIGH_DETAIL_TABLES = frozenset({"clonotypes", "clonotype_members"})

THRESHOLD_RELEASE_DIRECTORY = "igh_clone_threshold"
THRESHOLD_TABLES = {
    "_clone_threshold.tsv": "scope_threshold",
    "_donor_clone_thresholds.tsv": "donor_fits",
    "_clone_distance_diagnostics.tsv": "donor_distance_histograms",
}


def release_path(relative: str) -> str:
    """The published path of a scope-relative allow-list entry."""
    path = PurePosixPath(relative)
    if path.parent == PurePosixPath("metrics") and path.stem in RELEASE_DETAIL_TABLES:
        directory = "metrics/detail/igh" if path.stem in IGH_DETAIL_TABLES else (
            "metrics/detail/non_igh" if path.stem in NONIGH_DETAIL_TABLES else "metrics/detail")
        return str(PurePosixPath(directory) / path.name)
    if path.parent == PurePosixPath("clones"):
        for suffix, table in THRESHOLD_TABLES.items():
            if path.name.endswith(suffix):
                return f"{THRESHOLD_RELEASE_DIRECTORY}/{table}.tsv"
    return relative


def release_provenance(record: Mapping[str, Any]) -> dict[str, Any]:
    """Return a scope provenance record whose files map names release paths.

    provenance.py records the scope working-tree path of each file in `files`.
    The release and its report carry this copy, in which release_path() gives
    each path, so both name the files a recipient finds in the folder.
    """
    released = copy.deepcopy(dict(record))
    released["output_schema_version"] = OUTPUT_SCHEMA_VERSION
    released["files"] = {name: release_path(path) for name, path in record["files"].items()}
    return released


def _library_file(path: Path) -> bool:
    """True for libraries/<id>/<subdir>/<id>_<suffix>, judged by position.

    A per-library file is named after its library, so a library called clone
    writes libraries/clone/metrics/clone_isotype.tsv, which shares its name with
    the scope table metrics/clone_isotype.tsv. Position separates them without
    depending on the directory a copied release happens to sit in.
    """
    parts = path.parts
    return (len(parts) >= 4 and parts[-4] == "libraries"
            and path.name.startswith(parts[-3] + "_"))


# The publisher's list of the release's files, with its dictionary table.
RELEASE_FILES_MANIFEST = "manifest/release_files.tsv"
RELEASE_FILES_TABLE = "release_files"


def table_for(path: Path) -> str | None:
    """The dictionary table a released file belongs to, or None if unrecognised."""
    name = path.name
    library = _library_file(path)
    in_detail = ((path.parent.name == "detail" and path.parent.parent.name == "metrics") or
                 (path.parent.name in ("igh", "non_igh") and path.parent.parent.name == "detail"
                  and path.parent.parent.parent.name == "metrics"))
    if (not library and path.parent.name == "metrics" and name.endswith(".tsv")
            and path.stem in SCOPE_METRIC_TABLES | SCOPE_WORKING_TABLES):
        return path.stem
    if not library and in_detail and name.endswith(".tsv") and path.stem in RELEASE_DETAIL_TABLES:
        return path.stem
    if name.endswith("_rearrangements.tsv.gz"):
        return "exported_AIRR"
    # Audit shards are the one family matched inside the name rather than by
    # suffix, so their directory is checked too. A scope or library whose own name
    # contained _qc_audit_ would otherwise be filed as an audit table.
    if ("_qc_audit_" in name and name.endswith(".tsv")
            and (path.parent.name == "qc" or "intermediate" in path.parts)):
        return "qc_audit"
    if name.endswith("_nonigh_grouping.tsv"):
        return "nonigh_grouping"
    if name.endswith(("_igh_grouping.tsv", "_IGH_grouping.tsv")):
        return "igh_grouping"
    if name.endswith("_seq_depth.tsv"):
        return "native_input_depth"
    if name.endswith("_input_depth.tsv"):
        return "input_depth"
    library_metrics = {
        "_rarefaction_curves.tsv": "rarefaction_curves",
        "_diversity.tsv": "diversity",
        "_isotype_summary.tsv": "isotype_summary",
        "_group_abundance.tsv": "group_abundance",
        "_isotype_usage.tsv": "isotype_usage",
        "_yield.tsv": "yield",
        "_clonotype_gene_usage.tsv": "gene_usage",
        "_clonotype_vj_pair_usage.tsv": "vj_pair_usage",
        "_clonotypes.tsv": "clonotypes",
        "_clonotype_members.tsv": "clonotype_members",
    }
    if library:
        for suffix, table in library_metrics.items():
            if name.endswith(suffix):
                return table
    if name.endswith("_airr_qc_pass.tsv"):
        return "reannotated_AIRR"
    # The pre-QC reannotation carries the same schema as the QC-pass subset.
    if name.endswith("_igblast.tsv"):
        return "reannotated_AIRR"
    if name.endswith("_igh_consensus_votes.tsv"):
        return "IGH_consensus_votes"
    if name.endswith("_isotype_calls.tsv"):
        return "IGH_isotype_calls"
    if name.endswith("_igh_features.tsv"):
        return "IGH_features"
    if name.endswith("_airr.tsv"):
        return "native_TRUST4_AIRR"
    if name.endswith("_IGH_clone_membership.tsv"):
        return "IGH_clone_membership"
    if name == "scope_manifest.tsv":
        return "scope_manifest"
    if name == "scope_assignment.tsv":
        return "scope_assignment"
    if name == "release_files.tsv":
        return RELEASE_FILES_TABLE
    if (path.parent.name == THRESHOLD_RELEASE_DIRECTORY and name.endswith(".tsv")
            and path.stem in THRESHOLD_TABLES.values()):
        return path.stem
    for suffix, table in THRESHOLD_TABLES.items():
        if name.endswith(suffix):
            return table
    return None


def is_exempt(path: Path) -> bool:
    return any(path.name.endswith(suffix) for suffix in EXEMPT_SUFFIXES)


def is_header_only(path: Path) -> bool:
    return any(path.name.endswith(suffix) for suffix in HEADER_ONLY_SUFFIXES)


GENERATED_FILES = {"README.md", "data_dictionary.csv", RELEASE_FILES_MANIFEST,
                   "provenance.json", "effective_config.yaml", "airr/validation.json"}


def _relative(value: str) -> Path:
    path = PurePosixPath(str(value))
    if (not str(value) or path.is_absolute() or any(part in ("", ".", "..") for part in str(value).split("/"))
            or "\\" in str(value) or any(character in str(value) for character in "\t\r\n")):
        raise ValueError(f"Release destination must be a plain relative path: {value!r}")
    return Path(*path.parts)


def _scope_source(scope_dir: Path, relative: str) -> Path:
    path = _relative(relative)
    source = scope_dir
    for part in path.parts:
        source = source / part
        if source.is_symlink():
            raise ValueError(f"Scope release sources must not be symlinks: {source}")
    return source


def _place(source: Path, destination: Path, file_mode: int) -> None:
    """Link a regular source when mode and group match; otherwise copy it.

    Copying permits destination permissions without changing a source inode.
    Producers must replace linked source files atomically on later runs.
    """
    if source.is_symlink() or not source.is_file():
        raise ValueError(f"Release input must be a regular file: {source}")
    destination.parent.mkdir(parents=True, exist_ok=True)
    source_stat = source.stat()
    if (stat.S_IMODE(source_stat.st_mode) == file_mode
            and source_stat.st_gid == destination.parent.stat().st_gid):
        try:
            os.link(source, destination)
            return
        except OSError:
            pass
    shutil.copyfile(source, destination)
    destination.chmod(file_mode)


def _manifest(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    keys = ("scope_id", "library_id", "donor_id", "tissue_id")
    if not rows or any(not isinstance(row.get(key), str) or not row[key] for row in rows for key in keys):
        raise ValueError("Release requires a nonempty scope manifest with string identities")
    if len({row["scope_id"] for row in rows}) != 1 or len({row["library_id"] for row in rows}) != len(rows):
        raise ValueError("Release manifest must contain one scope and unique libraries")
    return rows


def _check_table(path: Path, manifest: list[dict], library_id: str | None = None) -> None:
    """Check shipped identity columns, streaming large compressed AIRR tables."""
    opener = gzip.open if path.name.endswith(".gz") else open
    allowed = {key: {row[key] for row in manifest} for key in ("scope_id", "library_id", "donor_id", "tissue_id")}
    with opener(path, "rt", newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        header = reader.fieldnames or []
        if not header or len(header) != len(set(header)):
            raise ValueError(f"{path.name}: missing or repeated columns")
        if library_id is not None and "library_id" not in header:
            raise ValueError(f"{path.name}: AIRR source requires library_id")
        for row in reader:
            if None in row or None in row.values():
                raise ValueError(f"{path.name}: malformed row")
            for key, values in allowed.items():
                for column in (key, key + "_1", key + "_2"):
                    if column in row and row[column] not in values:
                        raise ValueError(f"{path.name}: out-of-scope {column} {row[column]!r}")
            if library_id is not None and row["library_id"] != library_id:
                raise ValueError(f"{path.name}: AIRR source belongs to another library")


# Dictionary fields for maintainers only. The release copy leaves them out.
MAINTAINER_DICTIONARY_FIELDS = ("implementation",)


def _dictionary(source: Path, tables: set[str], destination: Path) -> None:
    with source.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        fields = [field for field in reader.fieldnames or []
                  if field not in MAINTAINER_DICTIONARY_FIELDS]
        if not {"table", "column"}.issubset(fields):
            raise ValueError("Dictionary requires table and column fields")
        rows = [row for row in reader if row["table"] in tables and not row["column"].startswith("sort_")]
    missing = tables - {row["table"] for row in rows}
    if missing:
        raise ValueError(f"Released tables have no dictionary entries: {sorted(missing)}")
    keys = {(row["table"], row["column"]) for row in rows}
    if len(keys) != len(rows):
        raise ValueError("Dictionary contains repeated table/column entries")
    with destination.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n",
                                extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def _write_json(path: Path, value: Any) -> None:
    path.write_text(json.dumps(value, ensure_ascii=True, allow_nan=False, indent=2) + "\n", encoding="utf-8")


def _write_transaction(path: Path, value: Any) -> None:
    pending = path.with_suffix(".pending")
    _write_json(pending, value)
    os.replace(pending, path)


def _transaction_paths(bundle: Path) -> tuple[Path, Path]:
    return (bundle.with_name("." + bundle.name + ".transaction.json"),
            bundle.with_name("." + bundle.name + ".previous"))


def _remove(path: Path) -> None:
    if path.is_dir() and not path.is_symlink():
        shutil.rmtree(path)
    elif path.exists() or path.is_symlink():
        path.unlink()


def recover_release_bundle(bundle: Path) -> None:
    """Restore an interrupted publication, or finish cleaning a published folder.

    The replacement moves the previous folder aside before the staged folder
    takes its place, so an interruption can leave the release absent. The
    transaction record says whether a previous folder existed and whether the
    replacement completed, which is enough to restore either state.
    """
    bundle = Path(bundle)
    marker, previous = _transaction_paths(bundle)
    if not marker.exists():
        if previous.exists():
            raise ValueError("A previous release exists without its publication record")
        return
    record = json.loads(marker.read_text())
    if record["phase"] not in {"publishing", "published"}:
        raise ValueError("Unknown release publication phase")
    if record["phase"] == "publishing":
        if previous.exists():
            _remove(bundle)
            os.replace(previous, bundle)
        elif not record["had_bundle"]:
            _remove(bundle)
    else:
        _remove(previous)
    marker.unlink()


def _publish(staged: Path, bundle: Path) -> None:
    marker, previous = _transaction_paths(bundle)
    record = {"phase": "publishing", "had_bundle": bundle.exists()}
    _write_transaction(marker, record)
    try:
        if record["had_bundle"]:
            os.replace(bundle, previous)
        os.replace(staged, bundle)
        record["phase"] = "published"
        _write_transaction(marker, record)
    except BaseException:
        recover_release_bundle(bundle)
        raise
    recover_release_bundle(bundle)


def _release_inputs(scope_dir: Path, bundle: Path, manifest: list[dict], *,
                    scope_files: Sequence[str], scope_directories: Sequence[str],
                    airr_files: Mapping[str, Mapping[Path, str]],
                    table_names: Mapping[str, str]) -> dict[Path, Path]:
    """Collect allowed files and validate their identities before staging.

    Reject symlinks, destination collisions and inputs inside the destination.
    AIRR tables must match their named library. Scope tables must contain only
    manifest identities, except the cohort's run-wide scope assignment table.
    """
    scope_id = manifest[0]["scope_id"]
    inputs = {}

    def add(source: Path, relative: str) -> None:
        destination = _relative(relative)
        if str(destination) in GENERATED_FILES or destination in inputs:
            raise ValueError(f"Repeated or reserved release destination: {relative}")
        if Path(source).is_symlink() or not Path(source).is_file():
            raise ValueError(f"Release input must be a regular file: {source}")
        inputs[destination] = Path(source)

    for relative in scope_files:
        add(_scope_source(scope_dir, relative), release_path(relative))
    if Path("manifest/scope_manifest.tsv") not in inputs:
        raise ValueError("The release allow-list must include its scope manifest")
    for relative in scope_directories:
        directory = _scope_source(scope_dir, relative)
        if not directory.is_dir():
            raise ValueError(f"Release directory is missing: {directory}")
        for source in sorted(directory.rglob("*")):
            if source.name == ".snakemake_timestamp":
                continue
            if source.is_symlink():
                raise ValueError(f"Release directory contains a symlink: {source}")
            if source.is_file():
                add(source, str(source.relative_to(scope_dir)))
    checked_airr = set()
    for library_id, sources in airr_files.items():
        if not sources:
            raise ValueError(f"No AIRR annotation supplied for {library_id!r}")
        for source, relative in sources.items():
            source = Path(source)
            if not source.name.endswith(".tsv.gz"):
                raise ValueError("Per-library AIRR annotation must be gzip-compressed TSV")
            add(source, relative)
            _check_table(source, manifest, library_id)
            checked_airr.add(source)
    if any(bundle == source.resolve() or bundle in source.resolve().parents for source in inputs.values()):
        raise ValueError("The release destination must not replace its own inputs")
    names = {str(path) for path in inputs}
    if set(table_names) - names:
        raise ValueError("Dictionary table mapping names files outside the release allow-list")
    for relative, source in inputs.items():
        if source.name.endswith((".tsv", ".tsv.gz")):
            if str(relative) not in table_names:
                raise ValueError(f"Released TSV has no dictionary table: {relative}")
            # The internal cohort register deliberately includes unanalysed
            # libraries. A project can never release that run-wide register.
            if str(relative) == "manifest/scope_assignment.tsv":
                if scope_id != "cohort":
                    raise ValueError("Only the internal cohort can release scope_assignment.tsv")
            elif source not in checked_airr:
                _check_table(source, manifest)
    return inputs


DOCUMENTATION_URL = "https://github.com/swarbricklab/airrborne/blob/{revision}/docs/{page}"
# Loci in the order the README lists rarefaction outcomes.
README_LOCI = ("IGH", "IGK", "IGL", "TRA", "TRB", "TRD", "TRG")


def _readme_context(staged: Path, manifest: list[dict], revision: str | None) -> dict[str, str]:
    """Sentences and links the README fills from this release.

    Counts donors with more than one library, and libraries whose rarefaction
    status is available at each locus with a configured target, from the staged
    diversity table when it is released. Describes manifest/scope_assignment.tsv
    only when the staged folder holds it. Documentation links name the source
    revision, or main when provenance records none.
    """
    libraries_by_donor: dict[str, int] = {}
    for row in manifest:
        libraries_by_donor[row["donor_id"]] = libraries_by_donor.get(row["donor_id"], 0) + 1
    repeated = sum(count > 1 for count in libraries_by_donor.values())
    if repeated:
        noun = "donor contributes" if repeated == 1 else "donors contribute"
        donors = (f"Libraries from one donor are not independent observations, and {repeated} "
                  f"{noun} more than one library.")
    else:
        donors = "Each donor contributes one library."
    diversity = staged / "metrics/diversity.tsv"
    rows = []
    if diversity.is_file():
        with diversity.open(newline="", encoding="utf-8") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
    outcomes = []
    for locus in README_LOCI:
        at_locus = [row for row in rows if row["locus"] == locus and row["rarefaction_target"]]
        if at_locus:
            reached = sum(row["rarefaction_status"] == "available" for row in at_locus)
            outcomes.append(f"{reached} of {len(at_locus)} libraries at {locus} "
                            f"({at_locus[0]['rarefaction_target']} support units)")
    listed = " and ".join(filter(None, (", ".join(outcomes[:-1]), *outcomes[-1:])))
    targets = (f"Rarefaction targets were reached by {listed}."
               if outcomes else "No rarefaction target is configured.")
    assignment = ""
    if (staged / "manifest/scope_assignment.tsv").is_file():
        assignment = ("This release also includes `manifest/scope_assignment.tsv` "
                      "(`scope_assignment`), keyed by `library_id` and\n`project`. A supplied "
                      "`fail_qc` decision excludes that project membership. Every\n"
                      "master-samplesheet library still enters `cohort`.\n")
    revision = revision or "main"
    return {"donor_repeats": donors, "target_summary": targets, "scope_assignment": assignment,
            "guide_url": DOCUMENTATION_URL.format(revision=revision, page="interpretation.md"),
            "running_url": DOCUMENTATION_URL.format(revision=revision, page="running.md"),
            "qc_matrix_url": DOCUMENTATION_URL.format(revision=revision, page="qc_matrix.md")}


def _write_readme(staged: Path, scope_id: str, manifest: list[dict],
                  table_names: Mapping[str, str], methods: Sequence[str],
                  readme_template: Path | None,
                  readme_values: Mapping[str, str] | None,
                  revision: str | None = None) -> None:
    """Write the release README from the template and validated scope context.

    Supplied values supplement the scope counts, run-specific sentences and
    generated methods. Stop before publication if any template field remains
    unresolved.
    """
    values = {
        **(readme_values or {}),
        **_readme_context(staged, manifest, revision),
        "scope": scope_id,
        "n_libraries": str(len(manifest)),
        "n_tissues": str(len({row["tissue_id"] for row in manifest})),
        "n_donors": str(len({row["donor_id"] for row in manifest})),
        "n_tables": str(len(set(table_names.values()))),
        "methods": "\n\n".join(methods),
    }
    readme = Path(readme_template).read_text() if readme_template else (
        "# Immune repertoire: {{scope}}\n\n"
        "This scope contains {{n_libraries}} libraries from {{n_tissues}} tissues and {{n_donors}} donors. "
        "The manifest identifies these observations; repeated libraries remain distinct.\n\n"
        "Open `report/index.html` in a browser and keep its folder together. "
        "Start an analysis with `metrics/headline.tsv`; the detailed tables retain their source measurements. "
        "`data_dictionary.csv` defines the {{n_tables}} shipped scientific tables. Optional project metadata appears only in the report.\n\n"
        "The per-library AIRR files are gzip-compressed TSVs. `provenance.json` records "
        "the source, tool and resource identities, and `effective_config.yaml` records "
        "the settings used.\n\n{{methods}}\n")
    for key, value in values.items():
        readme = readme.replace("{{" + key + "}}", str(value))
    if "{{" in readme:
        raise ValueError("Release README has unresolved template fields")
    (staged / "README.md").write_text(readme, encoding="utf-8")


def _file_sha256(path: Path) -> str:
    """Hash one file in 1 MiB reads, because released files reach hundreds of MB."""
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()


def _write_release_files(staged: Path) -> str:
    """Write RELEASE_FILES_MANIFEST over the staged folder and return its SHA256.

    The manifest has one row per file, sorted by release-relative path as text,
    and leaves out itself and the receipt .complete.json, which is written after
    publication. A staged file can be a hard link to a source table, so this step
    only reads the staged files and creates the manifest as a new file.
    """
    manifest = staged / RELEASE_FILES_MANIFEST
    relative = sorted(path.relative_to(staged).as_posix() for path in staged.rglob("*") if path.is_file())
    with manifest.open("x", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(("path", "n_bytes", "sha256"))
        for name in relative:
            writer.writerow((name, (staged / name).stat().st_size, _file_sha256(staged / name)))
    return _file_sha256(manifest)


def build_release_bundle(scope_dir: Path, bundle: Path, *,
                         scope_files: Sequence[str], scope_directories: Sequence[str] = (),
                         airr_files: Mapping[str, Mapping[Path, str]], dictionary: Path,
                         table_names: Mapping[str, str],
                         provenance: Path | Mapping[str, Any], effective_config: Path | str,
                         calibration: Path,
                         readme_template: Path | None = None,
                         readme_values: Mapping[str, str] | None = None,
                         validate_staged: Callable[[Path], None] | None = None) -> dict[str, int]:
    """Publish only named scope files, report directories and library exports.

    scope_files are scope-relative sources; release_path() gives each one's
    destination. table_names maps each shipped TSV destination to its
    dictionary table.
    airr_files maps each manifest library to its source/destination pairs and
    must cover the manifest exactly. Paths within the folder are relative, so a
    recipient can move the whole folder and keep every link working.

    The publisher also writes data_dictionary.csv, filtered to the shipped
    tables, provenance.json, effective_config.yaml and README.md, whose methods
    come from release_methods.methods_paragraphs(). validate_staged runs on the
    staged folder, and manifest/release_files.tsv is then written over it, before
    it replaces the release. Returns the counts of staged files, dictionary tables
    and libraries, and the manifest's SHA256. The file count includes the manifest.
    """
    scope_dir = Path(scope_dir).resolve()
    if Path(bundle).is_symlink():
        raise ValueError("The release destination must not be a symlink")
    bundle = Path(bundle).resolve()
    if bundle == scope_dir or bundle in scope_dir.parents:
        raise ValueError("The release destination must not replace its scope inputs")
    bundle.parent.mkdir(parents=True, exist_ok=True)
    recover_release_bundle(bundle)
    manifest = _manifest(_scope_source(scope_dir, "manifest/scope_manifest.tsv"))
    scope_id = manifest[0]["scope_id"]
    if set(airr_files) != {row["library_id"] for row in manifest}:
        raise ValueError("AIRR exports must cover exactly the in-scope libraries")
    inputs = _release_inputs(
        scope_dir, bundle, manifest, scope_files=scope_files,
        scope_directories=scope_directories, airr_files=airr_files,
        table_names=table_names)
    table_names = {**table_names, RELEASE_FILES_MANIFEST: RELEASE_FILES_TABLE}
    # Imported here rather than at module scope. The release contract checkers
    # import this module only for its allow-list and filename rules, on hosts
    # that need not carry PyYAML.
    import yaml
    from workflow.lib.release_methods import methods_paragraphs

    record = release_provenance(json.loads(Path(provenance).read_text())
                                if isinstance(provenance, (Path, str)) else provenance)
    if record["scope_id"] != scope_id:
        raise ValueError("Release provenance belongs to another scope")
    config_text = effective_config.read_text() if isinstance(effective_config, Path) else effective_config
    with calibration.open(newline="") as handle:
        thresholds = list(csv.DictReader(handle, delimiter="\t"))
    if len(thresholds) != 1:
        raise ValueError("Release methods require one operative threshold")
    methods = methods_paragraphs(record, yaml.safe_load(config_text), manifest, thresholds[0])
    staged = Path(tempfile.mkdtemp(prefix="." + bundle.name + ".stage-", dir=bundle.parent))
    publication_umask = os.umask(0)
    os.umask(publication_umask)
    file_mode = 0o666 & ~publication_umask
    directory_mode = 0o777 & ~publication_umask
    try:
        for relative, source in inputs.items():
            _place(source, staged / relative, file_mode)
        _dictionary(Path(dictionary), set(table_names.values()), staged / "data_dictionary.csv")
        _write_json(staged / "provenance.json", record)
        (staged / "effective_config.yaml").write_text(config_text, encoding="utf-8")
        _write_readme(staged, scope_id, manifest, table_names, methods,
                      readme_template, readme_values, record.get("source", {}).get("revision"))
        if validate_staged is not None:
            validate_staged(staged)
        release_files_sha256 = _write_release_files(staged)
        for path in [staged, *staged.rglob("*")]:
            target_mode = (directory_mode | (stat.S_IMODE(path.stat().st_mode) & stat.S_ISGID)
                           if path.is_dir() else file_mode)
            if stat.S_IMODE(path.stat().st_mode) != target_mode:
                path.chmod(target_mode)
        summary = {"output_schema_version": OUTPUT_SCHEMA_VERSION, "files": sum(path.is_file() for path in staged.rglob("*")),
                   "tables": len(set(table_names.values())), "libraries": len(manifest),
                   "release_files_sha256": release_files_sha256}
        _publish(staged, bundle)
        return summary
    finally:
        _remove(staged)
