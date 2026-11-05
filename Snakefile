# Validate run inputs, select scope membership and load the workflow rules.
# Library jobs are shared across scopes. Donor-pooled IGH jobs run separately
# within each scope, and the final targets request each selected release receipt.

import pandas as pd
import re
from pathlib import Path

from snakemake.utils import min_version, validate

min_version("9.27.1")

# Configured paths use the launch directory. Bundled code, schemas, report assets
# and environments use this Snakefile's directory, including in a parent workflow.
WORKFLOW_DIR = workflow.basedir
SCHEMA_DIR = Path(WORKFLOW_DIR) / "workflow" / "schemas"

# Validate configuration before reading its fields. The schema supplies defaults
# and rejects unknown keys, invalid values and paths that TRUST4 cannot accept.
# File contents and membership relationships are checked below.
validate(config, str(SCHEMA_DIR / "config.schema.yaml"))
ANALYSIS_ONLY = config["analysis_only"]

# Read one samplesheet row per sequencing library. Keep identifiers as strings:
# leading zeroes and literal NA values are labels. Strip surrounding whitespace
# before samplesheet.schema.yaml checks the columns and path-safe identifiers.
# The checks here also reject an empty sheet, repeated libraries and tissues
# assigned to more than one donor. They run before scope selection.
REQUIRED_COLUMNS = ("library_id", "tissue_id", "donor_id", "path")

samplesheet = pd.read_csv(
    config['deps']['samplesheet'], dtype="string", keep_default_na=False)
if samplesheet.empty:
    raise ValueError("samplesheet must contain at least one library")

_missing = [c for c in REQUIRED_COLUMNS if c not in samplesheet.columns]
if _missing:
    raise ValueError(f"deps.samplesheet must contain column(s) {_missing}. Add them. "
                     "docs/running.md lists the required samplesheet columns.")

for _col in REQUIRED_COLUMNS:
    samplesheet[_col] = samplesheet[_col].str.strip()

validate(samplesheet, str(SCHEMA_DIR / "samplesheet.schema.yaml"), set_default=False)

_duplicate_libraries = sorted(set(samplesheet.loc[samplesheet['library_id'].duplicated(), 'library_id'].tolist()))
if _duplicate_libraries:
    raise ValueError(f"samplesheet: library_id must be unique; duplicates: {_duplicate_libraries}")

# A biological specimen has one donor even when several libraries describe it.
_tissue_donors = samplesheet.groupby("tissue_id")["donor_id"].nunique()
if (_tissue_donors > 1).any():
    raise ValueError("samplesheet: each tissue_id must map to one donor_id")

# Library identifiers drive scope-invariant targets; donor and tissue metadata
# remain attached through the samplesheet and scope manifests.
LIBRARY_IDS = samplesheet['library_id'].tolist()
LOCI = ['IGH', 'IGK', 'IGL', 'TRA', 'TRB', 'TRD', 'TRG']

# A within-donor overlap comparison needs at least two libraries in that scope.
# scope.smk uses SCOPE_MULTI_PAIRS to request these comparisons.
MIN_OVERLAP_SAMPLES = 2

# Configured output paths remain relative to the launch directory. The
# resources directory is shared by every run and configuration that names it,
# so it resolves to one absolute path: Snakemake then keeps one record for
# each obtained image, receipt and reference file however a configuration
# spells the path.
results_dir = Path(config['outs']['results_dir'])
resources_dir = Path(config['deps']['resources']).resolve()
log_dir = Path(config['outs']['log_dir'])

# Read one membership decision per project and library. usage=fail_qc excludes
# the library from that project before donor pooling and release. The reserved
# cohort scope includes every samplesheet library, regardless of these decisions.
# An omitted or empty membership_csv setting gives a cohort-only run.
membership_csv = config.get('membership_csv') or ""
MEMBERSHIP_COLUMNS = ("project", "library_id", "usage")
FAILED_QC = "fail_qc"
# Reject unknown usage values before filtering: a misspelt fail_qc would otherwise
# survive the exclusion and add the library to the project's analyses.
MEMBERSHIP_USAGE = ("analysed", FAILED_QC)
if membership_csv:
    _membership_path = Path(membership_csv)
    if not _membership_path.is_file():
        raise ValueError(
            f"membership_csv={membership_csv} is set but is not an existing file")
    _membership = pd.read_csv(_membership_path, dtype="string", keep_default_na=False)
    _absent = [c for c in MEMBERSHIP_COLUMNS if c not in _membership.columns]
    if _absent:
        raise ValueError(
            f"membership_csv={membership_csv} is missing required column(s): "
            f"{_absent}; found columns: {list(_membership.columns)}")
    for _col in MEMBERSHIP_COLUMNS:
        _membership[_col] = _membership[_col].str.strip()
    # Report record positions including the header, so the first data row is 2.
    _membership['_csv_line'] = _membership.index + 2
    _blank_project = _membership.loc[_membership['project'] == "", '_csv_line'].tolist()
    if _blank_project:
        raise ValueError(f"membership_csv={membership_csv} has a blank project on "
                         f"line(s) {_blank_project}")
    # Validate every decision before dropping fail_qc rows, including repeated
    # decisions for the same project and library. A fail_qc row may name a library
    # absent from the samplesheet; only analysed rows need matching library data.
    _unknown_usage = _membership.loc[
        ~_membership['usage'].isin(MEMBERSHIP_USAGE), '_csv_line'].tolist()
    if _unknown_usage:
        raise ValueError(
            f"membership_csv={membership_csv} has a usage outside "
            f"{list(MEMBERSHIP_USAGE)} on line(s) {_unknown_usage}")
    _rows_per_library = {}
    for _project, _library, _usage, _line in zip(
            _membership['project'], _membership['library_id'],
            _membership['usage'], _membership['_csv_line']):
        _rows_per_library.setdefault((_project, _library), []).append((_line, _usage))
    _repeated = sorted(pair for pair, seen in _rows_per_library.items() if len(seen) > 1)
    if _repeated:
        _detail = "; ".join(
            "{}/{} on line(s) {}".format(_project, _library, ", ".join(
                "{} ({})".format(_line, _usage)
                for _line, _usage in _rows_per_library[(_project, _library)]))
            for _project, _library in _repeated)
        raise ValueError(
            f"membership_csv={membership_csv} records more than one decision for a "
            f"library in a project: {_detail}")
    _analysed = _membership[_membership['usage'] != FAILED_QC]
    # Reject projects emptied by fail_qc rather than silently omitting their scope.
    _emptied = sorted(set(_membership['project']) - set(_analysed['project']))
    if _emptied:
        raise ValueError(
            f"membership_csv={membership_csv}: project(s) {_emptied} have no library "
            f"left after dropping usage='{FAILED_QC}' rows")
    if _analysed.empty:
        logger.warning(
            "membership_csv=%s assigns no libraries to a project; proceeding with the "
            "reserved 'cohort' scope only.", membership_csv)
    _project_groups = list(_analysed.groupby('project', sort=True))
else:
    _project_groups = []

# Build all declared memberships before selecting targets. Assembly, annotation
# and non-IGH analysis run once per library. IGH clone inference pools libraries
# by donor independently within each scope. A library may join several projects.
RESERVED_SCOPE = "cohort"
_SCOPE_NAME_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_-]*$")
# Scope names become result directories and clone-threshold filename prefixes.
# Keep them distinct from the cohort and the directory names used in a release.
_RESERVED_SCOPE_NAMES = {name.lower() for name in (
    RESERVED_SCOPE, "release", "intermediate", "libraries", "scopes", "shared",
    "metrics", "manifest", "clones", "report", "airr", "trust4")}

# scope_id -> set(library_id); cohort = every library in the master samplesheet.
SCOPE_MEMBERS = {RESERVED_SCOPE: set(LIBRARY_IDS)}
SCOPE_CSVS = {}  # Map each project scope to its source CSV path.
# Reject case-only differences so scopes stay distinct on case-insensitive disks.
_scope_lower = {RESERVED_SCOPE.lower(): RESERVED_SCOPE}
for _scope, _rows in _project_groups:
    # Check scope paths and library membership before building donor input lists.
    _lines = sorted(_rows['_csv_line'].tolist())
    if not _SCOPE_NAME_RE.match(_scope):
        raise ValueError(f"scope name '{_scope}' ({membership_csv} line(s) {_lines}) is "
                         f"not path-safe (letters/digits/_/-, must start alphanumeric)")
    if _scope.lower() in _RESERVED_SCOPE_NAMES:
        raise ValueError(f"scope name '{_scope}' ({membership_csv} line(s) {_lines}) is "
                         f"reserved: it names a directory of the results tree or the "
                         f"internal full-cohort scope")
    if _scope.lower() in _scope_lower:
        raise ValueError(f"scope name '{_scope}' case-fold collides with "
                         f"'{_scope_lower[_scope.lower()]}'")
    _scope_lower[_scope.lower()] = _scope
    _col = _rows['library_id']
    _blank = _rows.loc[_col == "", '_csv_line'].tolist()
    if _blank:
        raise ValueError(f"scope '{_scope}' ({membership_csv}) has a blank library_id on "
                         f"line(s) {_blank}")
    _unknown = sorted(set(_col) - SCOPE_MEMBERS[RESERVED_SCOPE])
    if _unknown:
        raise ValueError(f"scope '{_scope}' ({membership_csv}) lists library_ids absent from the "
                         f"master samplesheet: {_unknown}")
    SCOPE_MEMBERS[_scope] = set(_col)
    SCOPE_CSVS[_scope] = str(_membership_path)

# --scope narrows requested releases and donor inputs, while validation above
# still checks the complete samplesheet and membership file.
SCOPES = sorted(SCOPE_MEMBERS)
if config.get("scope"):
    if config["scope"] not in SCOPE_MEMBERS:
        raise ValueError(f"Unknown scope: {config['scope']}")
    SCOPES = [config["scope"]]

# Retain every membership decision for build_scope_assignment in inputs.smk.
# in_samplesheet marks whether the library is in the samplesheet. analysed marks
# inclusion in that project after fail_qc filtering. Downstream record QC is
# separate.
# Samplesheet libraries with no project row receive an empty-project row, and
# every samplesheet library still belongs to the cohort.
_claimed_libraries = set()
SCOPE_ASSIGNMENT_ROWS = []
if membership_csv:
    for _project, _library, _usage in zip(
            _membership['project'], _membership['library_id'], _membership['usage']):
        _claimed_libraries.add(_library)
        _in_samplesheet = _library in SCOPE_MEMBERS[RESERVED_SCOPE]
        SCOPE_ASSIGNMENT_ROWS.append(
            (_library, _project, _usage, _in_samplesheet,
             _in_samplesheet and _usage != FAILED_QC))
SCOPE_ASSIGNMENT_ROWS += [
    (_library, "", "", True, True)
    for _library in SCOPE_MEMBERS[RESERVED_SCOPE] - _claimed_libraries]

# Expand only observed donor pairs in selected scopes. Keep each scope's own
# rows in scope_manifest_by_scope, which inputs.smk tracks as manifest parameters.
# Editing unrelated membership therefore leaves this scope's manifest unchanged.
scope_donor_to_libraries = {}
scope_manifest_by_scope = {s: [] for s in SCOPES}  # scope -> [(library_id, donor_id, tissue_id), ...]
for _scope in SCOPES:
    _sub = samplesheet[samplesheet['library_id'].isin(SCOPE_MEMBERS[_scope])]
    for _donor, _grp in _sub.groupby('donor_id'):
        # Sort independently of CSV row order so pooled inputs and overlap output
        # are deterministic.
        scope_donor_to_libraries[(_scope, str(_donor))] = sorted(
            (str(r.tissue_id), str(r.library_id)) for r in _grp.itertuples(index=False)
        )
    for r in _sub.itertuples(index=False):
        scope_manifest_by_scope[_scope].append((str(r.library_id), str(r.donor_id), str(r.tissue_id)))

SCOPE_DONOR_PAIRS  = sorted(scope_donor_to_libraries)
SCOPE_MULTI_PAIRS  = [(s, d) for (s, d) in SCOPE_DONOR_PAIRS
                      if len(scope_donor_to_libraries[(s, d)]) >= MIN_OVERLAP_SAMPLES]

def _scope_targets(subpath):
    """Request one output path per selected scope."""
    return expand(str(results_dir / ("scopes/{scope}/" + subpath)), scope=SCOPES)

_unknown_metadata_scopes = set(config["scope_metadata"]) - set(SCOPE_MEMBERS)
if _unknown_metadata_scopes:
    raise ValueError(f"scope_metadata names undeclared scopes: {sorted(_unknown_metadata_scopes)}")

# Load shared paths, reference preparation and scope manifest writers first.
# analysis_only omits assembly.smk, including FASTQ discovery and depth counting.
# Downstream rules still require retained _airr.tsv, _cdr3.out, _final.out and
# _seq_depth.tsv files for selected libraries. Missing inputs fail during DAG
# creation. Image and reference preparation remain available through inputs.smk.
include: "workflow/rules/common.smk"
include: "workflow/rules/inputs.smk"
if not ANALYSIS_ONLY:
    include: "workflow/rules/assembly.smk"
include: "workflow/rules/library.smk"
include: "workflow/rules/igh.smk"
include: "workflow/rules/scope.smk"
include: "workflow/rules/provenance.smk"
include: "workflow/rules/exports.smk"
include: "workflow/rules/report.smk"

# Both targets request every selected scope's release receipt. run_module.sh
# --post-trust4 names post_trust4_all. The analysis_only setting, not the target,
# is what leaves FASTQ discovery, depth counting and assembly out of the DAG.
rule all:
    input:
        *_scope_targets("release/.complete.json"),
        *_scope_targets("intermediate/README.md"),
    message:
        "The airrborne workflow successfully completed for all samples."


rule post_trust4_all:
    input:
        *_scope_targets("release/.complete.json"),
        *_scope_targets("intermediate/README.md"),
    message:
        "Post-TRUST4 workflow successfully completed for all samples."
