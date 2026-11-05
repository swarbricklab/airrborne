# Locate shared code, images and reference files for the rules in Snakefile.
# Each rule declares the code and settings its job consumes. The processing-code
# list below also defines which source files identify each library AIRR export.
# Scientific decisions live in the scripts and modules those rules call.

import os
import csv
import json
import re
import shlex
import sys
from pathlib import Path

# scope_release_source_files() supplies the publisher's required input paths.
# release_bundle.py uses the same allow-list when validating the handover.
sys.path.insert(0, str(Path(WORKFLOW_DIR)))
from workflow.lib.provenance import image_receipt_path, image_sif_path
from workflow.lib.release_bundle import RESERVED_COHORT_SCOPE, scope_release_source_files

# The reserved scope is named twice: here it is a directory of the results tree,
# there it decides which bundle carries the run-wide assignment table.
assert RESERVED_COHORT_SCOPE == RESERVED_SCOPE, (
    "release_bundle.RESERVED_COHORT_SCOPE and the Snakefile's RESERVED_SCOPE name "
    "different scopes: {} and {}".format(RESERVED_COHORT_SCOPE, RESERVED_SCOPE))

# Resolve bundled code relative to the workflow, including when a parent
# repository launches it. Rules declare consumed source files as inputs.
SCRIPTS_DIR = Path(WORKFLOW_DIR) / "workflow" / "scripts"
LIB_DIR = Path(WORKFLOW_DIR) / "workflow" / "lib"
CONTAINER_EXEC = str(SCRIPTS_DIR / "container_exec.sh")


def named_values(mapping):
    """Render a locus-keyed setting as NAME=value pairs for one CLI option.

    Stages read the result with parse_named_values() in workflow/lib/runtime.R.
    An empty mapping renders as an empty string, which that parser accepts.
    """
    # Snakemake 7.32.4 inserts spaces into f-strings in .smk files under
    # Python 3.12, so this string is built without one.
    return ",".join("{}={}".format(name, value) for name, value in sorted(mapping.items()))


# obtain_analysis_image records the pinned image's runtime identity once, in a
# receipt beside the SIF under deps.resources. Analysis jobs declare that
# receipt. Apptainer and Singularity run the obtained SIF; Docker uses the image
# in its daemon. provenance.py owns image_sif_path(), image_receipt_path() and
# image validation. A receipt records more than the configured digest.
# The receipt lives with the obtained resources, not in the results directory,
# so a fresh results directory reuses it. obtain_imgt_reference declares the
# receipt and writes shared files under deps.resources, so a receipt written
# for each new results directory would rebuild them for every fresh run and
# stop concurrent runs on the working-directory lock.
ANALYSIS_IMAGE_REFERENCE = config["deps"]["immcantation_docker"]
ANALYSIS_SIF = image_sif_path(resources_dir, ANALYSIS_IMAGE_REFERENCE)
ANALYSIS_IMAGE_RECEIPT = str(image_receipt_path(resources_dir, ANALYSIS_IMAGE_REFERENCE))

# The pinned assembly image holding TRUST4, obtained the same way when a run
# assembles. Compute nodes can lack network access, and Singularity resolves a
# docker:// reference through its registry even with a populated cache (see
# test/README.md, Test on Gadi). Assembly jobs therefore run the SIF the main
# process pulls, not the configured docker:// URI.
ASSEMBLY_IMAGE_REFERENCE = config["deps"]["trust4_container"].removeprefix("docker://")
ASSEMBLY_SIF = image_sif_path(resources_dir, ASSEMBLY_IMAGE_REFERENCE)
ASSEMBLY_IMAGE_RECEIPT = str(image_receipt_path(resources_dir, ASSEMBLY_IMAGE_REFERENCE))

# The germline tree inside the analysis image, laid out as Immcantation's
# fetch_imgtdb.sh writes it: the configured vdj directory, its sibling constant
# directory and the download record IMGT.yaml two levels up. MakeDb and
# CreateGermlines read the vdj directory. Reference curation reads vdj and
# constant, and the audit reads all three. The receipt records what the audit
# found. Scope provenance declares it; annotation does not. An edit to
# audit_imgt_reference.py therefore reruns the audit without requiring IgBLAST.
GERMLINE_REFERENCE = str(LIB_DIR / "germline_reference.py")
GERMLINE_DIR = config["params"]["igblast"]["germline_dir"]
CONSTANT_DIR = str(Path(GERMLINE_DIR).parent / "constant")
IMGT_SNAPSHOT = str(Path(GERMLINE_DIR).parents[1] / "IMGT.yaml")
IMGT_REFERENCE_RECEIPT = str(results_dir / "shared/provenance/imgt_reference.json")
# IgBLAST's auxiliary data: NCBI's file inside the image, and the supplemented
# copy written beside the reference under the same name. igblastn reads the
# written copy, which supplements missing IMGT J alleles where
# supplement_auxiliary_data() finds a qualifying source entry.
SHIPPED_AUXILIARY_DATA = str(Path(config["params"]["igblast"]["igblast_db"]) / "optional_file" / "human_gl.aux")
AUXILIARY_DATA_FILE = resources_dir / "human_gl.aux"

# The user-facing release bundle: the allow-list of what leaves the pipeline, the
# documentation shipped beside it, and the data dictionary those docs are cut from.
RELEASE_BUNDLE = str(Path(WORKFLOW_DIR) / "workflow/lib/release_bundle.py")
RELEASE_DIR = Path(WORKFLOW_DIR) / "workflow" / "release"
RELEASE_README = str(RELEASE_DIR / "README_template.md")
DATA_DICTIONARY = str(Path(WORKFLOW_DIR) / "docs" / "data_dictionary.csv")
# runtime.R locates the shared R modules and parses their command-line options.
RUNTIME = str(LIB_DIR / "runtime.R")
# Declare these paths in the rules that consume each module, including imports
# made by another module. Keep independent analyses' code inputs separate.
TABLE_IO = str(LIB_DIR / "table_io.R")
AIRR_IO = str(LIB_DIR / "airr_io.R")
BASELINE_SELECTION = str(LIB_DIR / "baseline_selection.R")
IGH_SUMMARIES = str(LIB_DIR / "igh_summaries.R")
DIVERSITY = str(LIB_DIR / "diversity.R")
DONOR_DISTANCE_THRESHOLD = str(LIB_DIR / "donor_distance_threshold.R")
GENE_CALLS = str(LIB_DIR / "gene_calls.R")
GENE_USAGE = str(LIB_DIR / "gene_usage.R")
REANNOTATED_AIRR = str(LIB_DIR / "reannotated_airr.R")
ASSEMBLY_BACKBONE = str(LIB_DIR / "assembly_backbone.R")
IGH_FEATURES = str(LIB_DIR / "igh_features.R")
IGH_CLONES = str(LIB_DIR / "igh_clones.R")
IGH_GROUPING = str(LIB_DIR / "igh_grouping.R")
ISOTYPE_CALLS = str(LIB_DIR / "isotype_calls.R")
ISOTYPE_METRICS = str(LIB_DIR / "isotype_metrics.R")
ISOTYPE_RESOLUTION = str(LIB_DIR / "isotype_resolution.R")
OVERLAP_METRICS = str(LIB_DIR / "overlap_metrics.R")
POOLED_RECORDS = str(LIB_DIR / "pooled_records.R")
QC_AUDIT = str(LIB_DIR / "qc_audit.R")
REARRANGEMENT_KEYS = str(LIB_DIR / "rearrangement_keys.R")
REARRANGEMENT_GROUPS = str(LIB_DIR / "rearrangement_groups.R")
RECORD_QC = str(LIB_DIR / "record_qc.R")
REPERTOIRE_SUPPORT = str(LIB_DIR / "repertoire_support.R")
SCOPE_CLONE_THRESHOLD = str(LIB_DIR / "scope_clone_threshold.R")
SHM_REGIONS = str(LIB_DIR / "shm_regions.R")
SHM_SUMMARIES = str(LIB_DIR / "shm_summaries.R")
CONSENSUS_VOTES = str(LIB_DIR / "trust4_consensus_votes.R")
TRUST4_PARSING = str(LIB_DIR / "trust4_parsing.R")
SCOPE_TABLES_MODULE = str(LIB_DIR / "scope_tables.py")

# Keep each wildcard within one path component so directory names cannot be
# absorbed into identifiers.
wildcard_constraints:
    scope = r"[^/]+",
    donor_id = r"[^/]+",
    library_id = r"[^/]+",

# export_library_airr passes this list to library_processing_identity() in
# provenance.py. Editing any listed file changes each library's processing ID.
# Rule files and scope methods stay outside that ID; build_scope_provenance
# records the broader source state separately.
LIBRARY_PROCESSING_CODE = [str(LIB_DIR / (name + ".R")) for name in (
    "airr_io", "assembly_backbone", "gene_calls", "igh_grouping",
    "igh_clones", "isotype_calls", "isotype_resolution", "pooled_records",
    "qc_audit", "reannotated_airr",
    "rearrangement_groups", "rearrangement_keys", "record_qc", "runtime",
    "table_io", "trust4_parsing"
)] + [str(SCRIPTS_DIR / name) for name in (
    "container_exec.sh", "run_igblast_all_loci.sh", "parse_trust4_airr.R",
    "build_reannotated_airr.R", "filter_airr.R", "prepare_igh_grouping.R",
    "build_clonotypes.R", "qualify_isotype_calls.R", "export_library_airr.py"
)] + [str(LIB_DIR / name) for name in ("airr_exports.py", "alignment_export.py", "provenance.py")]
assert len(LIBRARY_PROCESSING_CODE) == len(set(LIBRARY_PROCESSING_CODE))
