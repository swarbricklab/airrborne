"""Build scope AIRR metadata and validate staged links with the official schema."""

from collections import Counter
from collections.abc import Iterable, Mapping, Sequence
from pathlib import Path, PurePosixPath
import re
from typing import Any

from .airr_exports import (
    AIRR_VERSION, SCHEMA_SHA256, _airr, _identifier, export_ids, iter_library_airr,
    validate_library_airr,
)

SEMANTIC_LOCI = ("IGH", "IGK", "IGL", "TRA", "TRB", "TRD", "TRG")
# The counts airr/validation.json holds for each locus.
SEMANTIC_COUNTS = ("n_records", "n_member_representatives", "n_members_with_several_contig_sources",
                   "n_members_with_support_above_source_record", "n_blank_germline_alignment",
                   "n_blank_germline_alignment_d_mask", "n_noncanonical_junction_aa_records",
                   "n_noncanonical_junction_aa_members")
CANONICAL_JUNCTION_AA = re.compile(r"^C.*[FW]$")
# Statements airr/validation.json carries beside those counts.
SEMANTIC_STATEMENTS = {
    "query_orientation": (
        "sequence is stored in V orientation and standard rev_comp is F. "
        "changeo_rev_comp preserves the source transformation; query coordinates index sequence."),
    "v_germline_coordinates": (
        "Standard V coordinates index the ungapped selected reference. "
        "changeo_v_germline_start/end preserve native IMGT reconstruction coordinates."),
    "cigar_query_spans_scope": (
        "Query spans are checked separately from reference spans and bases. "
        "Unsupported reverse-D standard fields are empty with native evidence retained."),
    "germline_alignment_d_mask": (
        "Germlines are reconstructed for the IG loci only, so germline_alignment_d_mask is "
        "empty on every TCR record, and on an IG record whose germline could not be reconstructed."),
    "noncanonical_junction_aa": (
        f"A noncanonical junction peptide fails the regular expression {CANONICAL_JUNCTION_AA.pattern}. "
        "This screen is a diagnostic, and no record is filtered on it."),
}


def _required_record(name: str) -> dict[str, Any]:
    """Keep required AIRR keys; callers add only identifiers they can populate."""
    schema = _airr().schema.Schema(name)
    fields = set(schema.required)
    return {key: None for key in schema.properties if key in fields}


def _relative_file(value: str) -> str:
    path = PurePosixPath(value)
    if path.is_absolute() or ".." in path.parts or not value.endswith(".tsv.gz"):
        raise ValueError(f"Rearrangement file must be a release-relative .tsv.gz path: {value!r}")
    return str(path)


def build_scope_repertoires(scope_id: str, libraries: Sequence[Mapping[str, Any]],
                            exported_files: Mapping[str, str], processing_identity: str, *,
                            processing_fields: Mapping[str, str | None] | None = None,
                            library_processing_fields: Mapping[str, Mapping[str, str | None]] | None = None) -> dict[str, Any]:
    """Describe the supplied scope libraries, asserting only what this pipeline fixes.

    Species, template class, single-cell status, physical linkage, repertoire
    type and subject synthetic status are constants of this pipeline: human
    germline reference, bulk RNA, no cell barcodes, no chain pairing, observed
    repertoires of real donors. Every export is the primary annotation, and each
    Repertoire is named by its library_id. Preparation details are unavailable
    from the manifest. Required unavailable values stay null, with a shared explanation
    in Info.description. Optional unavailable fields are omitted. Subjects are
    donors and samples are tissues.

    Processing fields describe the recorded computation, not the packaging run.
    Optional project annotations belong to the web report only.
    """
    _identifier(scope_id, "scope_id")
    processing_fields = dict(processing_fields or {})
    library_processing_fields = dict(library_processing_fields or {})
    allowed_processing = set(_airr().schema.Schema("DataProcessing").required)
    if set(processing_fields) - allowed_processing:
        raise ValueError("Processing fields must name the declared AIRR DataProcessing descriptions")
    if any(set(fields) - allowed_processing for fields in library_processing_fields.values()):
        raise ValueError("Library processing fields must name the declared AIRR DataProcessing descriptions")
    if set(library_processing_fields) - {row["library_id"] for row in libraries}:
        raise ValueError("Library processing fields contain an out-of-scope library")
    records = []
    seen = set()
    tissues = {}
    for library in libraries:
        library_id, tissue_id, donor_id = (_identifier(library[key], key)
                                          for key in ("library_id", "tissue_id", "donor_id"))
        if library_id in seen:
            raise ValueError(f"Duplicate scope library: {library_id!r}")
        seen.add(library_id)
        if tissue_id in tissues and tissues[tissue_id] != donor_id:
            raise ValueError(f"Tissue {tissue_id!r} maps to several donors")
        tissues[tissue_id] = donor_id
        if library_id not in exported_files:
            raise ValueError(f"Missing Rearrangement file for library {library_id!r}")
        ids = export_ids(library_id, processing_identity)
        study = _required_record("Study")
        subject = _required_record("Subject")
        subject.update(subject_id=donor_id, synthetic=False,
                       species={"id": "NCBITAXON:9606", "label": "Homo sapiens"})
        # The manifest does not supply receptor-targeted PCR or protocol extent.
        # Those fields remain unavailable and are reported as metadata gaps.
        sample = _required_record("SampleProcessing")
        sample.update(sample_id=tissue_id, sample_processing_id=ids["sample_processing_id"],
                      single_cell=False, template_class="RNA", physical_linkage="none")
        processing = _required_record("DataProcessing")
        processing.update(processing_fields)
        processing.update(library_processing_fields.get(library_id, {}))
        processing.update(data_processing_id=ids["data_processing_id"], primary_annotation=True,
                          data_processing_files=[_relative_file(exported_files[library_id])])
        record = {"repertoire_id": ids["repertoire_id"], "repertoire_name": library_id,
                  "repertoire_type": "observed", "study": study, "subject": subject,
                  "sample": [sample], "data_processing": [processing]}
        records.append(record)
    if not records:
        raise ValueError("A scope must contain at least one library")
    # File maps may describe shared products for more libraries than this scope.
    return {"Info": {"title": "Bulk RNA repertoire library exports", "version": AIRR_VERSION,
                     "description": (
                         "Study, clinical, library-preparation and sequencing-run metadata are not supplied to this export. "
                         "Required unavailable values remain null; they do not denote negative findings. "
                         "Paired-end FASTQs and library, tissue and donor identifiers do not establish these values. "
                         "Optional project annotations and timepoints are used only in the web report. "
                         "The primary processing identity describes shared library annotation. "
                         "Clone assignments are scope-specific additions identified by scope_id "
                         "and the scope provenance; they are not a second primary annotation.")},
            "scope_id": scope_id, "processing_identity": processing_identity, "Repertoire": records}


def _metadata_gaps(schema: Any, obj: Mapping[str, Any], path: str = "") -> list[dict[str, str]]:
    """List nullability and required-key gaps using the official schema itself."""
    gaps = []
    for field, spec in schema.properties.items():
        if spec.get("x-airr", {}).get("deprecated"):
            continue
        location = path + field
        value = obj.get(field)
        if value is None:
            if field in schema.required and field not in obj:
                gaps.append({"field": location, "reason": "required field missing"})
            elif spec.get("x-airr", {}).get("nullable") is False:
                reason = "nonnullable field unavailable" if field in obj else "nonnullable field omitted"
                gaps.append({"field": location, "reason": reason})
            continue
        if "$ref" in spec:
            nested = _airr().schema.Schema(spec["$ref"].split("/")[-1])
            gaps.extend(_metadata_gaps(nested, value, location + "."))
        elif spec.get("type") == "array" and "$ref" in spec.get("items", {}):
            nested = _airr().schema.Schema(spec["items"]["$ref"].split("/")[-1])
            for index, item in enumerate(value):
                gaps.extend(_metadata_gaps(nested, item, f"{location}[{index}]."))
    return gaps


def validate_repertoire_metadata(data: Mapping[str, Any],
                               libraries: Sequence[Mapping[str, str]],
                               exported_files: Mapping[str, str],
                               processing_identity: str) -> dict[str, Any]:
    """Validate AIRR structure and expected identities without opening exports.

    Required clinical and preparation metadata are unavailable in these inputs.
    AIRR 2.0's check_nullable=False validates the supplied structure without
    inventing values. It is not a full schema-validity or MiAIRR completeness
    claim. The separate strict result and gaps make that limit machine-readable.
    """
    airr = _airr()
    if data.get("processing_identity") != processing_identity:
        raise ValueError("Repertoire processing identity differs from the library exports")
    schema = airr.schema.Schema("Repertoire")
    records = data.get("Repertoire")
    if not isinstance(records, list) or not records:
        raise ValueError("Repertoire must be a nonempty array")
    expected = {row["library_id"]: row for row in libraries}
    if len(expected) != len(libraries):
        raise ValueError("Duplicate libraries in the scope crosswalk")
    seen = set()
    files = set()
    gaps = []
    strict_errors = []
    for record in records:
        schema.validate_object(record, missing=False, nonairr=False, check_nullable=False)
        record_gaps = _metadata_gaps(schema, record)
        missing_keys = [gap["field"] for gap in record_gaps if gap["reason"] == "required field missing"]
        if missing_keys:
            raise ValueError(f"Missing AIRR required metadata keys: {missing_keys}")
        library_id = record.get("repertoire_name")
        if library_id not in expected or library_id in seen:
            raise ValueError(f"Unexpected or duplicate Repertoire library: {library_id!r}")
        seen.add(library_id)
        ids = export_ids(library_id, processing_identity)
        samples, processes = record.get("sample"), record.get("data_processing")
        if not samples or len(samples) != 1 or not processes or len(processes) != 1:
            raise ValueError("Each library Repertoire must have one sample and one processing record")
        sample, processing = samples[0], processes[0]
        if (record.get("repertoire_id") != ids["repertoire_id"]
                or sample.get("sample_processing_id") != ids["sample_processing_id"]
                or processing.get("data_processing_id") != ids["data_processing_id"]
                or record.get("subject", {}).get("subject_id") != expected[library_id]["donor_id"]
                or sample.get("sample_id") != expected[library_id]["tissue_id"]):
            raise ValueError(f"Repertoire identity or donor/specimen linkage differs for {library_id!r}")
        links = processing.get("data_processing_files")
        if not isinstance(links, list) or len(links) != 1:
            raise ValueError("Each processing record must link one shared library Rearrangement file")
        relative = _relative_file(links[0])
        if library_id not in exported_files or relative != _relative_file(exported_files[library_id]):
            raise ValueError(f"Rearrangement file link differs for library {library_id!r}")
        if relative in files:
            raise ValueError("Rearrangement file is linked by several libraries")
        files.add(relative)
        gaps.extend({"library_id": library_id, **gap} for gap in record_gaps)
        try:
            schema.validate_object(record, missing=False, nonairr=False, check_nullable=True)
        except airr.ValidationError as error:
            strict_errors.append({"library_id": library_id, "error": str(error)})
    if seen != set(expected):
        raise ValueError("Repertoire records do not cover every scope library")
    # structural_valid is a literal: schema.validate_object() above raises on a
    # structural failure, so a bundle that exists always passed. The strict result
    # and the gap list below are the values that vary.
    return {"schema_version": AIRR_VERSION, "schema_sha256": SCHEMA_SHA256,
            "validator": "airr-python", "validator_version": airr.__version__,
            "structural_valid": True, "check_nullable": False,
            "strict_valid": not strict_errors, "strict_errors": _group_metadata_reasons(strict_errors, ("error",)),
            "strict_metadata_gaps": _group_metadata_reasons(gaps, ("field", "reason")), "linked_libraries": len(seen),
            "miairr_complete": False}


def _group_metadata_reasons(rows: Sequence[Mapping[str, str]], keys: tuple[str, ...]) -> list[dict]:
    """Record a missing-field reason once with all affected library IDs."""
    groups = {}
    for row in rows:
        key = tuple(row[field] for field in keys)
        groups.setdefault(key, []).append(row["library_id"])
    return [{**dict(zip(keys, key)), "library_ids": sorted(libraries)}
            for key, libraries in sorted(groups.items())]


def _count_semantics(rows: Iterable[Mapping[str, str]], by_locus: dict[str, dict[str, int]]) -> int:
    """Add one library's scope-copy rows to the per-locus counts and return the row count.

    A contig representative names its member in member_representative_sequence_id.
    Those links are held for one library at a time, because a member belongs to
    one library.
    """
    contig_sources = Counter()
    n_rows = 0
    for row in rows:
        n_rows += 1
        counts = by_locus[row["locus"]]
        member = row["is_member_representative"] == "T"
        noncanonical = CANONICAL_JUNCTION_AA.search(row["junction_aa"]) is None
        counts["n_records"] += 1
        counts["n_member_representatives"] += member
        counts["n_members_with_support_above_source_record"] += (
            member and float(row["member_consensus_count"]) > float(row["consensus_count"]))
        counts["n_blank_germline_alignment"] += not row["germline_alignment"]
        counts["n_blank_germline_alignment_d_mask"] += not row["germline_alignment_d_mask"]
        counts["n_noncanonical_junction_aa_records"] += noncanonical
        counts["n_noncanonical_junction_aa_members"] += member and noncanonical
        if row["is_contig_representative"] == "T":
            contig_sources[row["locus"], row["member_representative_sequence_id"]] += 1
    for (locus, _), n_sources in contig_sources.items():
        by_locus[locus]["n_members_with_several_contig_sources"] += n_sources > 1
    return n_rows


def validate_scope_repertoires(data: Mapping[str, Any], release_root: Path,
                               libraries: Sequence[Mapping[str, str]],
                               exported_files: Mapping[str, str],
                               processing_identity: str, *,
                               member_assignments: bool = False,
                               reference: Mapping[str, str] | None = None) -> dict[str, Any]:
    """Validate the actual staged AIRR files and count their linked records.

    `member_assignments` marks scope copies, which carry clone_id. Their receipt
    also holds semantic_checks: SEMANTIC_STATEMENTS and each locus's
    SEMANTIC_COUNTS, read while each file is validated. Every locus appears, with
    zero counts where the scope has no record.
    """
    result = validate_repertoire_metadata(data, libraries, exported_files, processing_identity)
    root = Path(release_root).resolve()
    counts = {}
    by_locus = {locus: dict.fromkeys(SEMANTIC_COUNTS, 0) for locus in SEMANTIC_LOCI}
    for library in libraries:
        name = library["library_id"]
        path = (root / _relative_file(exported_files[name])).resolve()
        if not path.is_relative_to(root):
            raise ValueError("Rearrangement file escapes the bundle")
        if member_assignments:
            counts[name] = _count_semantics(
                iter_library_airr(path, name, processing_identity, member_assignments=True, reference=reference), by_locus)
        else:
            counts[name] = validate_library_airr(path, name, processing_identity)
    result["linked_records"] = counts
    result["cigar_query_spans_valid"] = True
    result["cigar_reference_spans_and_bases_valid"] = reference is not None
    if member_assignments:
        result["semantic_checks"] = {**SEMANTIC_STATEMENTS, "by_locus": by_locus}
    return result
