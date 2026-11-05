"""Library exports and shared AIRR schema and linkage primitives.

The official AIRR 2.0.0 schema is copied unchanged under workflow/schemas:
https://github.com/airr-community/airr-standards/blob/v2.0.0/specs/airr-schema.yaml
It is distributed by the AIRR Community under CC BY 4.0. Validation uses the
matching official Python package (airr 2.0.0, pinned in envs/release.yaml), not
the airr 1.6.1 package in the analysis image. AIRR Python library: Vander
Heiden et al. 2018, doi:10.3389/fimmu.2018.02206.

Repertoire IDs identify this pipeline's library observations, not donors or
projects. The data-processing ID also includes the caller's recorded processing
identity. The shared library export carries no scope-dependent clonal
assignment; annotate_member_assignments() adds one scope's clone_id to a copy.
"""

import csv
import gzip
import hashlib
import importlib.resources
import io
import itertools
import re
from collections.abc import Iterator, Mapping, Sequence
from pathlib import Path
from typing import Any
from urllib.parse import quote
from uuid import NAMESPACE_URL, uuid5


from .alignment_export import ALIGNMENT_FIELDS, normalise_alignment, validate_alignment

AIRR_VERSION = "2.0.0"
SCHEMA_SHA256 = "8e627adb91f14219c04a6c24194ac7cb30ef71165f8fdef174f5e329fbf66de5"
SCHEMA_PATH = Path(__file__).resolve().parents[1] / "schemas" / "airr-2.0.0.yaml"
ID_NAMESPACE = "https://github.com/swarbricklab/airrborne/"
# Columns the export adds to the source columns. source_sequence_id is an
# extension holding the original library-local ID. The other three are the AIRR
# Rearrangement fields that link each row to its Repertoire.json record.
LINK_COLUMNS = ("source_sequence_id", "repertoire_id", "sample_processing_id",
                "data_processing_id")
# Source columns the export requires, by their post-QC names. library_id and
# the four gene columns are extensions that the post-QC annotation already
# carries.
SOURCE_COLUMNS = ("library_id", "consensus_count", "v_gene", "j_gene",
                  "v_gene_status", "j_gene_status")

# Extension columns copied from the library's qualified isotype-call sidecar.
# They are filled on IGH rows and empty on every other locus. isotype_class and
# isotype_subclass are the resolve_isotype() names of the qualified call, as the
# isotype tables and clone members report them. The sidecar's qualified_c_call
# stays out of the export, because it equals c_call where c_call_callable is T
# and is empty elsewhere. The grouping extensions come from the grouping tables'
# own headers. The exported_AIRR rows of docs/data_dictionary.csv define every
# extension.
QUALIFIED_CALL_COLUMNS = ("c_call_present", "c_call_locus_match",
                          "c_call_enough_evidence", "c_call_callable",
                          "isotype_class", "isotype_subclass")
# Source columns the export renames to the release vocabulary. Internal tables
# keep their own names. The grouping tables name the member a row belongs to
# with the pipeline's internal inference_* columns, and the export names them
# after the member, as the members tables do. Internal grouping calls the
# TRUST4 contig a record comes from its backbone, and the export calls it a
# contig. The post-QC annotation writes each entity's key as *_gene and its
# genes as *_gene_label, and the export writes them as *_gene_key and
# *_gene_names, as the metrics tables do.
EXPORT_NAMES = {"inference_consensus_count": "member_consensus_count",
                "inference_representative_sequence_id": "member_representative_sequence_id",
                "is_inference_representative": "is_member_representative",
                "assembly_backbone_id": "contig_id",
                "backbone_representative_sequence_id": "contig_representative_sequence_id",
                "is_backbone_representative": "is_contig_representative",
                "v_gene": "v_gene_key", "j_gene": "j_gene_key",
                "v_gene_label": "v_gene_names", "j_gene_label": "j_gene_names"}
# Internal grouping names the status of a record whose contig's dominant variant
# fails after the contig's internal name, and the export names it after the
# contig.
EXPORT_VALUES = {"grouping_status": {"ineligible_backbone_representative":
                                     "ineligible_contig_representative"}}
# Pipeline tables write TRUE and FALSE. The export writes these extensions as
# T and F, the AIRR encoding its standard boolean fields already use.
BOOLEAN_ANNOTATIONS = ("grouping_eligible", "is_contig_representative",
                       "is_member_representative", "c_call_present",
                       "c_call_locus_match", "c_call_enough_evidence", "c_call_callable")
# Standard AIRR fields that MakeDb wrote and the post-QC annotation does not
# carry. Change-O CreateGermlines reads them to rebuild a clone's germline.
GERMLINE_FIELDS = ("np1_length", "np2_length", "d_germline_start", "d_germline_end",
                   "j_germline_start", "j_germline_end")
# Standard AIRR V-domain region fields that MakeDb wrote and the post-QC
# annotation does not carry. Change-O MakeDb writes them IMGT-gapped, so they
# keep its gap characters and can contain N. The post-QC annotation already
# carries cdr3, which equals MakeDb's.
REGION_FIELDS = ("fwr1", "cdr1", "fwr2", "cdr2", "fwr3", "fwr4")
# Characters of an IMGT-gapped region that are not sequenced bases. MakeDb
# writes '.' for an IMGT numbering gap and for each position before the first
# sequenced base, and '-' for a deletion. A region made only of these, such as
# the FWR1 of a record whose sequence starts in FWR3, has no sequenced base.
REGION_GAP_CHARACTERS = frozenset(".-")
# The column order of an export, which docs/data_dictionary.csv repeats in its
# exported_AIRR rows. A column the export lacks is skipped: the shared library
# export has no clone_id until annotate_member_assignments() writes it, and an
# export without annotations has no grouping, C-region or MakeDb columns. Any
# column outside this list keeps its source order after the listed columns.
EXPORT_ORDER = (
    # Identifiers and processing links.
    "sequence_id", "library_id", "source_sequence_id", "repertoire_id",
    "sample_processing_id", "data_processing_id",
    # Locus and productivity.
    "locus", "productive", "vj_in_frame", "stop_codon", "rev_comp", "changeo_rev_comp",
    # Gene calls, their resolved keys and names, and the qualified isotype.
    "v_call", "d_call", "j_call", "c_call", "v_gene_key", "v_gene_names",
    "j_gene_key", "j_gene_names", "isotype_class", "isotype_subclass",
    # Junction.
    "junction", "junction_aa", "junction_length", "cdr3", "cdr3_aa",
    # Support and membership.
    "consensus_count", "member_consensus_count", "clone_id", "contig_id",
    "contig_representative_sequence_id", "is_contig_representative",
    "member_representative_sequence_id", "is_member_representative",
    "grouping_eligible", "grouping_status",
    # Sequence and alignments.
    "sequence", "sequence_alignment", "germline_alignment", "germline_alignment_d_mask",
    # V-domain regions.
    "fwr1", "cdr1", "fwr2", "cdr2", "fwr3", "fwr4",
    # CIGARs, coordinates, NP lengths and identities.
    "v_cigar", "d_cigar", "j_cigar", "v_sequence_start", "v_sequence_end",
    "v_germline_start", "v_germline_end", "d_sequence_start", "d_sequence_end",
    "d_germline_start", "d_germline_end", "j_sequence_start", "j_sequence_end",
    "j_germline_start", "j_germline_end", "np1_length", "np2_length",
    "v_identity", "j_identity",
    # Call qualification.
    "v_gene_status", "j_gene_status", "v_call_candidates_complete",
    "j_call_candidates_complete", "c_call_present", "c_call_locus_match",
    "c_call_enough_evidence", "c_call_candidates_complete", "c_call_callable",
    "c_informative_positions",
    *ALIGNMENT_FIELDS,
)
# Post-QC columns the export leaves out: reference functionality and orphon
# context, which the pinned reference metadata holds per allele, and the D-call
# candidate completeness that no released analysis reads.
EXPORT_DROPPED = ("v_reference_functionality", "d_reference_functionality",
                  "j_reference_functionality", "c_reference_functionality",
                  "v_reference_orphon", "d_reference_orphon", "j_reference_orphon",
                  "c_reference_orphon", "d_call_candidates_complete")
AIRR_BOOLEANS = {"TRUE": "T", "FALSE": "F", "": ""}
# Joins a library ID to a library-local ID, as pooled_records.R does.
QUALIFIED_ID_SEPARATOR = "~"


def airr_boolean(value: str, field: str) -> str:
    """Write one pipeline boolean in the AIRR T/F encoding; blank stays blank."""
    if value not in AIRR_BOOLEANS:
        raise ValueError(f"{field} is not TRUE, FALSE or blank: {value!r}")
    return AIRR_BOOLEANS[value]


def validate_cigar_query_spans(row: Mapping[str, str]) -> None:
    """Check AIRR operations against oriented-query coordinates for each segment.

    Schema validation checks field types, but cannot detect swapped I/D
    operators. The query span consumes M, =, X and I. S/N encode unaligned
    offsets at the ends. Germline coordinates may include IMGT spacers and do
    not enter this check. Empty CIGARs retain AIRR's missing-evidence meaning.
    """
    for segment in ("v", "d", "j"):
        cigar = row.get(f"{segment}_cigar", "")
        if not cigar:
            continue
        context = f"AIRR {row['sequence_id']} {segment}_cigar"
        if not re.fullmatch(r"(?:[1-9][0-9]*S)?(?:[1-9][0-9]*N)?"
                            r"(?:[1-9][0-9]*[=XMID])+"
                            r"(?:[1-9][0-9]*S)?(?:[1-9][0-9]*N)?", cigar):
            raise ValueError(f"{context} has invalid operations: {cigar}")
        tokens = [(int(n), op) for n, op in re.findall(r"([0-9]+)([=XMIDNS])", cigar)]
        query_start = tokens[0][0] + 1 if tokens[0][1] == "S" else 1
        query_span = sum(n for n, op in tokens if op in "=XMI")
        try:
            start = int(row[f"{segment}_sequence_start"])
            end = int(row[f"{segment}_sequence_end"])
        except (KeyError, ValueError) as error:
            raise ValueError(f"{context} lacks query coordinates") from error
        if start != query_start or end != start + query_span - 1 or end > len(row["sequence"]):
            raise ValueError(f"{context} contradicts query coordinates {start}..{end}: {cigar}")


def read_library_annotations(grouping_files: Sequence[Path], calls_file: Path,
                             library_id: str) -> tuple[dict[str, dict[str, str]], list[str]]:
    """Join producer decisions for every source row, preserving textual values.

    Boolean decisions are rewritten as T and F; every other value is kept as written.
    """
    annotations, loci = {}, {}
    columns = None
    for path in grouping_files:
        with path.open(newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            header = reader.fieldnames
            if not header or not {"library_id", "sequence_id", "locus"} <= set(header):
                raise ValueError("Grouping annotations require source identities")
            if columns is not None and header != columns:
                raise ValueError("Grouping branch headers differ")
            columns = header
            source_fields = [field for field in columns if field not in ("library_id", "sequence_id", "locus")]
            fields = [EXPORT_NAMES.get(field, field) for field in source_fields]
            for row in reader:
                key = row["sequence_id"]
                if row["library_id"] != library_id or key in annotations:
                    raise ValueError("Grouping annotations contain a duplicate or another library")
                extra = {field: row[field] for field in source_fields}
                # These identifiers refer to this library's assemblies and rows,
                # and the exported sequence_id is library-qualified. Qualify them
                # the same way so a recipient can join them inside the file. An
                # empty value stays empty.
                for field in ("assembly_backbone_id", "backbone_representative_sequence_id",
                              "inference_representative_sequence_id"):
                    if extra[field]:
                        extra[field] = library_id + QUALIFIED_ID_SEPARATOR + extra[field]
                extra = {EXPORT_NAMES.get(field, field): value for field, value in extra.items()}
                extra.update(dict.fromkeys(QUALIFIED_CALL_COLUMNS, ""))
                annotations[key], loci[key] = extra, row["locus"]
    if columns is None:
        raise ValueError("No grouping annotations supplied")
    seen = set()
    with calls_file.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames or not {"library_id", "sequence_id", *QUALIFIED_CALL_COLUMNS} <= set(reader.fieldnames):
            raise ValueError("Qualified calls lack their declared columns")
        for row in reader:
            key = row["sequence_id"]
            if row["library_id"] != library_id or key in seen or loci.get(key) != "IGH":
                raise ValueError("Qualified calls contain a duplicate or unmatched IGH key")
            seen.add(key)
            annotations[key].update({field: row[field] for field in QUALIFIED_CALL_COLUMNS})
    for extra in annotations.values():
        for field in BOOLEAN_ANNOTATIONS:
            if field in extra:
                extra[field] = airr_boolean(extra[field], field)
        for field, values in EXPORT_VALUES.items():
            if field in extra:
                extra[field] = values.get(extra[field], extra[field])
    if seen != {key for key, locus in loci.items() if locus == "IGH"}:
        raise ValueError("Qualified calls do not cover every common-QC IGH record")
    return annotations, fields + list(QUALIFIED_CALL_COLUMNS)


def add_makedb_fields(annotations: dict[str, dict[str, str]], fields: list[str],
                      makedb_file: Path) -> list[str]:
    """Add MakeDb's np and D/J germline coordinates and V-domain regions to every record.

    `makedb_file` is the library's gzip-compressed MakeDb output, keyed by the same
    library-local sequence_id as the post-QC annotation. Every annotated record
    must appear there. Germline coordinates are copied as written. A region is
    copied as written when it holds a sequenced base, any character outside
    REGION_GAP_CHARACTERS, and is written empty otherwise, so an empty region
    means none of its positions was sequenced. A library without post-QC records
    has a MakeDb output with neither column group, and that file is not read.
    Returns the field order with GERMLINE_FIELDS and REGION_FIELDS appended.
    """
    added = (*GERMLINE_FIELDS, *REGION_FIELDS)
    if not annotations:
        return fields + list(added)
    found = set()
    with gzip.open(makedb_file, "rt", newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames or not {"sequence_id", *added} <= set(reader.fieldnames):
            raise ValueError(f"MakeDb output lacks germline coordinates or V-domain regions: {makedb_file}")
        for row in reader:
            key = row["sequence_id"]
            if key in annotations:
                if key in found:
                    raise ValueError(f"MakeDb output repeats {key!r}")
                found.add(key)
                annotations[key].update({field: row[field] for field in GERMLINE_FIELDS})
                annotations[key].update({field: row[field] if set(row[field]) - REGION_GAP_CHARACTERS
                                         else "" for field in REGION_FIELDS})
    if found != set(annotations):
        raise ValueError("MakeDb output does not cover every post-QC record")
    return fields + list(added)


def _airr() -> Any:
    """Refuse a validator whose bundled schema differs from the declared one."""
    import airr

    installed = importlib.resources.files("airr").joinpath("specs/airr-schema.yaml")
    if airr.__version__ != AIRR_VERSION:
        raise RuntimeError(f"AIRR exports require airr=={AIRR_VERSION}; found {airr.__version__}")
    for name, content in (("installed", installed.read_bytes()),
                          ("vendored", SCHEMA_PATH.read_bytes())):
        if hashlib.sha256(content).hexdigest() != SCHEMA_SHA256:
            raise RuntimeError(f"The {name} AIRR schema differs from the pinned {AIRR_VERSION} schema")
    return airr


def _identifier(value: str, name: str) -> str:
    if not isinstance(value, str) or not value or value != value.strip():
        raise ValueError(f"{name} must be a nonempty string without surrounding whitespace")
    if any(character in value for character in "\t\r\n"):
        raise ValueError(f"{name} contains a table delimiter")
    return value


def export_ids(library_id: str, processing_identity: str) -> dict[str, str]:
    """Return scope-independent identifiers for one library.

    The repertoire and sample-processing identifiers are derived from library_id
    alone, so they are the same in every scope and survive a processing change.
    The data-processing identifier includes the processing identity, so changed
    code, settings or resources give it a new value.
    """
    library_id = _identifier(library_id, "library_id")
    processing_identity = _identifier(processing_identity, "processing_identity")
    base = ID_NAMESPACE + "library/" + quote(library_id, safe="")
    names = {"repertoire_id": base,
             "sample_processing_id": base + "/sample-processing",
             "data_processing_id": base + "/processing/" + quote(processing_identity, safe="")}
    return {key: uuid5(NAMESPACE_URL, value).urn for key, value in names.items()}


def _header(reader: csv.DictReader, name: Path) -> list[str]:
    fields = reader.fieldnames
    if not fields or any(not field for field in fields) or len(fields) != len(set(fields)):
        raise ValueError(f"Missing or duplicate columns in {name}")
    _airr().schema.RearrangementSchema.validate_header(fields)
    return fields


def _export_rows(reader: csv.DictReader, library_id: str, processing_identity: str,
                 annotations: Mapping[str, Mapping[str, str]] | None,
                 annotation_fields: Sequence[str]):
    ids = export_ids(library_id, processing_identity)
    seen = set()
    schema = _airr().schema.RearrangementSchema
    for row in reader:
        if None in row or any(value is None for value in row.values()):
            raise ValueError(f"Malformed AIRR row {reader.line_num}")
        native_id = _identifier(row["sequence_id"], "source sequence_id")
        if native_id in seen:
            raise ValueError(f"Duplicate source sequence_id: {native_id!r}")
        seen.add(native_id)
        if row.get("library_id") != library_id:
            raise ValueError(f"Source row does not belong to library {library_id!r}")
        if row.get("clone_id", ""):
            raise ValueError("Shared library annotation must not contain clone assignments")
        exported = {**{EXPORT_NAMES.get(field, field): value for field, value in row.items()
                       if field not in EXPORT_DROPPED},
                    "sequence_id": library_id + QUALIFIED_ID_SEPARATOR + native_id,
                    "source_sequence_id": native_id, **ids}
        if annotations is not None:
            if native_id not in annotations:
                raise ValueError(f"Missing authority annotation for {native_id!r}")
            extra = annotations[native_id]
            if set(extra) != set(annotation_fields):
                raise ValueError(f"Authority columns differ for {native_id!r}")
            if any(not isinstance(value, str) for value in extra.values()):
                raise ValueError("Authority values must preserve TSV string representations")
            exported.update(extra)
        if row["rev_comp"] not in ("T", "F", ""):
            raise ValueError("Source rev_comp must use AIRR boolean encoding")
        exported["changeo_rev_comp"] = row["rev_comp"]
        exported["rev_comp"] = "F"
        validate_cigar_query_spans(exported)
        if "d_alignment_status" in exported:
            normalise_alignment(exported)
            validate_alignment(exported)
        schema.validate_row(exported)
        validate_cigar_query_spans(exported)
        yield exported
    if annotations is not None and seen != set(annotations):
        raise ValueError("Authority annotations contain records absent from post-QC annotation")


def order_export_columns(columns: Sequence[str]) -> list[str]:
    """Write the supplied columns in EXPORT_ORDER, then any others in their source order."""
    return ([column for column in EXPORT_ORDER if column in columns]
            + [column for column in columns if column not in EXPORT_ORDER])


def _export_fields(fields: Sequence[str], annotation_fields: Sequence[str]) -> list[str]:
    missing = set(SOURCE_COLUMNS) - set(fields)
    if missing:
        raise ValueError(f"Post-QC annotation lacks canonical fields: {sorted(missing)}")
    if "clone_id" in annotation_fields:
        raise ValueError("Shared library annotation must not contain clone assignments")
    kept = [EXPORT_NAMES.get(field, field) for field in fields if field not in EXPORT_DROPPED]
    extra = list(LINK_COLUMNS) + ["changeo_rev_comp"] + list(annotation_fields)
    if len(extra) != len(set(extra)) or set(extra).intersection(kept):
        raise ValueError("Export columns would overwrite source or linkage columns")
    return order_export_columns(kept + extra)


def export_library_airr(source: Path, destination: Path, library_id: str,
                        processing_identity: str, *,
                        annotations: Mapping[str, Mapping[str, str]] | None = None,
                        annotation_fields: Sequence[str] = ()) -> dict[str, str]:
    """Compress source cells with qualified IDs and optional authoritative fields.

    EXPORT_DROPPED columns are left out and EXPORT_NAMES columns renamed; every
    other source cell is kept except rev_comp. MakeDb stores the query in V
    orientation, so the standard flag is F and changeo_rev_comp preserves the
    historical transformation. order_export_columns() gives the column order.

    Authority fields must cover every post-QC record exactly. They may describe
    eligibility and transferred support, but cannot alter the source annotation.
    Gzip timestamps and filenames are omitted, so the same input compressed by the
    same library produces the same bytes.
    """
    if bool(annotation_fields) != (annotations is not None):
        raise ValueError("Supply authority annotations and their column order together")
    source, destination = Path(source), Path(destination)
    ids = export_ids(library_id, processing_identity)
    destination.parent.mkdir(parents=True, exist_ok=True)
    with source.open(newline="", encoding="utf-8") as source_handle:
        reader = csv.DictReader(source_handle, delimiter="\t")
        fields = _export_fields(_header(reader, source), annotation_fields)
        with destination.open("wb") as raw:
            with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as compressed:
                with io.TextIOWrapper(compressed, encoding="utf-8", newline="") as handle:
                    writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
                    writer.writeheader()
                    writer.writerows(_export_rows(reader, library_id, processing_identity,
                                                  annotations, annotation_fields))
    return ids


def iter_library_airr(path: Path, library_id: str, processing_identity: str, *,
                      member_assignments: bool = False,
                      reference: Mapping[str, str] | None = None) -> Iterator[dict[str, str]]:
    """Yield each export row after checking the official schema, query CIGAR spans and identities.

    The shared library export carries no clone_id. A scope copy written by
    annotate_member_assignments() passes `member_assignments`: every member
    representative then carries a clone_id and no other row does.
    """
    schema = _airr().schema.RearrangementSchema
    ids = export_ids(library_id, processing_identity)
    seen = set()
    with gzip.open(path, "rt", newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = _header(reader, Path(path))
        if not set(LINK_COLUMNS).union(EXPORT_NAMES.get(column, column)
                                       for column in SOURCE_COLUMNS).issubset(fields):
            raise ValueError(f"Missing export linkage columns in {path}")
        if member_assignments and not {"clone_id", "is_member_representative"} <= set(fields):
            raise ValueError(f"Scope AIRR copy lacks its member assignments: {path}")
        for row in reader:
            if None in row or any(value is None for value in row.values()):
                raise ValueError(f"Malformed AIRR row {reader.line_num} in {path}")
            schema.validate_row(row)
            if row.get("rev_comp") != "F" or row.get("changeo_rev_comp") not in ("T", "F", ""):
                raise ValueError("Export orientation must describe the stored V-oriented query")
            validate_cigar_query_spans(row)
            if "d_alignment_status" in row:
                validate_alignment(row, reference)
            native_id = _identifier(row["source_sequence_id"], "source_sequence_id")
            if row["sequence_id"] != library_id + QUALIFIED_ID_SEPARATOR + native_id or row["sequence_id"] in seen:
                raise ValueError(f"Invalid or duplicate qualified sequence_id in {path}")
            seen.add(row["sequence_id"])
            if row["library_id"] != library_id or any(row[key] != value for key, value in ids.items()):
                raise ValueError(f"AIRR row linkage differs from Repertoire metadata in {path}")
            if not member_assignments:
                if row.get("clone_id", ""):
                    raise ValueError("Shared library annotation must not contain clone assignments")
            elif bool(row["clone_id"]) != (row["is_member_representative"] == "T"):
                raise ValueError(f"clone_id must name exactly the member representatives in {path}")
            yield row


def validate_library_airr(path: Path, library_id: str, processing_identity: str, *,
                          member_assignments: bool = False) -> int:
    """Check every row of one export with iter_library_airr() and return the row count."""
    return sum(1 for _ in iter_library_airr(path, library_id, processing_identity,
                                            member_assignments=member_assignments))


def annotate_member_assignments(source: Path, destination: Path,
                                assignments: Mapping[str, str]) -> int:
    """Write one library export with each member's scope clone or clonotype in clone_id.

    `assignments` maps a qualified sequence_id to the clone_id of the IGH clone, or
    the clonotype_id of the clonotype, that the member it names belongs to in one
    scope. Every member representative must be assigned and no other row may be:
    contig siblings and records outside grouping keep an empty clone_id,
    because a sibling's junction can differ from its member's. clone_id takes
    its EXPORT_ORDER position after member_consensus_count; every other cell is
    copied unchanged. Gzip timestamps and filenames are omitted, as in
    export_library_airr(). Returns the number of rows written.
    """
    source, destination = Path(source), Path(destination)
    destination.parent.mkdir(parents=True, exist_ok=True)
    assigned, count = set(), 0
    with gzip.open(source, "rt", newline="", encoding="utf-8") as source_handle:
        reader = csv.DictReader(source_handle, delimiter="\t")
        fields = list(reader.fieldnames or [])
        if "clone_id" in fields or "is_member_representative" not in fields:
            raise ValueError(f"{source} is not a shared library export")
        with destination.open("wb") as raw:
            with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as compressed:
                with io.TextIOWrapper(compressed, encoding="utf-8", newline="") as handle:
                    writer = csv.DictWriter(handle, fieldnames=order_export_columns([*fields, "clone_id"]),
                                            delimiter="\t", lineterminator="\n")
                    writer.writeheader()
                    for row in reader:
                        key = row["sequence_id"]
                        member = row["is_member_representative"] == "T"
                        if member != (key in assignments):
                            raise ValueError(f"Member assignment differs from the member flag for {key!r}")
                        if member:
                            assigned.add(key)
                        writer.writerow({**row, "clone_id": assignments.get(key, "")})
                        count += 1
    if assigned != set(assignments):
        raise ValueError("Member assignments name records absent from the library export")
    return count


def validate_export_content(source: Path, exported: Path, library_id: str,
                            processing_identity: str, *,
                            annotations: Mapping[str, Mapping[str, str]] | None = None,
                            annotation_fields: Sequence[str] = ()) -> int:
    """Compare every decompressed cell after the declared export additions and drops."""
    count = 0
    with Path(source).open(newline="", encoding="utf-8") as source_handle:
        with gzip.open(exported, "rt", newline="", encoding="utf-8") as export_handle:
            before = csv.DictReader(source_handle, delimiter="\t")
            after = csv.DictReader(export_handle, delimiter="\t")
            expected_fields = _export_fields(_header(before, Path(source)), annotation_fields)
            if after.fieldnames != expected_fields:
                raise ValueError("Decompressed export columns differ from declared source additions")
            expected = _export_rows(before, library_id, processing_identity, annotations, annotation_fields)
            for count, (left, right) in enumerate(itertools.zip_longest(expected, after), 1):
                if left != right:
                    raise ValueError(f"Decompressed export differs from its source at data row {count}")
    return count
