"""Publish reference coordinates without changing Change-O reconstruction inputs.

IgBLAST hit evidence identifies the selected reference, including independently
reversed D hits that Change-O 1.3.4 cannot represent. Raw hits and native MakeDb
coordinates remain evidence; unsupported standard D fields are empty.
"""

import gzip
import re
from collections.abc import Mapping
from pathlib import Path

NATIVE_FIELDS = ("changeo_v_germline_start", "changeo_v_germline_end",
                 "changeo_d_germline_start", "changeo_d_germline_end", "changeo_d_cigar")
D_HIT_FIELDS = ("subject_id", "query_start", "query_end", "subject_start", "subject_end",
                "query_alignment", "germline_alignment", "btop")
ALIGNMENT_FIELDS = (*NATIVE_FIELDS, "igblast_v_subject_id", "igblast_j_subject_id",
                    *("igblast_d_" + field for field in D_HIT_FIELDS),
                    "d_alignment_status", "germline_reconstruction_status")
CIGAR = re.compile(r"(?:[1-9][0-9]*S)?(?:[1-9][0-9]*N)?(?:[1-9][0-9]*[=XMID])+(?:[1-9][0-9]*S)?(?:[1-9][0-9]*N)?")


def read_reference(path: Path) -> dict[str, str]:
    """Read the pinned ungapped reference by its exact allele identifiers."""
    result = {}
    name = None
    with Path(path).open() as handle:
        for line in handle:
            if line.startswith(">"):
                name = line[1:].split()[0]
                if name in result:
                    raise ValueError(f"Duplicate reference allele: {name}")
                result[name] = ""
            elif line.strip():
                if name is None:
                    raise ValueError("Reference sequence precedes its FASTA header")
                result[name] += line.strip().upper().replace(".", "").replace("-", "")
    return result


def add_alignment_evidence(annotations: dict[str, dict[str, str]], fields: list[str],
                           hits: Path) -> list[str]:
    """Join the first V/D/J hits from the pinned 17-column IgBLAST format."""
    for row in annotations.values():
        row.update(dict.fromkeys(ALIGNMENT_FIELDS, ""))
    seen = set()
    with gzip.open(hits, "rt") as handle:
        for line in handle:
            if not line.startswith(("V\t", "D\t", "J\t")):
                continue
            values = line.rstrip("\n").split("\t")
            if len(values) != 17:
                raise ValueError("IgBLAST hit must have 17 fields")
            segment, sequence_id = values[0].lower(), values[1].removeprefix("reversed|")
            key = (sequence_id, segment)
            if sequence_id not in annotations or key in seen:
                continue
            seen.add(key)
            row = annotations[sequence_id]
            if segment == "d":
                for field, index in zip(D_HIT_FIELDS, (2, 8, 9, 10, 11, 14, 15, 16)):
                    row["igblast_d_" + field] = values[index]
            else:
                row[f"igblast_{segment}_subject_id"] = values[2]
    return [*fields, *ALIGNMENT_FIELDS]


def cigar_interval(cigar: str) -> tuple[int, int, int, int, list[tuple[int, str]]]:
    """Return aligned query/reference intervals, excluding trailing offsets."""
    if CIGAR.fullmatch(cigar) is None:
        raise ValueError(f"Invalid AIRR CIGAR: {cigar}")
    operations = [(int(n), op) for n, op in re.findall(r"([0-9]+)([=XMIDNS])", cigar)]
    query, reference = 1, 1
    for n, op in operations:
        if op not in "SN":
            break
        if op == "S":
            query += n
        else:
            reference += n
    aligned = [(n, op) for n, op in operations if op not in "SN"]
    return (query, query + sum(n for n, op in aligned if op in "=XMI") - 1,
            reference, reference + sum(n for n, op in aligned if op in "=XMD") - 1, aligned)


def normalise_alignment(row: dict[str, str]) -> None:
    """Change only the export view; reconstruction and analytical rows stay native."""
    for field in NATIVE_FIELDS:
        source = field.removeprefix("changeo_")
        # add_makedb_fields supplies the unmodified native D CIGAR separately.
        if field != "changeo_d_cigar":
            row[field] = row.get(source, "")
    if row.get("v_cigar"):
        _, _, start, end, _ = cigar_interval(row["v_cigar"])
        row["v_germline_start"], row["v_germline_end"] = str(start), str(end)
    else:
        row["v_germline_start"] = row["v_germline_end"] = ""
    if row.get("igblast_d_subject_id"):
        start, end = (int(row["igblast_d_subject_" + edge]) for edge in ("start", "end"))
        status = "unavailable_reverse_orientation" if start > end else (
            "available" if row.get("d_cigar") else "unavailable_incomplete_evidence")
    elif row.get("d_call") or row.get("d_cigar"):
        status = "unavailable_incomplete_evidence"
    else:
        status = "not_applicable" if row.get("locus") in ("IGK", "IGL", "TRA", "TRG") else "no_hit"
    row["d_alignment_status"] = status
    if status != "available":
        for field in ("d_cigar", "d_germline_start", "d_germline_end"):
            row[field] = ""
    if not row.get("locus", "").startswith("IG"):
        reconstruction = "not_applicable"
    elif row.get("germline_alignment") and row.get("germline_alignment_d_mask"):
        reconstruction = "available"
    else:
        reconstruction = ("unavailable_reverse_d" if status == "unavailable_reverse_orientation"
                          else "unavailable_other")
    row["germline_reconstruction_status"] = reconstruction


def validate_alignment(row: Mapping[str, str], reference: Mapping[str, str] | None = None) -> None:
    """Validate standard intervals and, when supplied, their actual reference bases."""
    for segment in ("v", "d", "j"):
        cigar = row.get(segment + "_cigar", "")
        if not cigar:
            continue
        qs, qe, rs, re_, operations = cigar_interval(cigar)
        expected = (qs, qe, rs, re_)
        fields = tuple(f"{segment}_{axis}_{edge}" for axis in ("sequence", "germline")
                       for edge in ("start", "end"))
        actual = tuple(int(row[field]) for field in fields)
        if actual != expected or qe > len(row["sequence"]):
            raise ValueError(f"{row['sequence_id']}: {segment} alignment coordinates disagree with CIGAR")
        if reference is None:
            continue
        allele = row.get(f"igblast_{segment}_subject_id", "")
        if allele not in reference or allele not in row.get(segment + "_call", "").split(","):
            raise ValueError(f"{row['sequence_id']}: selected {segment} reference is absent from call/reference")
        germline = reference[allele]
        if re_ > len(germline):
            raise ValueError(f"{row['sequence_id']}: {segment} reference interval exceeds {allele}")
        q, r = qs - 1, rs - 1
        for n, op in operations:
            if op in "=XM":
                query, subject = row["sequence"][q:q+n].upper(), germline[r:r+n]
                if op == "=" and query != subject or op == "X" and any(a == b for a, b in zip(query, subject)):
                    raise ValueError(f"{row['sequence_id']}: {segment} CIGAR contradicts reference bases")
            q += n if op in "=XMI" else 0
            r += n if op in "=XMD" else 0
    status = row.get("d_alignment_status")
    if status and status != "available" and any(row.get(field) for field in
            ("d_cigar", "d_germline_start", "d_germline_end")):
        raise ValueError(f"{row['sequence_id']}: unavailable D alignment contains standard evidence")
    if reference is not None and row.get("igblast_d_subject_id"):
        allele = row["igblast_d_subject_id"]
        start, end = (int(row["igblast_d_subject_" + edge]) for edge in ("start", "end"))
        if allele not in reference or min(start, end) < 1 or max(start, end) > len(reference[allele]):
            raise ValueError(f"{row['sequence_id']}: native D reference interval is out of bounds")
        bases = reference[allele][min(start, end)-1:max(start, end)]
        if start > end:
            bases = bases.translate(str.maketrans("ACGTN", "TGCAN"))[::-1]
        if bases != row["igblast_d_germline_alignment"].replace("-", "").upper():
            raise ValueError(f"{row['sequence_id']}: native D traversal contradicts reference bases")
        qs, qe = (int(row["igblast_d_query_" + edge]) for edge in ("start", "end"))
        if not 1 <= qs <= qe <= len(row["sequence"]) or row["sequence"][qs-1:qe].upper() != row["igblast_d_query_alignment"].replace("-", "").upper():
            raise ValueError(f"{row['sequence_id']}: native D traversal contradicts query bases")

        if row.get("d_sequence_start") and row.get("d_sequence_end"):
            start, end = int(row["d_sequence_start"]), int(row["d_sequence_end"])
            if not qs <= start <= end <= qe:
                raise ValueError(f"{row['sequence_id']}: exported D query interval escapes its native hit")
