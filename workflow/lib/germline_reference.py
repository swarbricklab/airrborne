"""Build and check the reference genes used to recognise receptor sequences.

A germline reference describes receptor gene segments before rearrangement or
somatic mutation. An allele is a named version of a gene. This module reads
IMGT/GENE-DB alleles from the pinned analysis image so assembly and annotation
can be checked against the same source sequences.

build_reference_metadata() decides which complete V/J reference sequences are
identical. gene_calls.R uses those classes to interpret assignments naming
several genes. build_reference() chooses one record per allele for TRUST4.
build_auxiliary_data() supplies J-allele frame and junction-boundary entries
needed by IgBLAST. audit_reference() checks the written files and IgBLAST
databases against their sources and returns the receipt read by scope provenance.

References:
  IMGT/GENE-DB: Lefranc et al. 2015, doi:10.1093/nar/gku1056.
  IgBLAST 1.22.0 and its auxiliary data: Ye et al. 2013, doi:10.1093/nar/gkt382.
  TRUST4 BuildImgtAnnot.pl:
    https://github.com/liulab-dfci/TRUST4/blob/v1.1.10/BuildImgtAnnot.pl
  Immcantation reference and database scripts in the pinned image under
  /usr/local/sbin: fetch_imgtdb.sh, imgt2igblast.sh, clean_imgtdb.py.
"""

import csv
import io
import fnmatch
import hashlib
import subprocess
from collections.abc import Callable, Mapping, Sequence
from pathlib import Path
from typing import Any

RECEIPT_SCHEMA_VERSION = "1.1"
SPECIES = "human"
# IgBLAST family names and the locus prefixes their databases combine, as
# imgt2igblast.sh globs them: imgt_human_IG?V.fasta and so on.
FAMILY_LOCI = {"ig": "IG", "tr": "TR"}
SEGMENTS = ("v", "d", "j", "c")
# Paths below params.igblast.igblast_db: the databases run_igblast_all_loci.sh
# passes to igblastn, and NCBI's auxiliary data file the supplement starts from.
AUXILIARY_DATA = "optional_file/human_gl.aux"
DATABASE_DIR = "database"
# The supplemented auxiliary data file, written beside the reference under
# IgBLAST's own name for it.
AUXILIARY_DATA_NAME = "human_gl.aux"
_CASE = str.maketrans("acgtn", "ACGTN")


def read_imgt_entries(path: Path) -> list[dict[str, str]]:
    """Read an IMGT/GENE-DB FASTA as allele, functionality and sequence, in file order.

    The header is pipe-delimited, as in
    >M99641|IGHV1-18*01|Homo_sapiens|F|V-REGION|188..483|... The allele is field
    2 and the functionality (F, ORF or P, sometimes in brackets) is field 4.
    Sequence lines are joined unchanged, so the result keeps IMGT gap dots and
    lower case.
    """
    entries: list[dict[str, str]] = []
    header: list[str] | None = None
    parts: list[str] = []
    with Path(path).open() as handle:
        for line in handle:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if header is not None:
                    entries.append(_imgt_entry(header, parts))
                header = line[1:].split("|")
                if len(header) < 4 or not header[1]:
                    raise ValueError(f"{path} must use IMGT pipe-delimited headers with the allele in "
                                     f"field 2 and the functionality in field 4. Found: {line}")
                parts = []
            elif header is not None:
                parts.append(line)
    if header is not None:
        entries.append(_imgt_entry(header, parts))
    return entries


def _imgt_entry(header: list[str], parts: list[str]) -> dict[str, str]:
    """Keep the region and completeness field from the original IMGT header."""
    completeness = header[13].strip() if len(header) >= 15 else "unknown"
    return {"allele": header[1], "functionality": header[3].strip(),
            "region": header[4].strip() if len(header) > 4 else "unknown",
            "completeness": completeness or "complete", "sequence": "".join(parts)}


REFERENCE_METADATA_COLUMNS = ("allele", "gene", "region", "functionality", "completeness",
                              "identity_id", "identity_label", "orphon")


def build_reference_metadata(files: Sequence[Path]) -> list[dict[str, str]]:
    """Describe reference alleles and group identical complete V/J sequences.

    Keep the first entry per allele in input-file order. Only a complete
    V-REGION or J-REGION made entirely of A/C/G/T after gap removal receives
    an identity_id. Equal full sequences within a region share that ID, even
    when their gene names differ. Other entries retain their metadata but
    have no identity class. Functionality is recorded without filtering.

    reannotated_airr.R reads these classes before gene_calls.R interprets
    multi-gene V/J calls. Equality of one observed alignment fragment does
    not establish equality of the reference alleles.
    """
    kept: dict[str, dict[str, str]] = {}
    for file in files:
        for entry in read_imgt_entries(file):
            kept.setdefault(entry["allele"], entry)
    rows = []
    for allele, entry in sorted(kept.items()):
        sequence = clean_sequence(entry["sequence"])
        region = entry["region"]
        eligible = (region in {"V-REGION", "J-REGION"} and entry["completeness"] == "complete"
                    and bool(sequence) and set(sequence) <= set("ACGT"))
        identity = region[0] + "SEQ_" + hashlib.sha256(sequence.encode()).hexdigest() if eligible else ""
        gene = allele.split("*", 1)[0]
        rows.append({"allele": allele, "gene": gene, "region": region,
                     "functionality": entry["functionality"], "completeness": entry["completeness"],
                     "identity_id": identity, "identity_label": "", "orphon": "T" if "/OR" in gene else "F"})
    labels: dict[str, set[str]] = {}
    for row in rows:
        if row["identity_id"]:
            labels.setdefault(row["identity_id"], set()).add(row["gene"])
    for row in rows:
        if row["identity_id"]:
            row["identity_label"] = ",".join(sorted(labels[row["identity_id"]]))
    return rows


def render_reference_metadata(rows: Sequence[Mapping[str, str]]) -> str:
    """Write the allele evidence used for grouping and reference-functionality context."""
    output = io.StringIO()
    writer = csv.DictWriter(output, fieldnames=REFERENCE_METADATA_COLUMNS, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
    return output.getvalue()


def reference_identity_classes(rows: Sequence[Mapping[str, str]]) -> list[dict[str, Any]]:
    """Expose every eligible sequence class and its exact allele membership in the receipt."""
    classes: dict[str, dict[str, Any]] = {}
    for row in rows:
        if row["identity_id"]:
            entry = classes.setdefault(row["identity_id"], {"identity_id": row["identity_id"],
                "label": row["identity_label"], "region": row["region"], "alleles": []})
            entry["alleles"].append(row["allele"])
    return [classes[key] for key in sorted(classes)]


def read_imgt_records(path: Path) -> list[tuple[str, str]]:
    """Read an IMGT/GENE-DB FASTA as (allele, sequence) pairs in file order."""
    return [(entry["allele"], entry["sequence"]) for entry in read_imgt_entries(path)]


def clean_sequence(sequence: str) -> str:
    """Spell a sequence the way IgBLAST databases and the R identity classes do."""
    return sequence.replace(".", "").upper()


def keep_first(records: Sequence[tuple[str, str, str]]) -> tuple[list[tuple[str, str]], list[dict[str, Any]]]:
    """Keep the first record for each allele name across (allele, sequence, file) triples.

    Kept sequences have acgtn uppercased and their gap dots retained. Each
    skipped record is described by its allele, its file, the file of the record
    kept instead and whether the two spell the same nucleotides once the gap
    dots are removed.
    """
    kept: dict[str, tuple[str, str]] = {}
    order: list[str] = []
    skipped: list[dict[str, Any]] = []
    for allele, sequence, file in records:
        if allele in kept:
            kept_sequence, kept_file = kept[allele]
            skipped.append({"allele": allele, "file": file, "kept_file": kept_file,
                            "identical_sequence": clean_sequence(sequence) == clean_sequence(kept_sequence)})
            continue
        kept[allele] = (sequence.translate(_CASE), file)
        order.append(allele)
    return [(allele, kept[allele][0]) for allele in order], skipped


def render_fasta(records: Sequence[tuple[str, str]]) -> str:
    """Write one two-line record per allele with a bare allele header."""
    return "".join(f">{allele}\n{sequence}\n" for allele, sequence in records)


def germline_files(germline_dir: Path, constant_dir: Path) -> list[Path]:
    """List the vdj FASTAs and then the constant FASTAs, each in name order.

    fetch_imgtdb.sh writes human/vdj with one FASTA per V, D and J locus and
    human/constant with one per constant locus. MakeDb and CreateGermlines
    read human/vdj directly. IMGT.yaml at the tree root records the source
    and download date; read_imgt_snapshot() reads that record.
    """
    files = sorted(Path(germline_dir).glob("*.fasta")) + sorted(Path(constant_dir).glob("*.fasta"))
    if not files:
        raise ValueError(f"No IMGT germline FASTA found under {germline_dir} or {constant_dir}. "
                         "params.igblast.germline_dir must name the analysis image's human/vdj "
                         "directory, beside human/constant.")
    return files


def build_reference(germline_dir: Path, constant_dir: Path) -> tuple[str, dict[str, Any]]:
    """Curate the reference text and describe the files and records it came from.

    Keep the first record per allele across the V/D/J files and then the
    constant files. Return the FASTA text and a description of input files,
    record counts, SHA-256 values and every skipped duplicate record.

    The FASTA follows TRUST4 BuildImgtAnnot.pl's format: bare allele headers,
    acgtn uppercased and IMGT gap dots retained. TRUST4 removes those dots in
    SeqSet::InputRefFa. BuildImgtAnnot.pl concatenates entries sharing a name
    and accession; this module instead records each skipped duplicate so the
    keep-first choice can be checked. For constant alleles, that choice keeps
    the first spliced form IMGT lists. Released C calls come from IgBLAST.
    """
    files = germline_files(germline_dir, constant_dir)
    triples: list[tuple[str, str, str]] = []
    described = []
    for file in files:
        records = read_imgt_records(file)
        triples.extend((allele, sequence, file.name) for allele, sequence in records)
        described.append({"path": str(file), "records": len(records),
                          "sha256": hashlib.sha256(file.read_bytes()).hexdigest()})
    kept, skipped = keep_first(triples)
    if not kept:
        raise ValueError(f"No IMGT record found in {len(files)} germline FASTA files under "
                         f"{germline_dir} and {constant_dir}.")
    curation = {"germline_files": described, "records_read": len(triples), "alleles_written": len(kept),
                "duplicate_records_skipped": skipped, "allele_metadata": build_reference_metadata(files)}
    return render_fasta(kept), curation


def write_reference(text: str, output: Path) -> None:
    """Replace the reference atomically so a failed write leaves no partial file."""
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_name(output.name + ".tmp")
    temporary.write_text(text)
    temporary.replace(output)


def read_imgt_snapshot(path: Path) -> dict[str, Any]:
    """Read the source, download date and species that fetch_imgtdb.sh recorded in IMGT.yaml."""
    import yaml

    record = yaml.safe_load(Path(path).read_text())
    if not isinstance(record, dict) or not {"source", "date", "species"} <= set(record):
        raise ValueError(f"{path} must record source, date and species, as fetch_imgtdb.sh writes them.")
    return {"path": str(path), "source": str(record["source"]), "date": str(record["date"]),
            "species": [str(value) for value in record["species"]]}


def reference_alleles(text: str) -> dict[str, str]:
    """Map each allele of a curated reference to its cleaned sequence.

    Headers are bare allele names; ignore text after the first whitespace.
    Repeated names fail because the audit requires one sequence per allele.
    The assignment classes are read separately from reference metadata by
    workflow/lib/reannotated_airr.R.
    """
    alleles: dict[str, str] = {}
    allele = None
    parts: list[str] = []
    for line in text.splitlines():
        if line.startswith(">"):
            if allele is not None:
                alleles[allele] = clean_sequence("".join(parts))
            allele = line[1:].strip().split()[0] if line[1:].strip() else ""
            if not allele or allele in alleles:
                raise ValueError(f"Curated reference names an allele twice or leaves a header blank: {line}")
            parts = []
        elif allele is not None:
            parts.append(line.strip())
    if allele is not None:
        alleles[allele] = clean_sequence("".join(parts))
    if not alleles:
        raise ValueError("Curated reference holds no FASTA record.")
    return alleles


def database_name(family: str, segment: str) -> str:
    """Name a database as run_igblast_all_loci.sh passes it to igblastn."""
    return f"imgt_{SPECIES}_{family}_{segment}"


def expected_database_alleles(files: Sequence[Path], family: str, segment: str) -> dict[str, str]:
    """Rebuild the allele set imgt2igblast.sh gave one database, from the same files.

    V, D and J databases combine the vdj files matching imgt_human_<LOCUS>?<SEGMENT>.fasta
    for the family's loci, and the C database the constant files matching
    imgt_human_<LOCUS>?C.fasta, each concatenated in name order. clean_imgtdb.py
    then keeps the first record per allele name, removes the gap dots and
    uppercases.
    """
    pattern = f"imgt_{SPECIES}_{FAMILY_LOCI[family]}?{segment.upper()}.fasta"
    triples = [(allele, sequence, file.name)
               for file in sorted(files, key=lambda path: path.name) if fnmatch.fnmatch(file.name, pattern)
               for allele, sequence in read_imgt_records(file)]
    kept, _ = keep_first(triples)
    return {allele: clean_sequence(sequence) for allele, sequence in kept}


def dump_blast_database(database: Path) -> dict[str, str]:
    """Read every allele and sequence from a BLAST nucleotide database with blastdbcmd."""
    result = subprocess.run(["blastdbcmd", "-db", str(database), "-entry", "all", "-outfmt", "%a\t%s"],
                            capture_output=True, text=True, check=False)
    if result.returncode != 0:
        raise ValueError(f"blastdbcmd cannot read {database}: {result.stderr.strip()}")
    alleles: dict[str, str] = {}
    for line in result.stdout.splitlines():
        if not line.strip():
            continue
        allele, _, sequence = line.partition("\t")
        if not sequence or allele in alleles:
            raise ValueError(f"{database} returned a blank or repeated allele: {line}")
        alleles[allele] = sequence.upper()
    if not alleles:
        raise ValueError(f"{database} holds no sequence.")
    return alleles


def compare_alleles(expected: Mapping[str, str], observed: Mapping[str, str]) -> dict[str, list[str]]:
    """List alleles missing from observed, unexpected in it, and those whose sequences differ."""
    shared = expected.keys() & observed.keys()
    return {"missing_alleles": sorted(expected.keys() - observed.keys()),
            "unexpected_alleles": sorted(observed.keys() - expected.keys()),
            "sequence_mismatches": sorted(allele for allele in shared if expected[allele] != observed[allele])}


def _auxiliary_lines(text: str) -> list[str]:
    return [line.strip() for line in text.splitlines() if line.strip() and not line.strip().startswith("#")]


def parse_auxiliary_data(text: str) -> dict[str, tuple[str, ...]]:
    """Map each J allele in an auxiliary data file to its entry fields.

    An entry is the allele name followed by the coding frame start (0-based),
    the chain type, the CDR3 end position (0-based) and the bases beyond the
    last complete codon, as NCBI's IgBLAST setup guide describes the file. A
    pseudogene entry can stop after the chain type. A name repeated with the
    same fields is kept once, and a name repeated with different fields is an
    error, because IgBLAST would then use whichever line it read last.
    """
    entries: dict[str, tuple[str, ...]] = {}
    for line in _auxiliary_lines(text):
        fields = line.split()
        name, values = fields[0], tuple(fields[1:])
        if len(values) not in (2, 4) or not values[0].isdigit() or (len(values) == 4 and not (
                values[2].isdigit() and values[3].isdigit())):
            raise ValueError(f"Auxiliary data entry must hold a frame, a chain type and optionally a CDR3 end "
                             f"and extra bases. Found: {line}")
        if name in entries and entries[name] != values:
            raise ValueError(f"Auxiliary data names {name} twice with different fields")
        entries.setdefault(name, values)
    return entries


def read_j_alleles(germline_dir: Path) -> dict[str, dict[str, str]]:
    """Map each J allele of the vdj directory to its cleaned sequence and functionality.

    The J files are the ones named imgt_human_<LOCUS>J.fasta, the files
    imgt2igblast.sh combines into the two J databases. The first record for
    each allele name is kept, as the databases keep it.
    """
    alleles: dict[str, dict[str, str]] = {}
    for file in sorted(Path(germline_dir).glob(f"imgt_{SPECIES}_*J.fasta")):
        for entry in read_imgt_entries(file):
            alleles.setdefault(entry["allele"], {"sequence": clean_sequence(entry["sequence"]),
                                                 "functionality": entry["functionality"], "file": file.name})
    if not alleles:
        raise ValueError(f"No J allele found in {germline_dir}: expected files named imgt_{SPECIES}_<LOCUS>J.fasta.")
    return alleles


def _gene(allele: str) -> str:
    return allele.split("*", 1)[0]


IDENTICAL_SEQUENCE = "identical_sequence"
SAME_GENE_EQUAL_LENGTH = "same_gene_equal_length"


def supplement_auxiliary_data(shipped_text: str, j_alleles: Mapping[str, Mapping[str, str]]) -> tuple[str, dict[str, Any]]:
    """Add entries for J alleles the shipped auxiliary data file lacks, by copying another allele's.

    An allele without an entry receives, in this order of preference:

    1. the entry of an allele whose cleaned sequence is identical, since the
       frame and CDR3 end are positions in that same sequence;
    2. the entry of a same-gene allele of equal length whose entry has all
       four fields, when the allele keeps that sibling's conserved codon (the
       three bases after the CDR3 end) and the extra-base count the frame
       implies.

    Every qualifying source must give the same entry. An allele with no
    qualifying source is left without an entry and reported with its IMGT
    functionality; NCBI's own file has no entry for the human J pseudogenes
    that have no functional allele.

    Before adding anything, both rules are checked on the shipped file itself:
    every pair of identical-sequence alleles, and every pair of same-gene,
    equal-length alleles, that both have full entries must carry identical
    entries. A file that breaks either rule stops the build.

    Returns the shipped text with the new lines appended, and a report naming
    the pairs checked, each added entry with its rule and source alleles, and
    each allele left without an entry.
    """
    entries = parse_auxiliary_data(shipped_text)
    full = {name: values for name, values in entries.items() if len(values) == 4 and name in j_alleles}
    by_gene: dict[str, list[str]] = {}
    by_sequence: dict[str, list[str]] = {}
    for name in sorted(j_alleles):
        by_gene.setdefault(_gene(name), []).append(name)
        by_sequence.setdefault(j_alleles[name]["sequence"], []).append(name)

    def check_pairs(groups: Mapping[str, list[str]], same_length_only: bool, label: str) -> int:
        checked = 0
        for names in groups.values():
            with_entry = [name for name in names if name in full]
            for index, first in enumerate(with_entry):
                for second in with_entry[index + 1:]:
                    if same_length_only and len(j_alleles[first]["sequence"]) != len(j_alleles[second]["sequence"]):
                        continue
                    checked += 1
                    if full[first] != full[second]:
                        raise ValueError(f"The shipped auxiliary data gives {first} and {second}, {label}, "
                                         f"different entries ({full[first]} and {full[second]}), so an entry "
                                         "cannot be copied between such alleles on this image.")
        return checked

    identical_checked = check_pairs(by_sequence, False, "alleles with identical sequences")
    sibling_checked = check_pairs(by_gene, True, "same-gene alleles of equal length")
    added: list[dict[str, Any]] = []
    without: list[dict[str, str]] = []
    for name in sorted(j_alleles):
        if name in entries:
            continue
        sequence = j_alleles[name]["sequence"]
        qualifying: dict[tuple[str, ...], list[str]] = {}
        rule = IDENTICAL_SEQUENCE
        for twin in by_sequence[sequence]:
            if twin != name and twin in full:
                qualifying.setdefault(full[twin], []).append(twin)
        if not qualifying:
            rule = SAME_GENE_EQUAL_LENGTH
            for sibling in by_gene[_gene(name)]:
                if sibling == name or sibling not in full:
                    continue
                other = j_alleles[sibling]["sequence"]
                frame, _, cdr3_stop, extra = full[sibling]
                stop = int(cdr3_stop)
                if (len(other) == len(sequence) and sequence[stop + 1:stop + 4] == other[stop + 1:stop + 4]
                        and (len(sequence) - int(frame)) % 3 == int(extra)):
                    qualifying.setdefault(full[sibling], []).append(sibling)
        if len(qualifying) > 1:
            raise ValueError(f"Source alleles give {name} different auxiliary entries: {qualifying}")
        if not qualifying:
            without.append({"allele": name, "functionality": j_alleles[name]["functionality"],
                            "reason": "no allele with an identical sequence, and no same-gene allele of equal "
                                      "length with a full entry and the same conserved codon"})
            continue
        (values, sources), = qualifying.items()
        added.append({"allele": name, "functionality": j_alleles[name]["functionality"], "rule": rule,
                      "copied_from": sources, "frame": int(values[0]), "chain_type": values[1],
                      "cdr3_stop": int(values[2]), "extra_bps": int(values[3])})
    # The added lines carry no comment. igblastn 1.22.0 in the pinned image
    # was killed for memory (over 3 GB on a four-sequence query, measured on
    # 28 September 2026 through Docker) when a 246-character comment line
    # followed NCBI's entries, while the same entries without it ran in 46 MB.
    # The receipt records what was added and why.
    text = shipped_text if shipped_text.endswith("\n") else shipped_text + "\n"
    text += "".join(f"{entry['allele']}\t{entry['frame']}\t{entry['chain_type']}\t{entry['cdr3_stop']}\t"
                    f"{entry['extra_bps']}\n" for entry in added)
    report = {"shipped_entries": len(_auxiliary_lines(shipped_text)), "shipped_names": len(entries),
              "identical_pairs_checked": identical_checked, "sibling_pairs_checked": sibling_checked,
              "supplemented": added, "j_alleles_without_entry": without}
    return text, report


def build_auxiliary_data(germline_dir: Path, shipped: Path) -> tuple[str, dict[str, Any]]:
    """Build J-allele frame and CDR3-boundary entries from NCBI's shipped file.

    IgBLAST reports no CDR3 when its first-listed top J hit lacks an entry.
    The image's NCBI file predates its IMGT download. supplement_auxiliary_data()
    copies entries from identical sequences or qualifying same-gene alleles;
    alleles without a qualifying source remain listed as missing in the report.
    Return supplemented text and its source and coverage report, without writing.
    """
    text, report = supplement_auxiliary_data(Path(shipped).read_text(), read_j_alleles(germline_dir))
    report["shipped_path"] = str(shipped)
    report["shipped_sha256"] = hashlib.sha256(Path(shipped).read_bytes()).hexdigest()
    return text, report


def audit_reference(*, reference: Path, auxiliary: Path, metadata: Path, germline_dir: Path, constant_dir: Path, snapshot: Path,
                    igblast_db: Path, image: Mapping[str, Any],
                    read_database: Callable[[Path], Mapping[str, str]] = dump_blast_database) -> dict[str, Any]:
    """Check the curated reference and auxiliary data against their sources and the IgBLAST databases.

    Rebuild the reference, metadata and auxiliary data and compare them with
    the files on disk. Dump each IgBLAST database through blastdbcmd to check
    its allele names and sequences against the reference and its source files.
    Return the audit receipt used by scope provenance.

    image is the analysis image identity from the analysis image receipt.
    auxiliary is the auxiliary data file written beside the reference.
    read_database reads one BLAST database; the default runs blastdbcmd.

    Raises ValueError, naming every disagreement, when the reference or the
    auxiliary data on disk differs from a fresh build from the image's files,
    or when any IgBLAST database holds a different allele set or sequence from
    the reference or from the germline files it was built from. J alleles left
    without an auxiliary entry are recorded, not failed on.
    """
    reference, auxiliary = Path(reference), Path(auxiliary)
    text, curation = build_reference(germline_dir, constant_dir)
    metadata = Path(metadata)
    metadata_bytes = metadata.read_bytes()
    metadata_matches = metadata_bytes == render_reference_metadata(curation["allele_metadata"]).encode()
    on_disk = reference.read_bytes()
    rebuild_matches = on_disk == text.encode()
    alleles = reference_alleles(on_disk.decode())
    files = [Path(entry["path"]) for entry in curation["germline_files"]]
    problems = []
    if not metadata_matches:
        problems.append(f"{metadata} differs from the original IMGT allele metadata.")
    if not rebuild_matches:
        problems.append(f"{reference} differs from a fresh curation of the germline files under "
                        f"{germline_dir} and {constant_dir}.")
    shipped_auxiliary = Path(igblast_db) / AUXILIARY_DATA
    auxiliary_text, auxiliary_report = build_auxiliary_data(germline_dir, shipped_auxiliary)
    auxiliary_on_disk = auxiliary.read_bytes()
    auxiliary_matches = auxiliary_on_disk == auxiliary_text.encode()
    if not auxiliary_matches:
        problems.append(f"{auxiliary} differs from a fresh supplement of {shipped_auxiliary} for the J alleles "
                        f"under {germline_dir}.")
    databases = {}
    for family in FAMILY_LOCI:
        for segment in SEGMENTS:
            path = Path(igblast_db) / DATABASE_DIR / database_name(family, segment)
            observed = read_database(path)
            # The reference holds every family and segment, so the check that
            # matters is that each database allele is in it with the same
            # sequence; alleles of other databases are expected to be absent.
            against_reference = {
                "missing_from_reference": sorted(observed.keys() - alleles.keys()),
                "sequence_mismatches": sorted(allele for allele in observed.keys() & alleles.keys()
                                              if observed[allele] != alleles[allele])}
            against_sources = compare_alleles(expected_database_alleles(files, family, segment), observed)
            entry = {"path": str(path), "alleles": len(observed), "against_reference": against_reference,
                     "against_germline_files": against_sources}
            databases[f"{family}_{segment}"] = entry
            for label, comparison in (("reference", against_reference), ("germline files", against_sources)):
                found = {key: value for key, value in comparison.items() if value}
                if found:
                    problems.append(f"{path} disagrees with the {label}: {found}")
    # The written file is what igblastn reads, so its coverage of the J
    # databases is what the receipt reports.
    written_entries = parse_auxiliary_data(auxiliary_on_disk.decode())
    j_alleles = {allele for family in FAMILY_LOCI
                 for allele in read_database(Path(igblast_db) / DATABASE_DIR / database_name(family, "j"))}
    without_entry = {entry["allele"] for entry in auxiliary_report["j_alleles_without_entry"]}
    if j_alleles - set(written_entries) != without_entry:
        problems.append(f"{auxiliary} leaves J database alleles without an entry that the supplement report "
                        f"does not name: {sorted((j_alleles - set(written_entries)) ^ without_entry)}")
    receipt = {
        "schema_version": RECEIPT_SCHEMA_VERSION,
        "image": dict(image),
        "reference": {"name": reference.name, "sha256": hashlib.sha256(on_disk).hexdigest(),
                      "size_bytes": len(on_disk)},
        "reference_metadata": {"name": metadata.name, "sha256": hashlib.sha256(metadata_bytes).hexdigest(),
                               "rebuild_matches": metadata_matches},
        "identity_classes": reference_identity_classes(curation["allele_metadata"]),
        "imgt_snapshot": read_imgt_snapshot(snapshot),
        "germline_files": curation["germline_files"],
        "curation": {"records_read": curation["records_read"], "alleles_written": curation["alleles_written"],
                     "duplicate_records_skipped": curation["duplicate_records_skipped"],
                     "rebuild_matches_reference": rebuild_matches},
        "igblast_databases": databases,
        "auxiliary_data": {
            "shipped": {"path": auxiliary_report["shipped_path"], "sha256": auxiliary_report["shipped_sha256"],
                        "entries": auxiliary_report["shipped_entries"], "names": auxiliary_report["shipped_names"]},
            "written": {"name": auxiliary.name, "sha256": hashlib.sha256(auxiliary_on_disk).hexdigest(),
                        "entries": len(_auxiliary_lines(auxiliary_on_disk.decode()))},
            "rebuild_matches_written": auxiliary_matches,
            "identical_pairs_checked": auxiliary_report["identical_pairs_checked"],
            "sibling_pairs_checked": auxiliary_report["sibling_pairs_checked"],
            "supplemented": auxiliary_report["supplemented"],
            "j_alleles_without_entry": auxiliary_report["j_alleles_without_entry"],
            "entries_without_j_allele": sorted(set(written_entries) - j_alleles)},
        "agreement": {"reference_rebuild": rebuild_matches, "metadata_rebuild": metadata_matches, "auxiliary_rebuild": auxiliary_matches,
                      "igblast_databases": not any(
                          value for entry in databases.values() for comparison in
                          (entry["against_reference"], entry["against_germline_files"])
                          for value in comparison.values())},
    }
    if problems:
        raise ValueError("The germline reference audit found disagreements:\n" + "\n".join(problems)
                         + "\nRebuild the reference and auxiliary data from the pinned image, or review the "
                           "image's germline tree and IgBLAST databases.")
    return receipt
