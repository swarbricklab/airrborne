"""Reconstruct annotation results independently from retained scientific inputs."""

from pathlib import Path
from test.support.reference_tables import check_fields, number, read_table, index_rows
import json


def verify_reference_receipt(root: Path) -> int:
    """Check the run's reference receipt against the fixture's recorded reference.

    The fixture was constructed from the reference whose SHA-256 is recorded in
    test/data/construction/construction.json. The receipt the audit_imgt_reference rule
    wrote must name that file and record agreement between the reference, a
    fresh curation and every IgBLAST database. Returns the number of databases
    compared.
    """
    receipt = json.loads((root / "shared/provenance/imgt_reference.json").read_text())
    construction = json.loads((Path(__file__).resolve().parents[3] / "data/construction/construction.json").read_text())
    assert receipt["reference"]["sha256"] == construction["reference_sha256"], \
        "The audited reference differs from the fixture's construction reference"
    assert receipt["agreement"] == {"reference_rebuild": True, "metadata_rebuild": True, "auxiliary_rebuild": True,
                                    "igblast_databases": True}, receipt["agreement"]
    classes = receipt["identity_classes"]
    allele_classes = {allele: row["identity_id"] for row in classes for allele in row["alleles"]}
    assert allele_classes["IGLJ2*01"] == allele_classes["IGLJ3*01"]
    assert allele_classes["IGHV3-30*02"] == allele_classes["IGHV3-30-5*02"]
    assert not {"IGHV4-59*09", "IGHV4-61*07"} & allele_classes.keys(), \
        "Partial reference entries must not establish sequence identity"
    assert receipt["reference_metadata"]["rebuild_matches"]
    assert len(allele_classes) == sum(len(row["alleles"]) for row in classes)
    auxiliary = receipt["auxiliary_data"]
    assert auxiliary["written"]["sha256"] != auxiliary["shipped"]["sha256"] or not auxiliary["supplemented"], \
        "A supplemented auxiliary file must differ from the shipped file"
    assert auxiliary["sibling_pairs_checked"] > 0, "The copy rule was not checked on the shipped file"
    for entry in auxiliary["supplemented"]:
        assert entry["copied_from"] and entry["allele"] not in entry["copied_from"], entry
    # In this image every J allele left without an entry is an IMGT pseudogene
    # of a gene with no other allele. A functional allele here needs review.
    for entry in auxiliary["j_alleles_without_entry"]:
        assert entry["functionality"] == "P", f"Functional J allele without an auxiliary entry: {entry}"
    databases = receipt["igblast_databases"]
    assert set(databases) == {f"{family}_{segment}" for family in ("ig", "tr") for segment in "vdjc"}
    for name, entry in databases.items():
        assert entry["alleles"] > 0, f"{name}: empty IgBLAST database"
        for comparison in ("against_reference", "against_germline_files"):
            assert not any(entry[comparison].values()), f"{name}: {entry[comparison]}"
    assert receipt["curation"]["alleles_written"] == sum(
        entry["alleles"] for name, entry in databases.items()), "Reference alleles differ from database alleles"
    assert receipt["imgt_snapshot"]["source"] and receipt["imgt_snapshot"]["date"]
    return len(databases)
