"""Check cross-branch counts and audit context without repeating scientific gates.

Run: python test/run_tests.py --component record_qc --level unit
"""

import csv
import gzip
from pathlib import Path
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[4]))
from workflow.lib.scope_qc_audit import (
    COLUMNS, LOCI, QC_COLUMNS, TRAIL_COLUMNS, TRAIL_STAGES, canonical_totals, checkpoint,
    collect_scope_qc_audit, contextualise, export_support, grouping_reasons,
    qc_table, qc_trail, read_rows, validate_rows, SCOPE_COLUMNS,
)

ORDERS = {"isotype_group": ["all_IGH", "IgM", "IgG", "IgA"],
          "v_region": ["whole_v", "fwr1", "cdr1", "fwr2", "cdr2", "fwr3"]}
MANIFEST = {"0460": dict(scope_id="s", library_id="0460", donor_id="004", tissue_id="007")}


def population(stage, index, locus, counts, *, analysis, kind, unit="rearrangement",
               output_unit=None, basis="annotated", isotype="", region=""):
    """One library audit row with units and support (in, out, in, out)."""
    row = checkpoint("0460", stage, "0", "0", locus=locus, analysis=analysis, index=index,
                     unit=unit, isotype=isotype, region=region)
    n_in, n_out, support_in, support_out = (str(value) for value in counts)
    row.update(kind=kind, locus_basis=basis, output_unit=output_unit or unit,
               n_units_in=n_in, n_units_out=n_out, consensus_count_in=support_in,
               consensus_count_out=support_out,
               support_conserved="TRUE" if kind == "grouping" else "FALSE")
    return row


def library_rows():
    """A library with IGH and TRB records, empty at the other loci.

    IGH: 100 annotated support, 90 passing productive and common QC, 60 grouping
    eligible, 81 class-callable and 45 subclass-callable. TRB: 40 annotated and
    common-QC support, 30 grouping eligible.
    """
    counts = {"IGH": ((10, 8, 100, 90), (8, 6, 90, 60)), "TRB": ((4, 4, 40, 40), (4, 3, 40, 30))}
    rows = [population("trust4_reported_airr", 3, "IGH", (10, 10, 100, 100), analysis="common",
                       kind="checkpoint", basis="native")]
    rows[0]["support_conserved"] = "TRUE"
    for locus in LOCI:
        common, grouping = counts.get(locus, ((0, 0, 0, 0), (0, 0, 0, 0)))
        members = (grouping[1], max(grouping[1] - 1, 0), grouping[3], grouping[3])
        rows += [population("qc_productive", 14, locus, common, analysis="common", kind="criterion"),
                 population("qc_pass_total", 20, locus, common, analysis="common", kind="eligibility"),
                 population("grouping_eligible", 30, locus, grouping, analysis="grouping",
                            kind="eligibility"),
                 population("exact_sequence_collapse", 32, locus, members, analysis="grouping",
                            kind="grouping", unit="backbone_representative",
                            output_unit="inference_representative"),
                 population("inference_representatives", 33, locus, (members[1],) * 2 + members[2:],
                            analysis="grouping", kind="checkpoint", unit="inference_representative")]
        unit = "library_clone_observation" if locus == "IGH" else "clonotype"
        groups = max(members[1] - 1, 0)
        rows.append(population("clone_observations" if locus == "IGH" else "exact_groups", 34, locus,
                               (groups, groups, members[3], members[3]), analysis="grouping",
                               kind="checkpoint", unit=unit))
    for level, index, callable_support in (("class", 46, 81), ("subclass", 47, 45)):
        rows.append(population(f"isotype_{level}_callable", index, "IGH", (8, 7, 90, callable_support),
                               analysis="isotype", kind="eligibility"))
    for group in ("IgM", "all_IGH"):
        rows.append(population("shm_eligible_members", 60, "IGH", (5, 2, 60, 30), analysis="shm",
                               kind="eligibility", unit="clone_member", isotype=group, region="whole_v"))
    return rows


def write_export(path, records):
    """Write an AIRR export holding (locus, consensus_count, grouping_status) records."""
    with gzip.open(path, "wt", newline="") as handle:
        handle.write("sequence_id\tlocus\tconsensus_count\tgrouping_status\n")
        for index, record in enumerate(records):
            handle.write("\t".join(("0460~" + str(index), *map(str, record))) + "\n")


# IGH: 60 eligible, 20 with an ambiguous V call, 5 with an ineligible junction
# and ambiguous J, and 5 blocked by their contig's dominant variant. TRB: 30
# eligible and 10 with a constant-region locus conflict.
EXPORT_RECORDS = [("IGH", 20, "eligible"), ("IGH", 40, "eligible"), ("IGH", 20, "v_call_ambiguous"),
                  ("IGH", 5, "ineligible_junction;j_call_ambiguous"),
                  ("IGH", 5, "ineligible_contig_representative"),
                  ("TRB", 30, "eligible"), ("TRB", 10, "c_locus_conflict")]


class ScopeAuditTests(unittest.TestCase):
    def grouping(self):
        rows = []
        for stage, stage_index, unit in (("grouping_eligible", 30, "rearrangement"),
                                         ("inference_representatives", 33, "inference_representative")):
            for index, locus in enumerate(LOCI):
                row = checkpoint("0460", stage, str(index), str(index * 10),
                                 locus=locus, analysis="grouping", index=stage_index, unit=unit)
                if stage == "grouping_eligible":
                    row.update(kind="eligibility", support_conserved="FALSE")
                rows.append(row)
        return rows

    def exact_grouping(self):
        rows = []
        for index, locus in enumerate(LOCI):
            row = checkpoint("0460", "exact_sequence_collapse", str(index), str(index * 10),
                             locus=locus, analysis="grouping", index=32,
                             unit="backbone_representative")
            row.update(kind="grouping", output_unit="inference_representative",
                       n_units_in=str(index * 2))
            rows.append(row)
        return rows

    def test_exact_sequence_totals_preserve_both_units_and_stage_context(self):
        rows = self.exact_grouping()
        # Backbone collapse applies only to BCR loci and has no seven-locus total.
        for locus in ("IGH", "IGK", "IGL"):
            row = checkpoint("0460", "assembly_backbone_collapse", "2", "10",
                             locus=locus, analysis="grouping", index=31, unit="rearrangement")
            row.update(kind="grouping", output_unit="backbone_representative")
            rows.append(row)
        canonical_totals(rows)
        totals = [row for row in rows if row["locus"] == "all"]
        self.assertEqual(len(totals), 1)
        total = totals[0]
        self.assertEqual((total["stage"], total["stage_index"], total["kind"]),
                         ("exact_sequence_collapse", "32", "grouping"))
        self.assertEqual((total["input_unit"], total["output_unit"]),
                         ("backbone_representative", "inference_representative"))
        self.assertEqual((total["n_units_in"], total["n_units_out"]), ("42", "21"))
        self.assertEqual((total["consensus_count_in"], total["consensus_count_out"]), ("210", "210"))
        validate_rows(rows)

    def test_exact_sequence_totals_reject_missing_loci_and_keep_explicit_zeros(self):
        rows = self.exact_grouping()
        with self.assertRaisesRegex(ValueError, "complete canonical-locus partition"):
            canonical_totals(rows[:-1])
        for row in rows:
            row.update(n_units_in="0", n_units_out="0", consensus_count_in="0", consensus_count_out="0")
        canonical_totals(rows)
        total, = [row for row in rows if row["locus"] == "all"]
        self.assertEqual([total[field] for field in ("n_units_in", "n_units_out", "consensus_count_in", "consensus_count_out")],
                         ["0", "0", "0", "0"])

    def test_split_totals_require_every_locus_and_preserve_zeros(self):
        rows = self.grouping()
        with self.assertRaisesRegex(ValueError, "canonical"):
            canonical_totals(rows[:-1])
        canonical_totals(rows)
        totals = [row for row in rows if row["locus"] == "all"]
        self.assertEqual(len(totals), 2)
        self.assertTrue(all(row["n_units_out"] == "21" for row in totals))
        self.assertTrue(all(row["consensus_count_out"] == "210" for row in totals))
        validate_rows(rows)

    def test_conflicting_totals_and_duplicate_context_fail(self):
        rows = self.grouping()
        canonical_totals(rows)
        rows[-1]["n_units_out"] = "22"
        with self.assertRaisesRegex(ValueError, "total"):
            canonical_totals(rows)
        rows = self.grouping()
        with self.assertRaisesRegex(ValueError, "Duplicate"):
            validate_rows(rows + rows[:1])

    def test_context_keeps_identifier_strings_and_checks_donor_grain(self):
        manifest = {"0460": dict(scope_id="NA", donor_id="004", tissue_id="007")}
        row = self.grouping()[0]
        self.assertEqual(contextualise(row, manifest, "NA")["donor_id"], "004")
        row.update(grain="donor", library_id="", donor_id="004")
        self.assertEqual(contextualise(row, manifest, "NA")["library_id"], "")
        row["library_id"] = "0460"
        with self.assertRaisesRegex(ValueError, "donor"):
            contextualise(row, manifest, "NA")

    def test_units_and_support_are_validated_independently(self):
        row = self.grouping()[0]
        row.update(n_units_in="1", n_units_out="1", consensus_count_in="10",
                   consensus_count_out="9", support_conserved="TRUE")
        with self.assertRaisesRegex(ValueError, "support"):
            validate_rows([row])
        row.update(consensus_count_out="10", n_units_in="0.5", n_units_out="0.5")
        with self.assertRaisesRegex(ValueError, "integer"):
            validate_rows([row])

    def collect(self, root, shard_rows, records=EXPORT_RECORDS):
        manifest = root / "manifest.tsv"
        manifest.write_text("scope_id\tlibrary_id\tdonor_id\ttissue_id\ns\t0460\t004\t007\n")
        shard = root / "shard.tsv"
        with shard.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, COLUMNS, delimiter="\t")
            writer.writeheader()
            writer.writerows(shard_rows)
        selection = root / "selection.tsv"
        selection.write_text("library_id\tisotype_group\tn_clone_consensus_sequences\tn_clone_consensus_sequences_used_cdr\tn_clone_consensus_sequences_used_fwr\n")
        export = root / ("0460" + "_rearrangements.tsv.gz")
        write_export(export, records)
        outputs = {name: root / f"{name}.tsv" for name in ("audit", "qc", "trail")}
        collect_scope_qc_audit(manifest, [shard], [], selection, outputs["audit"],
                               exports={"0460": export}, qc_output=outputs["qc"],
                               trail_output=outputs["trail"], value_orders=ORDERS)
        return outputs

    def test_collector_writes_every_table_with_manifest_context(self):
        with tempfile.TemporaryDirectory() as directory:
            outputs = self.collect(Path(directory), library_rows())
            # The audit leads with the identifiers; shards keep COLUMNS.
            self.assertEqual(outputs["audit"].read_text().splitlines()[0].split("\t"), list(SCOPE_COLUMNS))
            self.assertEqual(set(SCOPE_COLUMNS), set(COLUMNS))
            total, = [row for row in read_rows(outputs["audit"])
                      if row["stage"] == "exact_sequence_collapse" and row["locus"] == "all"]
            self.assertEqual((total["scope_id"], total["library_id"], total["stage_index"]),
                             ("s", "0460", "32"))
            self.assertEqual((total["n_units_in"], total["n_units_out"]), ("9", "7"))
            self.assertEqual((total["consensus_count_in"], total["consensus_count_out"]), ("90", "90"))
            qc = read_rows(outputs["qc"], QC_COLUMNS)
            self.assertEqual([row["locus"] for row in qc], list(LOCI))
            self.assertTrue(all(row["scope_id"] == "s" and row["donor_id"] == "004"
                                and row["tissue_id"] == "007" for row in qc))
            trail = read_rows(outputs["trail"], TRAIL_COLUMNS)
            self.assertTrue(all(row["scope_id"] == "s" and row["tissue_id"] == "007" for row in trail))
            # The canonical-locus totals the collector adds lead the library's trail.
            self.assertEqual([(row["locus"], row["stage"]) for row in trail[:2]],
                             [("all", "grouping_eligible"), ("all", "exact_sequence_collapse")])

    def test_failed_collection_preserves_existing_output(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            outputs = {name: root / f"{name}.tsv" for name in ("audit", "qc", "trail")}
            for path in outputs.values():
                path.write_bytes(b"previous snapshot\n")
            # Eligible export support must equal the audit's grouping_eligible output.
            records = [*EXPORT_RECORDS[:-2], ("TRB", 40, "eligible")]
            with self.assertRaisesRegex(ValueError, "AIRR export support of 0460 TRB"):
                self.collect(root, library_rows(), records)
            for path in outputs.values():
                self.assertEqual(path.read_bytes(), b"previous snapshot\n")
                self.assertFalse(path.with_suffix(".tsv.tmp").exists())
            with self.assertRaises(ValueError):
                self.collect(root, self.grouping()[:-1])
            self.assertEqual(outputs["audit"].read_bytes(), b"previous snapshot\n")


class QcTableTests(unittest.TestCase):
    def rows(self):
        rows = [contextualise(row, MANIFEST, "s") for row in library_rows()]
        validate_rows(rows)
        return rows

    def table(self, rows, records=EXPORT_RECORDS):
        with tempfile.TemporaryDirectory() as directory:
            export = Path(directory) / "0460_rearrangements.tsv.gz"
            write_export(export, records)
            return qc_table(rows, {"0460": export}, MANIFEST, "s")

    def test_reasons_follow_grouping_status_tokens(self):
        self.assertEqual(grouping_reasons("eligible"), frozenset())
        self.assertEqual(grouping_reasons("ineligible_junction;v_call_incomplete_candidates;j_call_missing"),
                         {"junction_ineligible", "v_call_unresolved", "j_call_unresolved"})
        self.assertEqual(grouping_reasons("c_locus_conflict;c_call_incomplete_candidates"),
                         {"c_call_ineligible"})
        self.assertEqual(grouping_reasons("ineligible_contig_representative"),
                         {"contig_representative_ineligible"})
        for status in ("", "v_call_", "eligible;v_call_ambiguous", "other_reason"):
            with self.subTest(status=status), self.assertRaisesRegex(ValueError, "Unknown grouping_status"):
                grouping_reasons(status)

    def test_export_support_sums_each_locus_and_reason(self):
        with tempfile.TemporaryDirectory() as directory:
            export = Path(directory) / "0460_rearrangements.tsv.gz"
            write_export(export, EXPORT_RECORDS)
            sums = export_support(export)
            write_export(export, [("TRX", 1, "eligible")])
            with self.assertRaisesRegex(ValueError, "locus 'TRX'"):
                export_support(export)
        self.assertEqual(set(sums), set(LOCI))
        self.assertEqual({key: int(value) for key, value in sums["IGH"].items()},
                         {"common_qc": 90, "grouping_eligible": 60, "v_call_unresolved": 20,
                          "j_call_unresolved": 5, "junction_ineligible": 5, "c_call_ineligible": 0,
                          "contig_representative_ineligible": 5})
        self.assertEqual(int(sums["TRB"]["c_call_ineligible"]), 10)
        self.assertTrue(all(value == 0 for value in sums["TRA"].values()))

    def test_fractions_use_their_row_denominators_and_empty_zero_denominators(self):
        table = self.table(self.rows())
        self.assertEqual([list(row) for row in table[:1]], [list(QC_COLUMNS)])
        igh, igk, _, _, trb, *_ = table
        self.assertEqual({column: igh[column] for column in QC_COLUMNS[7:QC_COLUMNS.index("n_grouped_igh_members")]}, {
            "igblast_consensus_count": "100", "nonproductive_fraction": "0.1",
            "common_qc_consensus_count": "90", "common_qc_fraction_of_annotated": "0.9",
            "grouped_consensus_count": "60", "grouping_eligible_fraction_of_common_qc": "0.666666666666667",
            "v_call_unresolved_fraction": "0.222222222222222", "j_call_unresolved_fraction": "0.0555555555555556",
            "junction_ineligible_fraction": "0.0555555555555556", "c_call_ineligible_fraction": "0",
            "contig_representative_ineligible_fraction": "0.0555555555555556",
            "class_callable_fraction": "0.9", "subclass_callable_fraction": "0.5"})
        self.assertEqual((trb["c_call_ineligible_fraction"], trb["class_callable_fraction"],
                          trb["nonproductive_fraction"]), ("0.25", "", "0"))
        self.assertEqual(igh["trust4_record_count"], "10")
        self.assertEqual(igh["trust4_consensus_count"], "100")
        # An observed empty locus has zero counts and no fraction.
        self.assertEqual([igk[column] for column in ("igblast_consensus_count", "common_qc_consensus_count",
                                                     "grouped_consensus_count")], ["0", "0", "0"])
        self.assertTrue(all(igk[column] == "" for column in QC_COLUMNS if column.endswith("_fraction")))

    def test_eligible_export_support_must_equal_the_audit(self):
        for records in ([*EXPORT_RECORDS[:-1], ("TRB", 10, "eligible")],
                        [*EXPORT_RECORDS, ("IGK", 1, "eligible")]):
            with self.subTest(records=records[-1]), self.assertRaisesRegex(ValueError, "AIRR export support"):
                self.table(self.rows(), records)

    def test_populations_must_follow_one_another(self):
        rows = self.rows()
        grouping = next(row for row in rows if row["stage"] == "grouping_eligible" and row["locus"] == "IGH")
        grouping["consensus_count_in"] = "100"
        with self.assertRaisesRegex(ValueError, "do not follow"):
            self.table(rows)
        rows = [row for row in self.rows() if row["stage"] != "isotype_subclass_callable"]
        with self.assertRaisesRegex(ValueError, "isotype_subclass_callable"):
            self.table(rows)


class QcTrailTests(unittest.TestCase):
    def rows(self):
        rows = [contextualise(row, MANIFEST, "s") for row in library_rows()]
        donor = population("donor_clones", 35, "IGH", (3, 3, 60, 60), analysis="grouping",
                           kind="checkpoint", unit="donor_clone")
        donor.update(grain="donor", library_id="", tissue_id="", donor_id="004", scope_id="s")
        return [*rows, donor]

    def test_release_view_renames_reorders_and_drops_repeated_rows(self):
        trail = qc_trail(self.rows(), ORDERS)
        self.assertEqual(list(trail[0]), list(TRAIL_COLUMNS))
        stages = {row["stage"] for row in trail}
        self.assertFalse(stages & {"inference_representatives", "donor_clones", "qc_pass_total"})
        self.assertTrue(all(row["library_id"] == "0460" for row in trail))
        igh = [(row["stage_order"], row["stage"], row["isotype_group"]) for row in trail if row["locus"] == "IGH"]
        self.assertEqual(igh, [("3", "trust4_records", ""), ("11", "qc_productive", ""), ("14", "common_qc", ""),
                               ("15", "grouping_eligible", ""), ("17", "exact_sequence_collapse", ""),
                               ("18", "receptor_units", ""), ("24", "class_callable", ""),
                               ("25", "subclass_callable", ""), ("28", "shm_eligible_members", "all_IGH"),
                               ("28", "shm_eligible_members", "IgM")])
        self.assertEqual([row["locus"] for row in trail][:12], ["IGH"] * 10 + ["IGK"] * 2)
        first, units = trail[0], next(row for row in trail if row["stage"] == "receptor_units"
                                       and row["locus"] == "IGH")
        self.assertEqual((first["locus_source"], first["input_unit"], first["kind"], first["analysis"]),
                         ("trust4", "airr_record", "checkpoint", "assembly"))
        self.assertEqual({row["locus_source"] for row in trail[1:]}, {"igblast"})
        # Receptor units start from the exact-sequence members and keep their support.
        self.assertEqual([units[column] for column in ("input_unit", "n_units_in", "output_unit", "n_units_out",
                                                       "consensus_count_in", "consensus_count_out")],
                         ["member", "5", "clone", "4", "60", "60"])
        collapse = next(row for row in trail if row["stage"] == "exact_sequence_collapse"
                        and row["locus"] == "IGH")
        self.assertEqual((collapse["input_unit"], collapse["output_unit"]), ("contig", "member"))
        shm = next(row for row in trail if row["stage"] == "shm_eligible_members")
        self.assertEqual((shm["input_unit"], shm["kind"]), ("member", "eligibility"))

    def test_stage_order_numbers_each_released_stage_once(self):
        names = list(dict.fromkeys(name for _, name, _, _ in TRAIL_STAGES))
        self.assertEqual(len(names), 29)
        self.assertEqual(sum(stage[1] == "receptor_units" for stage in TRAIL_STAGES), 2)

    def test_unknown_stages_units_strata_and_unsupported_groups_fail(self):
        for change, message in (({"stage": "new_population"}, "no released name"),
                                ({"output_unit": "effective_clone"}, "no released unit"),
                                ({"isotype_group": "IgE"}, "no listed order")):
            rows = self.rows()
            next(row for row in rows if row["stage"] == "shm_eligible_members").update(change)
            with self.subTest(change=change), self.assertRaisesRegex(ValueError, message):
                qc_trail(rows, ORDERS)
        rows = self.rows()
        next(row for row in rows if row["stage"] == "clone_observations").update(consensus_count_in="59")
        with self.assertRaisesRegex(ValueError, "Receptor units"):
            qc_trail(rows, ORDERS)

if __name__ == "__main__":
    unittest.main()
