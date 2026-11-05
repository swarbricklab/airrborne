"""The handover preserves source assemblies while deduplicating inference DNA.

Run: python test/run_tests.py --component grouping --level unit
"""

from pathlib import Path
import tempfile
import unittest

from test.components.release.integration.check_release import (
    check_clone_isotype, check_grouping, check_released_member_shm, check_shm_evidence,
    check_shm_reconstruction, check_clonotype_members, check_shared_representatives, rows_by_library, type7_quantile,
    LOCI, SHM_GROUPS, SHM_REGIONS)


def source_rows(locus="IGH"):
    # A_0 and B_0 share exact DNA. B_1 is a distinct sibling whose support
    # reaches A_0 through B_0, the representative of its contig.
    rows = []
    for identity, dna, support, contig, is_contig in (
            ("A_0", "ACGT", 100, "A_0", True),
            ("B_0", "ACGT", 100, "B_0", True),
            ("B_1", "TGCA", 10, "B_0", False)):
        rows.append(dict(library_id="L", sequence_id="L~" + identity,
                         sequence=dna, locus=locus, consensus_count=str(support),
                         contig_id="L~" + identity.split("_")[0],
                         grouping_eligible="T", grouping_status="eligible",
                         contig_representative_sequence_id="L~" + contig,
                         is_contig_representative="T" if is_contig else "F",
                         member_representative_sequence_id="L~A_0",
                         is_member_representative="T" if identity == "A_0" else "F",
                         member_consensus_count="210" if identity == "A_0" else "0"))
    # The scope export names the member's clone or clonotype in clone_id.
    for row in rows:
        row["clone_id"] = "S~D~1" if row["is_member_representative"] == "T" else ""
    return rows


def reviewed_trail(locus="IGH", source_count=3, backbones=2, support=210):
    rows = []
    for chain in (*LOCI, "all"):
        n, b, final, total = ((source_count, backbones, 1, support)
                              if chain in (locus, "all") else (0, 0, 0, 0))
        stages = [("grouping_eligible", 15, n, n),
                  ("exact_sequence_collapse", 17, b, final)]
        if chain in ("IGH", "IGK", "IGL"):
            stages.append(("contig_collapse", 16, n, b))
        if chain != "all":
            stages.append(("receptor_units", 18, final, final))
        for stage, order, incoming, outgoing in stages:
            rows.append(dict(library_id="L", locus=chain, stage=stage, stage_order=str(order),
                             isotype_group="", v_region="",
                             n_units_in=str(incoming), n_units_out=str(outgoing),
                             consensus_count_in=str(total), consensus_count_out=str(total)))
    return rows


class GroupingHandoverTests(unittest.TestCase):
    def test_nonigh_membership_joins_representatives_and_conserves_support(self):
        rows = source_rows("TRB")
        rows[0].update(junction="TGTGCC", junction_length="6", v_gene_key="TRBV1", j_gene_key="TRBJ1-1")
        key = "~".join(("TGTGCC", "TRBV1", "TRBJ1-1"))
        member = dict(sequence_id="L~A_0", locus="TRB", clonotype_id="L_TRB_ct00001",
                      member_consensus_count="210")
        group = dict(locus="TRB", clonotype_id="L_TRB_ct00001", clonotype_key=key, consensus_count="210",
                     representative_sequence_id="L~A_0", junction_length="6")
        check_clonotype_members([member], [group], rows, "L")
        for altered in ([{**member, "member_consensus_count": "100"}], [], [member, member],
                        [{**member, "sequence_id": "B_1"}], [{**member, "clonotype_id": "other"}]):
            with self.subTest(members=altered), self.assertRaises(AssertionError):
                check_clonotype_members(altered, [group], rows, "L")
        # The group key must follow the member's source, and the representative
        # must be a member of its own group, with its junction length.
        for altered in ({**group, "clonotype_key": "other"},
                        {**group, "representative_sequence_id": "L~B_0"},
                        {**group, "junction_length": "9"}):
            with self.subTest(group=altered), self.assertRaises(AssertionError):
                check_clonotype_members([member], [altered], rows, "L")

    def test_shared_representative_names_its_member_record(self):
        rows = source_rows("TRB")
        rows[0]["junction"] = "TGTGCC"
        shared = dict(locus="TRB", representative_junction="TGTGCC", representative_sequence_id="L~A_0")
        check_shared_representatives([shared], rows, "L")
        for altered in ({**shared, "representative_sequence_id": "L~B_0"},
                        {**shared, "representative_junction": "TGTGCA"},
                        {**shared, "locus": "TRA"}):
            with self.subTest(shared=altered), self.assertRaises(AssertionError):
                check_shared_representatives([altered], rows, "L")

    def test_distinct_sibling_transfers_via_its_contig_before_exact_merge(self):
        check_grouping(source_rows(), "L", reviewed_trail())

    def test_broken_links_duplicate_sequences_and_double_support_fail(self):
        mutations = [
            (1, "contig_representative_sequence_id", "L~A_0"),
            (1, "sequence", "TGCA"),
            (2, "member_representative_sequence_id", "L~B_1"),
            (2, "is_contig_representative", "T"),
            (0, "member_consensus_count", "220"),
            (2, "member_consensus_count", "10"),
        ]
        for index, column, value in mutations:
            with self.subTest(column=column, value=value):
                rows = source_rows()
                rows[index][column] = value
                with self.assertRaises(AssertionError):
                    check_grouping(rows, "L", reviewed_trail())
        rows = source_rows()
        rows[0]["member_consensus_count"] = "100"
        rows[1].update(is_member_representative="T", member_consensus_count="110")
        for row in rows[1:]:
            row["member_representative_sequence_id"] = "L~B_0"
        with self.assertRaisesRegex(AssertionError, "repeat an exact sequence"):
            check_grouping(rows, "L", reviewed_trail())

    def test_tcr_retains_its_own_contig_row_before_exact_merge(self):
        rows = source_rows("TRB")[:2]
        rows[0]["member_consensus_count"] = "200"
        check_grouping(rows, "L", reviewed_trail("TRB", 2, 2, 200))
        rows[1]["contig_representative_sequence_id"] = "L~A_0"
        rows[1]["is_contig_representative"] = "F"
        with self.assertRaisesRegex(AssertionError, "TCR contig"):
            check_grouping(rows, "L", reviewed_trail("TRB", 2, 2, 200))

    def test_trail_counts_contigs_sequences_and_groups_separately(self):
        for stage, column, value in (("contig_collapse", "n_units_out", "1"),
                                      ("exact_sequence_collapse", "n_units_in", "3"),
                                      ("exact_sequence_collapse", "stage_order", "5"),
                                      ("receptor_units", "n_units_out", "2")):
            rows = reviewed_trail()
            row = next(row for row in rows if row["stage"] == stage and row["locus"] == "IGH")
            row[column] = value
            with self.subTest(stage=stage, column=column), self.assertRaises(AssertionError):
                check_grouping(source_rows(), "L", rows)

    def test_shm_links_retain_only_identical_contig_evidence(self):
        # A donor shard row keeps the evidence IDs and every measured rate.
        member = dict(library_id="L", sequence_id="L~A_0")
        for region in ("whole_v", "fwr1", "cdr1", "fwr2", "cdr2", "fwr3"):
            member[f"shm_eligible_{region}"] = "TRUE"
            member[f"shm_evidence_sequence_id_{region}"] = "L~B_0"
        for region in ("whole_v", "cdr2", "fwr3"):
            member[f"shm_frequency_{region}"] = "0.02"
            member[f"mutations_{region}"] = "6"
        check_shm_evidence([member], source_rows())
        for evidence in ("", "L~B_1", "other~B_0"):
            altered = {**member, "shm_evidence_sequence_id_whole_v": evidence}
            with self.subTest(evidence=evidence), self.assertRaises(AssertionError):
                check_shm_evidence([altered], source_rows())
        member.update(shm_eligible_whole_v="FALSE", shm_evidence_sequence_id_whole_v="")
        check_shm_evidence([member], source_rows())

    def test_released_member_rates_are_empty_outside_eligibility(self):
        member = dict(library_id="L", sequence_id="L~A_0")
        row = dict(sequence_id="L~A_0", v_region="whole_v", shm_eligible="TRUE", n_mutations="6",
                   shm_frequency="0.02", median_trust4_consensus_base_votes="4")
        check_released_member_shm([member], [row], source_rows())
        row.update(shm_eligible="FALSE")
        for column in ("n_mutations", "shm_frequency", "median_trust4_consensus_base_votes"):
            with self.subTest(filled=column), self.assertRaises(AssertionError):
                check_released_member_shm([member], [{**row, "n_mutations": "", "shm_frequency": "",
                                                      "median_trust4_consensus_base_votes": "", column: "1"}], source_rows())
        row.update(n_mutations="", shm_frequency="", median_trust4_consensus_base_votes="")
        check_released_member_shm([member], [row], source_rows())
        # A member and region appears once, and every row names a representative.
        with self.assertRaisesRegex(AssertionError, "Duplicate"):
            check_released_member_shm([member], [row, row], source_rows())
        with self.assertRaisesRegex(AssertionError, "inference representative"):
            check_released_member_shm([member], [{**row, "sequence_id": "L~B_0"}], source_rows())
        with self.assertRaisesRegex(AssertionError, "inference representative"):
            check_released_member_shm([{**member, "sequence_id": "L~B_0"}], [row], source_rows())

    def test_type_7_quantiles_interpolate_between_order_statistics(self):
        # R gives quantile(c(0.1, 0.2, 0.6, 0.5), 0.25) = 0.175 and 0.75 = 0.525.
        self.assertAlmostEqual(type7_quantile([0.1, 0.2, 0.6, 0.5], 0.25), 0.175)
        self.assertAlmostEqual(type7_quantile([0.1, 0.2, 0.6, 0.5], 0.75), 0.525)
        self.assertEqual(type7_quantile([0.4], 0.25), 0.4)


def shm_rebuild_inputs():
    """Rows of one library: an IgG1 and an IgM member share clone C1, and C2 holds a tied IgG call.

    In whole V, m1 has 0 mutations at rate 0, m2 has 3 at rate 0.01 and m3 has 6 at
    rate 0.02. m2 is measured but ineligible in FWR1 and m3 has no FWR1 rate. No
    other region has a member row. The summaries hold the values worked by hand.
    """
    members = [dict(sequence_id="L~m1", clone_id="C1", isotype_class="IgG", isotype_subclass="IgG1", member_consensus_count="5"),
               dict(sequence_id="L~m2", clone_id="C1", isotype_class="IgM", isotype_subclass="IgM", member_consensus_count="7"),
               dict(sequence_id="L~m3", clone_id="C2", isotype_class="IgG", isotype_subclass="", member_consensus_count="11")]
    def row(member, region, eligible, mutations, rate, votes):
        return dict(sequence_id=f"L~{member}", clone_id="C2" if member == "m3" else "C1", v_region=region,
                    shm_eligible=eligible, n_mutations=mutations, shm_frequency=rate, median_trust4_consensus_base_votes=votes)
    member_shm = [row("m1", "whole_v", "TRUE", "0", "0", "4"),
                  row("m2", "whole_v", "TRUE", "3", "0.01", "6"),
                  row("m3", "whole_v", "TRUE", "6", "0.02", "10"),
                  row("m1", "fwr1", "TRUE", "0", "0", "2"),
                  row("m2", "fwr1", "FALSE", "", "", "")]
    empty = dict(n_members_measured="0", n_members_eligible="0", n_members_zero_shm="0", n_clones_eligible="0",
                 n_clones_zero_median_shm="0", fraction_clones_zero_median_shm="", mean_clone_median_shm="",
                 median_clone_median_shm="", q25_clone_median_shm="", q75_clone_median_shm="",
                 iqr_clone_median_shm="", median_trust4_consensus_base_votes="")
    summaries = {(group, region): dict(empty, isotype_group=group, v_region=region)
                 for group in SHM_GROUPS for region in SHM_REGIONS}
    def filled(measured, eligible, zero, clones, zero_clones, fraction, mean, median, q25, q75, votes):
        return dict(n_members_measured=measured, n_members_eligible=eligible, n_members_zero_shm=zero,
                    n_clones_eligible=clones, n_clones_zero_median_shm=zero_clones,
                    fraction_clones_zero_median_shm=fraction, mean_clone_median_shm=mean,
                    median_clone_median_shm=median, q25_clone_median_shm=q25, q75_clone_median_shm=q75,
                    iqr_clone_median_shm=str(float(q75) - float(q25)), median_trust4_consensus_base_votes=votes)
    # C1's median is 0.005 and C2's 0.02, so their mean and median are 0.0125.
    summaries[("all_IGH", "whole_v")].update(filled("3", "3", "1", "2", "0", "0", "0.0125", "0.0125", "0.00875", "0.01625", "6"))
    # IgG holds m1 alone in C1 and m3 in C2, so C1's median falls to zero.
    summaries[("IgG", "whole_v")].update(filled("2", "2", "1", "2", "1", "0.5", "0.01", "0.01", "0.005", "0.015", "7"))
    summaries[("IgM", "whole_v")].update(filled("1", "1", "0", "1", "0", "0", "0.01", "0.01", "0.01", "0.01", "6"))
    summaries[("IgG1", "whole_v")].update(filled("1", "1", "1", "1", "1", "1", "0", "0", "0", "0", "4"))
    summaries[("all_IGH", "fwr1")].update(filled("2", "1", "1", "1", "1", "1", "0", "0", "0", "0", "2"))
    summaries[("IgG", "fwr1")].update(filled("1", "1", "1", "1", "1", "1", "0", "0", "0", "0", "2"))
    summaries[("IgM", "fwr1")].update(n_members_measured="1")
    summaries[("IgG1", "fwr1")].update(filled("1", "1", "1", "1", "1", "1", "0", "0", "0", "0", "2"))
    # Each subset's members and support, which include its ineligible members.
    subsets = {("C1", "all_IGH"): ("2", "12"), ("C2", "all_IGH"): ("1", "11"), ("C1", "IgG"): ("1", "5"),
               ("C2", "IgG"): ("1", "11"), ("C1", "IgM"): ("1", "7"), ("C1", "IgG1"): ("1", "5")}
    def clone(clone_id, group, region, eligible, zero, positive, median, state, support):
        n_members, consensus_count = subsets[(clone_id, group)]
        return dict(clone_id=clone_id, isotype_group=group, v_region=region, n_members=n_members,
                    n_members_eligible=eligible, n_members_zero_shm=zero, n_members_positive_shm=positive,
                    median_shm_frequency=median, mutation_state=state, consensus_count=consensus_count,
                    eligible_consensus_count=support)
    clones = [clone("C1", "all_IGH", "whole_v", "2", "1", "1", "0.005", "mixed", "12"),
              clone("C2", "all_IGH", "whole_v", "1", "0", "1", "0.02", "all_positive", "11"),
              clone("C1", "IgG", "whole_v", "1", "1", "0", "0", "all_zero", "5"),
              clone("C2", "IgG", "whole_v", "1", "0", "1", "0.02", "all_positive", "11"),
              clone("C1", "IgM", "whole_v", "1", "0", "1", "0.01", "all_positive", "7"),
              clone("C1", "IgG1", "whole_v", "1", "1", "0", "0", "all_zero", "5"),
              clone("C1", "all_IGH", "fwr1", "1", "1", "0", "0", "all_zero", "5"),
              clone("C1", "IgG", "fwr1", "1", "1", "0", "0", "all_zero", "5"),
              clone("C1", "IgG1", "fwr1", "1", "1", "0", "0", "all_zero", "5")]
    return members, member_shm, list(summaries.values()), clones


class CloneIsotypeTests(unittest.TestCase):
    def inputs(self):
        members = [dict(clone_id="C1", isotype_class="IgG", isotype_subclass="IgG1", member_consensus_count="3"),
                   dict(clone_id="C1", isotype_class="IgG", isotype_subclass="", member_consensus_count="5"),
                   dict(clone_id="C1", isotype_class="IgM", isotype_subclass="IgM", member_consensus_count="2"),
                   dict(clone_id="C2", isotype_class="", isotype_subclass="", member_consensus_count="9")]
        def row(level, isotype, n_members, support, total):
            return dict(clone_id="C1", isotype_level=level, isotype=isotype, n_members=str(n_members),
                        consensus_count=str(support), callable_consensus_count=str(total),
                        support_fraction=str(support / total))
        rows = [row("class", "IgM", 1, 2, 10), row("class", "IgG", 2, 8, 10),
                row("subclass", "IgM", 1, 2, 5), row("subclass", "IgG1", 1, 3, 5)]
        return rows, members

    def test_rows_sum_to_the_members_within_each_clone_and_level(self):
        check_clone_isotype(*self.inputs(), "L")
        check_clone_isotype([], [], "L")

    def test_a_changed_cell_or_row_set_fails(self):
        for column, value in (("n_members", "9"), ("consensus_count", "9"),
                              ("callable_consensus_count", "9"), ("support_fraction", "0.9")):
            rows, members = self.inputs()
            rows[1][column] = value
            with self.subTest(column=column), self.assertRaisesRegex(AssertionError, column):
                check_clone_isotype(rows, members, "L")
        rows, members = self.inputs()
        with self.assertRaisesRegex(AssertionError, "differ from the resolved members"):
            check_clone_isotype(rows[1:], members, "L")
        with self.assertRaisesRegex(AssertionError, "differ from the resolved members"):
            check_clone_isotype([*rows, {**rows[0], "clone_id": "C2"}], members, "L")
        with self.assertRaisesRegex(AssertionError, "Duplicate"):
            check_clone_isotype([*rows, rows[0]], members, "L")


class LibraryStreamTests(unittest.TestCase):
    def table(self, text):
        directory = tempfile.TemporaryDirectory()
        self.addCleanup(directory.cleanup)
        path = Path(directory.name) / "table.tsv"
        path.write_text(text)
        return path

    def test_each_library_gets_its_rows_and_an_empty_library_none(self):
        path = self.table("library_id\tvalue\nA\t1\nA\t2\nC\t3\n")
        observed = [(library, [row["value"] for row in rows])
                    for library, rows in rows_by_library(path, {"C", "B", "A"})]
        self.assertEqual(observed, [("A", ["1", "2"]), ("B", []), ("C", ["3"])])
        self.assertEqual(list(rows_by_library(self.table("library_id\tvalue\n"), ["A"])), [("A", [])])

    def test_an_unsorted_or_outside_library_stops_the_stream(self):
        for text in ("library_id\tvalue\nB\t1\nA\t2\n", "library_id\tvalue\nZ\t1\n"):
            with self.subTest(text=text), self.assertRaisesRegex(AssertionError, "unsorted or names another"):
                list(rows_by_library(self.table(text), ["A", "B"]))


class ShmReconstructionTests(unittest.TestCase):
    def test_member_rows_rebuild_both_summaries(self):
        check_shm_reconstruction(*shm_rebuild_inputs(), "L")

    def test_each_released_summary_cell_is_compared(self):
        members, member_shm, summaries, clones = shm_rebuild_inputs()
        whole = next(row for row in summaries if (row["isotype_group"], row["v_region"]) == ("all_IGH", "whole_v"))
        for column in ("n_members_measured", "n_members_eligible", "n_members_zero_shm", "n_clones_eligible",
                       "n_clones_zero_median_shm", "fraction_clones_zero_median_shm", "mean_clone_median_shm",
                       "median_clone_median_shm", "q25_clone_median_shm", "q75_clone_median_shm",
                       "iqr_clone_median_shm", "median_trust4_consensus_base_votes"):
            original = whole[column]
            whole[column] = str(float(original) + 1)
            with self.subTest(table="shm", column=column), self.assertRaisesRegex(AssertionError, column):
                check_shm_reconstruction(members, member_shm, summaries, clones, "L")
            whole[column] = original
        for column, value in (("n_members", "9"), ("n_members_eligible", "9"), ("n_members_zero_shm", "9"),
                              ("n_members_positive_shm", "9"), ("median_shm_frequency", "0.9"),
                              ("mutation_state", "all_zero"), ("consensus_count", "99"),
                              ("eligible_consensus_count", "99")):
            original = clones[0][column]
            clones[0][column] = value
            with self.subTest(table="clone_shm", column=column), self.assertRaisesRegex(AssertionError, column):
                check_shm_reconstruction(members, member_shm, summaries, clones, "L")
            clones[0][column] = original

    def test_values_within_tolerance_and_written_precision_pass(self):
        members, member_shm, summaries, clones = shm_rebuild_inputs()
        whole = next(row for row in summaries if (row["isotype_group"], row["v_region"]) == ("all_IGH", "whole_v"))
        whole["mean_clone_median_shm"] = "0.0125000000005"
        check_shm_reconstruction(members, member_shm, summaries, clones, "L")
        whole["mean_clone_median_shm"] = "0.012500002"
        with self.assertRaisesRegex(AssertionError, "mean_clone_median_shm"):
            check_shm_reconstruction(members, member_shm, summaries, clones, "L")

    def test_row_sets_must_match(self):
        members, member_shm, summaries, clones = shm_rebuild_inputs()
        with self.assertRaisesRegex(AssertionError, "lacks a clone"):
            check_shm_reconstruction(members, member_shm, summaries, clones[1:], "L")
        extra = {**clones[0], "clone_id": "C9"}
        with self.assertRaisesRegex(AssertionError, "no eligible member"):
            check_shm_reconstruction(members, member_shm, summaries, [*clones, extra], "L")
        with self.assertRaisesRegex(AssertionError, "complete group and region grid"):
            check_shm_reconstruction(members, member_shm, summaries[1:], clones, "L")
        with self.assertRaisesRegex(AssertionError, "names no clone member"):
            check_shm_reconstruction(members[1:], member_shm, summaries, clones, "L")

    def test_a_library_with_no_members_rebuilds_zeros_and_empty_cells(self):
        empty = dict(n_members_measured="0", n_members_eligible="0", n_members_zero_shm="0", n_clones_eligible="0",
                     n_clones_zero_median_shm="0", fraction_clones_zero_median_shm="", mean_clone_median_shm="",
                     median_clone_median_shm="", q25_clone_median_shm="", q75_clone_median_shm="",
                     iqr_clone_median_shm="", median_trust4_consensus_base_votes="")
        summaries = [dict(empty, isotype_group=group, v_region=region) for group in SHM_GROUPS for region in SHM_REGIONS]
        check_shm_reconstruction([], [], summaries, [], "L")


if __name__ == "__main__":
    unittest.main()
