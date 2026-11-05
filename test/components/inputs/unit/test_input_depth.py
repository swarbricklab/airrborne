"""Retained depth supplies one exact read-pair denominator per library.

Run: python test/run_tests.py --component inputs --level unit
"""

from pathlib import Path
import tempfile
import unittest

from workflow.scripts.normalise_input_depth import normalise_input_depth
from test.components.inputs.integration.reference_calculations import verify_input_depth


class InputDepthTests(unittest.TestCase):
    header = "gussid\tn_fastq_pairs\ttotal_read_pairs\ttotal_reads\ttotal_bases\n"

    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.directory = Path(self.temporary.name)
        self.source = self.directory / "native.tsv"
        self.output = self.directory / "analytical.tsv"

    def test_valid_depth_and_zero_depth_keep_literal_library_ids(self):
        for library, pairs, bases in (("0460", 10, 80), ("NA", 0, 0)):
            with self.subTest(library=library):
                self.source.write_text(self.header + f"{library}\t2\t{pairs}\t{pairs * 2}\t{bases}\n")
                normalise_input_depth(self.source, self.output, library)
                self.assertEqual(self.output.read_text(),
                                 "library_id\tn_fastq_file_pairs\tinput_read_pairs\ttotal_bases\n"
                                 f"{library}\t2\t{pairs}\t{bases}\n")

    def test_invalid_input_never_creates_or_replaces_an_output(self):
        valid = "0460\t1\t10\t20\t80\n"
        cases = {
            "contradictory totals": self.header + "0460\t1\t10\t7\t80\n",
            "wrong library": self.header + valid.replace("0460", "other"),
            "missing row": self.header,
            "extra row": self.header + valid * 2,
            "wrong schema": self.header.replace("gussid", "library_id") + valid,
            "short row": self.header + "0460\t1\t10\t20\n",
        }
        for token in ("", "NA", "Inf", "1.5", "-1", " 2", "1e2", "two", "\u00b2"):
            cases[f"invalid count {token!r}"] = self.header + f"0460\t1\t10\t20\t{token}\n"
        for name, source in cases.items():
            for existing in (False, True):
                with self.subTest(case=name, existing=existing):
                    self.source.write_text(source)
                    self.output.unlink(missing_ok=True)
                    if existing:
                        self.output.write_text("previous valid result\n")
                    with self.assertRaisesRegex(ValueError, "sequencing-depth|Sequencing-depth"):
                        normalise_input_depth(self.source, self.output, "0460")
                    if existing:
                        self.assertEqual(self.output.read_text(), "previous valid result\n")
                    else:
                        self.assertFalse(self.output.exists())
                    self.assertFalse(self.output.with_suffix(".tmp").exists())


class InputDepthVerifierTests(unittest.TestCase):
    def test_expected_libraries_and_native_invariant_are_required(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            directory = root / "libraries/0460/metrics"
            directory.mkdir(parents=True)
            source = directory / "0460_seq_depth.tsv"
            output = directory / "0460_input_depth.tsv"
            native = InputDepthTests.header + "0460\t1\t10\t20\t80\n"
            analytical = ("library_id\tn_fastq_file_pairs\tinput_read_pairs\ttotal_bases\n"
                          "0460\t1\t10\t80\n")
            source.write_text(native)
            output.write_text(analytical)
            # A selected scope needs only its own population, even when an
            # unrelated library has no depth products.
            self.assertEqual(verify_input_depth(root, {"0460"}), 1)
            with self.assertRaisesRegex(FileNotFoundError, "other_seq_depth"):
                verify_input_depth(root, {"0460", "other"})
            for path in (source, output):
                content = path.read_text()
                path.unlink()
                with self.assertRaises(FileNotFoundError):
                    verify_input_depth(root, {"0460"})
                path.write_text(content)
            for content, message in (
                    (native.replace("\t20\t", "\t7\t"), "total_reads"),
                    (native.replace("gussid", "library_id"), "schema"),
                    (native.replace("total_bases\n", "total_bases\ttotal_bases\n")
                     .replace("\t80\n", "\t80\t80\n"), "schema"),
                    (native.replace("0460", "460"), "library"),
                    (native.replace("\t80\n", "\tNA\n"), "integer")):
                source.write_text(content)
                with self.assertRaisesRegex(AssertionError, message):
                    verify_input_depth(root, {"0460"})
            source.write_text(native)
            output.write_text(analytical.replace("input_read_pairs", "total_read_pairs"))
            with self.assertRaisesRegex(AssertionError, "schema"):
                verify_input_depth(root, {"0460"})


if __name__ == "__main__":
    unittest.main()
