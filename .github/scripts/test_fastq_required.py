"""The CI filter selects execution for native assembly and input changes."""

from pathlib import Path
import subprocess
import unittest


class FastqTriggerTests(unittest.TestCase):
    def test_actual_filter_for_native_inputs_and_unrelated_changes(self):
        root = Path(__file__).resolve().parents[2]
        for path, expected in (
                ("run_module.sh", "true"),
                ("workflow/scripts/obtain_image.py", "true"),
                ("workflow/scripts/record_assembly_provenance.pl", "true"),
                ("workflow/scripts/capture_assembly_context.py", "true"),
                ("workflow/lib/provenance.py", "true"),
                ("workflow/scripts/count_seq_depth.sh", "true"),
                ("workflow/lib/trust4_parsing.R", "true"),
                ("workflow/lib/table_io.R", "true"),
                ("test/data/reads/library/lane_R1.fastq.gz", "true"),
                ("test/data_builders/generate_receptors.py", "true"),
                ("test/system/run.sh", "true"),
                ("test/support/rscript.sh", "true"),
                ("test/run_tests.py", "true"),
                ("setup.sh", "true"),
                ("config/airrborne/config.yaml", "true"),
                ("test/test_env.sh", "true"),
                (".github/scripts/fastq_required.sh", "true"),
                ("docs/interpretation.md", "false"),
                ("workflow/report/report.js", "false"),
                ("workflow/lib/diversity.R", "false"),
                ("", "false")):
            with self.subTest(path=path):
                run = subprocess.run(["bash", str(root / ".github/scripts/fastq_required.sh")],
                                     input=path + "\n", capture_output=True, text=True, check=True)
                self.assertEqual(run.stdout.strip(), expected)


if __name__ == "__main__":
    unittest.main()
