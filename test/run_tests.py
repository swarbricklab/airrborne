#!/usr/bin/env python3
"""Run scientific component tests by pipeline component and test level.

Use --component shm --level unit for a focused run. R uses the pinned analysis
image; Python and JavaScript use the supplemental test environment. System
tests have their own command in test/system/run.sh.
"""

import argparse
import os
from pathlib import Path
import subprocess
import sys
import unittest

ROOT = Path(__file__).resolve().parents[1]
COMPONENTS = ROOT / "test/components"


def cases(suite: unittest.TestSuite):
    for item in suite:
        if isinstance(item, unittest.TestSuite):
            yield from cases(item)
        else:
            yield item


def run_python(directories: list[Path]) -> bool:
    """Run discovered modules and reject imports, omissions and required skips."""
    loader = unittest.TestLoader()
    suite = unittest.TestSuite()
    expected = set()
    for directory in directories:
        modules = list(directory.glob("test_*.py"))
        if not modules:
            continue
        expected.update(path.resolve() for path in modules)
        for path in sorted(modules):
            name = ".".join(path.relative_to(ROOT).with_suffix("").parts)
            suite.addTests(loader.loadTestsFromName(name))
    if not expected:
        return True
    discovered = {Path(sys.modules[type(case).__module__].__file__).resolve()
                  for case in cases(suite) if type(case).__module__ in sys.modules}
    missing = expected - discovered
    if loader.errors or missing:
        for error in loader.errors:
            print(error, file=sys.stderr)
        for path in sorted(missing):
            print(f"Undiscovered Python tests: {path.relative_to(ROOT)}", file=sys.stderr)
        return False
    print(f"Python: {len(expected)} modules, {suite.countTestCases()} cases", flush=True)
    result = unittest.TextTestRunner(verbosity=2).run(suite)
    if result.skipped:
        print("Required Python checks were skipped.", file=sys.stderr)
    return result.wasSuccessful() and not result.skipped


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--component", choices=sorted(p.name for p in COMPONENTS.iterdir()
                                                    if p.is_dir() and not p.name.startswith("_")))
    parser.add_argument("--level", choices=("unit", "integration", "all"), default="all")
    args = parser.parse_args()
    os.chdir(ROOT)
    sys.path.insert(0, str(ROOT))
    os.environ["AIRRBORNE_ROOT"] = str(ROOT)
    os.environ["RSCRIPT"] = str(ROOT / "test/support/rscript.sh")
    temporary = ROOT / ".snakemake/test-tmp"
    temporary.mkdir(parents=True, exist_ok=True)
    os.environ["TMPDIR"] = str(temporary)
    levels = ("unit", "integration") if args.level == "all" else (args.level,)
    components = [COMPONENTS / args.component] if args.component else sorted(
        path for path in COMPONENTS.iterdir() if path.is_dir() and not path.name.startswith("_"))
    directories = [component / level for component in components for level in levels
                   if (component / level).is_dir()]
    patterns = ("test_*.py", "test-*.R", "*.test.js", "test_*.sh")
    files = {pattern: sorted(path for directory in directories for path in directory.glob(pattern))
             for pattern in patterns}
    if not any(files.values()):
        parser.error("The selected component and level contain no tests.")
    for directory in directories:
        count = sum(len(list(directory.glob(pattern))) for pattern in patterns)
        if count:
            print(f"{directory.parent.name}/{directory.name}: {count} test files", flush=True)
    passed = run_python(directories)
    commands = []
    if files["test-*.R"]:
        commands.append([str(ROOT / "test/support/rscript.sh"),
                         str(ROOT / "test/support/run_r_tests.R"),
                         *(str(path) for path in files["test-*.R"])])
    if files["*.test.js"]:
        commands.append(["node", "--test", *(str(path) for path in files["*.test.js"])])
    commands.extend(["bash", str(path)] for path in files["test_*.sh"])
    for command in commands:
        try:
            passed = (subprocess.run(command, cwd=ROOT).returncode == 0) and passed
        except OSError as error:
            print(str(error), file=sys.stderr)
            passed = False
    print("Component tests passed." if passed else "Component tests failed.", flush=True)
    return 0 if passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
