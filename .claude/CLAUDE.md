# airrborne

A Snakemake workflow for human B-cell and T-cell receptor analysis from bulk
RNA-seq. Each scope produces a validated, portable `release/` directory.
Paths below are relative to the repository root.

## Start here

- `README.md` and `docs/running.md`: setup, inputs and launch commands.
- `docs/interpretation.md`, `docs/qc_matrix.md` and `docs/data_dictionary.csv`:
  scientific definitions, eligibility and output contracts.
- `docs/maintenance.md`: dependencies, provenance, publication and deployment.
- `test/README.md` and `test/COVERAGE.md`: test commands and scientific coverage.
- `test/data/README.md`: synthetic inputs and their regeneration procedure.

## Code and execution

- `Snakefile` and `workflow/rules/` define orchestration; `workflow/scripts/`
  contains entry points; `workflow/lib/` holds shared R and Python code.
- `workflow/report/` and `workflow/release/` contain report assets and templates.
- Launch through `setup.sh` and `run_module.sh`. Paths in configuration resolve
  from the launch directory, including when used as a submodule.
- Use the pinned Snakemake environment, container images and release environment.
  R tests run in the analysis image; Python tests use the supplemental test env.
- Declare consumed code and settings as workflow dependencies. Check rerun scope
  when changing imports or shared modules; never hide dependencies to avoid reruns.

## Scientific and release constraints

- Preserve analytical decisions during refactoring. Changes to methods or cutoffs
  need corresponding scientific tests and documentation.
- `consensus_count` is recovered read support, not a cell or molecule count.
  These data provide no receptor-chain pairing. Missing estimates are not zero.
- Preserve string identifiers, original calls, row units and denominators.
  IGH clone assignment is specific to donor and scope; transfer support once.
- Keep common QC distinct from eligibility for grouping, isotype, SHM and selection.
  Reports and exports must read the same decisions used by the analysis.
- Keep release publication transactional: validate staged files before publishing,
  and declare only `release/.complete.json` as the publisher's Snakemake output.
- Replace source tables atomically; released tables may be hard links to them.
  Preserve processing identities, effective settings and provenance.
- For table changes, update writers, readers, the data dictionary, release allow-list
  and interpretation guide together. Edit the dictionary with Python's CSV module.
- Update configuration schemas, profiles and guides when their interfaces change.

## Verification

Prepare the test environment from the repository root:

```bash
./setup.sh --config test/system/config.yaml --tests
source .snakemake/env.sh
python test/run_tests.py --component NAME --level unit
```

Choose the affected component and level from `test/README.md`. Complete pipeline
and release acceptance uses its full component and fresh system-run commands,
including FASTQ, post-TRUST4, empty, permuted and copied offline-report checks.
Local tests do not establish Gadi/PBS acceptance.

Committed test inputs are synthetic. Normal tests must not regenerate them.
Preserve reviewed expectations and independent reference calculations; change
fixtures only through the documented construction and promotion procedure.
Keep historical runs, reference caches and images. Do not commit generated data,
logs, credentials or container images. Working notes under `docs/` are ignored;
new public guides must be added to the `.gitignore` allow-list.
