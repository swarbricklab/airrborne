# Profiles

Profiles are configuration files for where jobs run and how many resources they request. `run_module.sh` selects the appropriate profiles. For a local first run, the bundled settings are ready to use.

On Gadi, the shared profile stores the release environment in the maintainers' project, a56. For another project, uncomment `conda-prefix` in `workflow/config.yaml` to use a directory in the clone. Leave the shared setting in `global/nci/config.yaml` unchanged. The [Gadi instructions](../docs/running.md#gadi) explain the setup.

| Directory | Use |
|---|---|
| `local/` | Runs every job on the launching host or inside one scheduler allocation. It holds the per-rule memory and time values for that case and an 8-core default that `--cores N` replaces. |
| `workflow/` | Requests cores, memory, walltime and disk for each rule when a site profile submits jobs to a scheduler. Its commented `conda-prefix` line keeps the release environment in the launch directory instead of the Gadi shared cache. |
| `global/` | A Git submodule shared by the lab's workflows. Its `nci/` directory is the PBS Pro site profile for NCI Gadi. Its READMEs cover the shared adapter across the lab's workflows. Follow this pipeline's running guide for airrborne setup. |

The [running guide](../docs/running.md#run-on-a-cluster) explains when each applies and how to write a site profile for another scheduler.
