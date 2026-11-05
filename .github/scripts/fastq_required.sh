#!/usr/bin/env bash
# Read changed repository paths on stdin, one per line, and print true when any
# of them needs fresh FASTQ acceptance, otherwise false. The patterns cover
# assembly, FASTQ discovery, native parsing, reference preparation, container
# execution, configuration, schemas, profiles, environments, fixture inputs, the
# fixture runner and this CI filter and workflow.
# unit-tests.yaml calls this filter for pull requests. Manual dispatch always
# selects FASTQ acceptance. Bash case wildcards also match nested path segments.
set -euo pipefail

run_fastq=false
while IFS= read -r path; do
    case "$path" in
        Snakefile|run_module.sh|workflow/scripts/obtain_image.py|workflow/rules/assembly.smk|workflow/rules/inputs.smk|workflow/rules/common.smk|workflow/rules/library.smk|\
        workflow/scripts/parse_trust4_airr.R|workflow/scripts/obtain_imgt_reference.py|workflow/lib/germline_reference.py|workflow/scripts/count_seq_depth.sh|\
        workflow/scripts/record_assembly_provenance.pl|workflow/scripts/capture_assembly_context.py|workflow/scripts/container_exec.sh|\
        workflow/lib/trust4_*|workflow/lib/runtime.R|workflow/lib/table_io.R|workflow/lib/airr_io.R|workflow/lib/qc_audit.R|workflow/lib/provenance.py|\
        config/*|workflow/schemas/*|envs/*|profiles/*|test/data/*|test/data_builders/*|test/system/*|test/support/*|test/run_tests.py|setup.sh|test/test_env.sh|\
        .github/scripts/fastq_required.sh|.github/workflows/unit-tests.yaml)
            run_fastq=true ;;
    esac
done
printf '%s\n' "$run_fastq"
