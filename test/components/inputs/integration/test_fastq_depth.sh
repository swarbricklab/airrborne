#!/usr/bin/env bash
# Check FASTQ structure and paired counts used as library depth denominators.
# Run: python test/run_tests.py --component inputs --level integration
# Run from the repository root. All constructed inputs and outputs live in one
# temporary directory, which this script removes on success or failure.
set -euo pipefail

test_dir="$(mktemp -d "${TMPDIR:-/tmp}/trust4-depth-test.XXXXXX")"
trap 'rm -rf "$test_dir"' EXIT

write_fastq() {
    destination="$1"
    records="$2"
    : > "${destination}.plain"
    for index in $(seq 1 "$records"); do
        printf '@read%s\nACGT\n+\nIIII\n' "$index" >> "${destination}.plain"
    done
    gzip -c "${destination}.plain" > "$destination"
}

write_fastq "$test_dir/sample_R1.fastq.gz" 2
write_fastq "$test_dir/sample_R2.fastq.gz" 2
# Two read pairs contain four reads. Each read has four bases.
bash workflow/scripts/count_seq_depth.sh \
    "$test_dir/depth.tsv" synthetic \
    --r1 "$test_dir/sample_R1.fastq.gz" \
    --r2 "$test_dir/sample_R2.fastq.gz" >/dev/null
awk -F '\t' 'NR == 2 && $3 == 2 && $4 == 4 { ok=1 } END { exit !ok }' \
    "$test_dir/depth.tsv"

# Multiple matched file pairs contribute to one library depth table.
write_fastq "$test_dir/extra_R1.fastq.gz" 1
write_fastq "$test_dir/extra_R2.fastq.gz" 1
bash workflow/scripts/count_seq_depth.sh \
    "$test_dir/multifile_depth.tsv" synthetic \
    --r1 "$test_dir/sample_R1.fastq.gz" "$test_dir/extra_R1.fastq.gz" \
    --r2 "$test_dir/sample_R2.fastq.gz" "$test_dir/extra_R2.fastq.gz" \
    >/dev/null
awk -F '\t' \
    'NR == 2 && $2 == 2 && $3 == 3 && $4 == 6 && $5 == 24 { ok=1 } END { exit !ok }' \
    "$test_dir/multifile_depth.tsv"

write_fastq "$test_dir/short_R2.fastq.gz" 1
if bash workflow/scripts/count_seq_depth.sh \
    "$test_dir/mismatch.tsv" synthetic \
    --r1 "$test_dir/sample_R1.fastq.gz" \
    --r2 "$test_dir/short_R2.fastq.gz" >/dev/null 2>&1; then
    echo "count_seq_depth accepted unequal mate record counts" >&2
    exit 1
fi
if [[ -e "$test_dir/mismatch.tsv" ]]; then
    echo "count_seq_depth wrote output after a mate-count failure" >&2
    exit 1
fi

# Equal library totals must not conceal lane/top-up pairs with unequal mates.
write_fastq "$test_dir/part_a_R1.fastq.gz" 1
write_fastq "$test_dir/part_a_R2.fastq.gz" 2
write_fastq "$test_dir/part_b_R1.fastq.gz" 2
write_fastq "$test_dir/part_b_R2.fastq.gz" 1
if bash workflow/scripts/count_seq_depth.sh \
    "$test_dir/compensating_mismatch.tsv" synthetic \
    --r1 "$test_dir/part_a_R1.fastq.gz" "$test_dir/part_b_R1.fastq.gz" \
    --r2 "$test_dir/part_a_R2.fastq.gz" "$test_dir/part_b_R2.fastq.gz" \
    >/dev/null 2>&1; then
    echo "count_seq_depth accepted per-pair mismatches with equal aggregate totals" >&2
    exit 1
fi
if [[ -e "$test_dir/compensating_mismatch.tsv" ]]; then
    echo "count_seq_depth wrote output after a per-pair mate-count failure" >&2
    exit 1
fi

# Each rejection is checked both before an output exists and with a previous
# valid result. The counter must finish validation before replacing that result.
reject_fastq() {
    local label="$1" input="$2" output="$test_dir/rejected.tsv"
    for existing in absent present; do
        rm -f "$output"
        if [[ "$existing" == present ]]; then
            cp "$test_dir/depth.tsv" "$output"
        fi
        if bash workflow/scripts/count_seq_depth.sh "$output" synthetic \
            --r1 "$input" --r2 "$input" >"$test_dir/rejected.log" 2>&1; then
            echo "count_seq_depth accepted $label" >&2
            exit 1
        fi
        if [[ "$existing" == present ]]; then
            cmp "$test_dir/depth.tsv" "$output"
        elif [[ -e "$output" ]]; then
            echo "count_seq_depth wrote output after $label" >&2
            exit 1
        fi
    done
}

# A complete leading record keeps these cases from passing by dropping a bad
# suffix and reporting only the valid population.
for malformed in incomplete header separator quality; do
    {
        printf '@valid\nACGT\n+\nIIII\n'
        case "$malformed" in
            incomplete) printf '@bad\nACGT\n+\n' ;;
            header) printf 'bad\nACGT\n+\nIIII\n' ;;
            separator) printf '@bad\nACGT\n-\nIIII\n' ;;
            quality) printf '@bad\nACGT\n+\nIII\n' ;;
        esac
    } | gzip -c > "$test_dir/$malformed.fastq.gz"
    reject_fastq "$malformed record" "$test_dir/$malformed.fastq.gz"
    # Structural diagnostics identify both the library and input file.
    grep -Fq synthetic "$test_dir/rejected.log"
    grep -Fq "$malformed.fastq.gz" "$test_dir/rejected.log"
done

# Remove the gzip trailer, retaining decompressed records that would otherwise
# look complete. pipefail must still reject the decompression failure.
size=$(wc -c < "$test_dir/sample_R1.fastq.gz")
dd if="$test_dir/sample_R1.fastq.gz" of="$test_dir/truncated.fastq.gz" \
    bs=1 count=$((size - 8)) 2>/dev/null
reject_fastq "truncated gzip" "$test_dir/truncated.fastq.gz"

# Four-line FASTQ allows a repeated header after + and a final quality line
# without a terminating newline. CRLF line endings do not contribute bases.
printf '@read1\r\nACGT\r\n+read1\r\nIIII' | gzip -c > "$test_dir/boundary.fastq.gz"
bash workflow/scripts/count_seq_depth.sh "$test_dir/boundary.tsv" synthetic \
    --r1 "$test_dir/boundary.fastq.gz" --r2 "$test_dir/boundary.fastq.gz" >/dev/null
awk -F '\t' 'NR == 2 && $3 == 1 && $4 == 2 && $5 == 8 { ok=1 } END { exit !ok }' \
    "$test_dir/boundary.tsv"

printf '' | gzip -c > "$test_dir/empty.fastq.gz"
bash workflow/scripts/count_seq_depth.sh "$test_dir/empty.tsv" synthetic \
    --r1 "$test_dir/empty.fastq.gz" --r2 "$test_dir/empty.fastq.gz" >/dev/null
awk -F '\t' 'NR == 2 && $2 == 1 && $3 == 0 && $4 == 0 && $5 == 0 { ok=1 } END { exit !ok }' \
    "$test_dir/empty.tsv"

echo "count_seq_depth paired-record and FASTQ-structure regressions passed"
