#!/usr/bin/env bash

set -euo pipefail

BASE_DIR="/pl/active/dashnowlab/work/ealiyev/SVTR_Analysis/DEVESON_STR/atarva"
OUTFILE="${BASE_DIR}/shaikh_atarva.final.new_zscore.tsv"

rm -f "$OUTFILE"

first_file=true

for sample_dir in "${BASE_DIR}"/*; do
    [ -d "$sample_dir" ] || continue

    sample=$(basename "$sample_dir")

    # Find LongTR final TSV inside the sample folder
    tsv=$(find "$sample_dir" -maxdepth 1 -type f -name "*.final.tsv" | head -n 1)

    if [ -z "$tsv" ]; then
        echo "WARNING: No *.atarva.final.tsv found in: $sample_dir" >&2
        continue
    fi

    if [ "$first_file" = true ]; then
        awk -v sample="$sample" 'BEGIN{OFS="\t"} NR==1 {print "sample", $0; next} {print sample, $0}' "$tsv" > "$OUTFILE"
        first_file=false
    else
        awk -v sample="$sample" 'BEGIN{OFS="\t"} NR>1 {print sample, $0}' "$tsv" >> "$OUTFILE"
    fi
done

echo "Merged TSV written to:"
echo "$OUTFILE"