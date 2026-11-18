#!/usr/bin/env python3

import argparse
import csv
import sys
from pathlib import Path

# Increase CSV field size limit for large annotations
max_int = sys.maxsize
while True:
    try:
        csv.field_size_limit(max_int)
        break
    except OverflowError:
        max_int = int(max_int / 10)


def main():
    parser = argparse.ArgumentParser(
        description="Merge multiple TSV reports by column name, handling column order differences gracefully."
    )
    parser.add_argument("-o", "--out", required=True, help="Output merged TSV file")
    parser.add_argument("inputs", nargs="+", help="Input TSV files to merge")
    args = parser.parse_args()

    input_paths = [Path(p) for p in args.inputs if Path(p).exists() and Path(p).stat().st_size > 0]
    if not input_paths:
        sys.exit("ERROR: No non-empty input TSV files found to merge.")

    # 1. Discover all unique column names, preserving order from the first file
    fieldnames = []
    seen = set()

    for p in input_paths:
        with open(p, "r", encoding="utf-8", errors="replace", newline="") as f:
            reader = csv.reader(f, delimiter="\t")
            try:
                header = next(reader)
            except StopIteration:
                continue
            for col in header:
                col = col.strip()
                if col and col not in seen:
                    seen.add(col)
                    fieldnames.append(col)

    if not fieldnames:
        sys.exit("ERROR: No valid header columns found in any input file.")

    # 2. Write merged output
    total_rows = 0
    with open(args.out, "w", encoding="utf-8", newline="") as fout:
        writer = csv.DictWriter(fout, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore", lineterminator="\n")
        writer.writeheader()

        for p in input_paths:
            with open(p, "r", encoding="utf-8", errors="replace", newline="") as fin:
                reader = csv.DictReader(fin, delimiter="\t")
                for row in reader:
                    clean_row = {k.strip(): v for k, v in row.items() if k is not None}
                    writer.writerow(clean_row)
                    total_rows += 1

    print(f"Successfully merged {len(input_paths)} files ({total_rows} total rows) into: {args.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
