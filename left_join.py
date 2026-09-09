#!/usr/bin/env python3

import pandas as pd
import argparse
import sys
from pathlib import Path


def read_table(path, sep):
    path = Path(path)

    if not path.exists():
        sys.exit(f"ERROR: File does not exist: {path}")

    if path.stat().st_size == 0:
        sys.exit(f"ERROR: File is empty: {path}")

    try:
        df = pd.read_csv(path, sep=sep, dtype=str, keep_default_na=False)
    except Exception as e:
        sys.exit(f"ERROR: Could not read file: {path}\n{e}")

    df.columns = [c.strip() for c in df.columns]

    return df


def check_keys(df, keys, label, file_path):
    missing = [k for k in keys if k not in df.columns]

    if missing:
        sys.exit(
            f"ERROR: Missing key column(s) in {label} file: {missing}\n"
            f"File: {file_path}\n"
            f"Available columns:\n{list(df.columns)}"
        )


def main():
    parser = argparse.ArgumentParser(
        description="Robust left join two CSV/TSV files by key columns."
    )

    parser.add_argument("left_file", help="Path to the left input file")
    parser.add_argument("right_file", help="Path to the right input file")
    parser.add_argument(
        "--left_key",
        required=True,
        help="Column name(s) in left file; comma-separated if multiple"
    )
    parser.add_argument(
        "--right_key",
        required=True,
        help="Column name(s) in right file; comma-separated if multiple"
    )
    parser.add_argument(
        "--out",
        default="joined_output.tsv",
        help="Output file name"
    )
    parser.add_argument(
        "--sep",
        default="\t",
        help="Delimiter for input/output files"
    )
    parser.add_argument(
        "--right_suffix",
        default="_AnnotSV",
        help="Suffix for overlapping right-side column names"
    )
    parser.add_argument(
        "--deduplicate_right",
        action="store_true",
        help="Keep only first row per right key"
    )

    args = parser.parse_args()

    left_keys = [k.strip() for k in args.left_key.split(",")]
    right_keys = [k.strip() for k in args.right_key.split(",")]

    if len(left_keys) != len(right_keys):
        sys.exit(
            "ERROR: Number of left and right keys must match.\n"
            f"Left keys: {left_keys}\n"
            f"Right keys: {right_keys}"
        )

    left_df = read_table(args.left_file, args.sep)
    right_df = read_table(args.right_file, args.sep)

    check_keys(left_df, left_keys, "left", args.left_file)
    check_keys(right_df, right_keys, "right", args.right_file)

    for k in left_keys:
        left_df[k] = left_df[k].astype(str).str.strip()

    for k in right_keys:
        right_df[k] = right_df[k].astype(str).str.strip()

    if args.deduplicate_right:
        before = len(right_df)
        right_df = right_df.drop_duplicates(subset=right_keys, keep="first")
        after = len(right_df)
        print(f"Right file deduplicated: {before} -> {after}", file=sys.stderr)

    duplicated_right = right_df.duplicated(subset=right_keys, keep=False).sum()

    if duplicated_right > 0 and not args.deduplicate_right:
        print(
            f"WARNING: right file has {duplicated_right} rows with duplicated keys. "
            f"This can increase output row count. Use --deduplicate_right if needed.",
            file=sys.stderr
        )

    merged_df = pd.merge(
        left_df,
        right_df,
        how="left",
        left_on=left_keys,
        right_on=right_keys,
        suffixes=("", args.right_suffix)
    )

    merged_df.to_csv(args.out, sep=args.sep, index=False)

    print(f"Left file rows:  {len(left_df)}", file=sys.stderr)
    print(f"Right file rows: {len(right_df)}", file=sys.stderr)
    print(f"Output rows:     {len(merged_df)}", file=sys.stderr)
    print(f"Left join completed: {args.out}", file=sys.stderr)


if __name__ == "__main__":
    main()