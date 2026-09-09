#!/usr/bin/env python3

import argparse
from collections import OrderedDict


def parse_info_field(info_str):
    info = {}
    if info_str == "." or info_str.strip() == "":
        return info

    for item in info_str.split(";"):
        if not item:
            continue
        if "=" in item:
            k, v = item.split("=", 1)
            info[k] = v
        else:
            info[item] = "True"
    return info


def normalize_id(chrom, pos, original_id, info):
    """
    Same logic as vcf_to_bed.py.

    LongTR:
      keep VCF ID

    Medaka:
      keep VCF ID, e.g. chr1_57367043_57367118

    ATaRVa:
      VCF ID is '.', use INFO/ID

    STRdust:
      VCF ID is '.', no INFO/ID, use chr_POS_END
    """

    if original_id and original_id != ".":
        return original_id

    if "ID" in info and info["ID"] not in ["", "."]:
        return info["ID"]

    end_for_id = info.get("END", pos)
    return f"{chrom}_{pos}_{end_for_id}"


def parse_meta_id(line, prefix):
    if not line.startswith(prefix):
        return None

    try:
        inner = line[line.index("<") + 1: line.rindex(">")]
        parts = inner.split(",")
        for p in parts:
            if p.startswith("ID="):
                return p.replace("ID=", "", 1)
    except Exception:
        return None

    return None


def discover_fields(vcf_path):
    info_fields = []
    format_fields = []
    info_seen = set()
    format_seen = set()
    sample_names = []

    extra_info = OrderedDict()
    extra_format = OrderedDict()

    with open(vcf_path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.rstrip("\n")

            if line.startswith("##INFO="):
                fid = parse_meta_id(line, "##INFO=")
                if fid and fid not in info_seen:
                    info_seen.add(fid)
                    info_fields.append(fid)

            elif line.startswith("##FORMAT="):
                fid = parse_meta_id(line, "##FORMAT=")
                if fid and fid not in format_seen:
                    format_seen.add(fid)
                    format_fields.append(fid)

            elif line.startswith("#CHROM"):
                cols = line.split("\t")
                if len(cols) > 9:
                    sample_names = cols[9:]
                continue

            elif line.startswith("#"):
                continue

            else:
                cols = line.split("\t")
                if len(cols) < 8:
                    continue

                info_dict = parse_info_field(cols[7])
                for k in info_dict:
                    if k not in info_seen and k not in extra_info:
                        extra_info[k] = None

                if len(cols) >= 9 and cols[8] != ".":
                    fmt_keys = cols[8].split(":")
                    for k in fmt_keys:
                        if k not in format_seen and k not in extra_format:
                            extra_format[k] = None

    info_fields.extend(list(extra_info.keys()))
    format_fields.extend(list(extra_format.keys()))

    return info_fields, format_fields, sample_names


def convert_vcf_to_tsv(vcf_path, out_path):
    info_fields, format_fields, sample_names = discover_fields(vcf_path)

    fixed_cols = ["#CHROM", "POS", "ID", "REF", "ALT_1", "ALT_2"]

    format_output_cols = []
    if sample_names:
        for sample in sample_names:
            for fmt in format_fields:
                format_output_cols.append(f"{sample}__{fmt}")
    else:
        for fmt in format_fields:
            format_output_cols.append(fmt)

    header = fixed_cols + info_fields + format_output_cols

    n_total = 0
    n_written = 0

    with open(vcf_path, "r", encoding="utf-8") as fin, \
         open(out_path, "w", encoding="utf-8") as fout:

        fout.write("\t".join(header) + "\n")

        for line in fin:
            line = line.rstrip("\n")

            if not line or line.startswith("#"):
                continue

            n_total += 1

            cols = line.split("\t")
            if len(cols) < 8:
                continue

            chrom = cols[0]
            pos = cols[1]
            original_id = cols[2]
            ref = cols[3]
            alt = cols[4]

            info_dict = parse_info_field(cols[7])
            vid = normalize_id(chrom, pos, original_id, info_dict)

            alts = alt.split(",") if alt != "." else []
            alt_1 = alts[0] if len(alts) >= 1 else "."
            alt_2 = alts[1] if len(alts) >= 2 else "."

            row = [chrom, pos, vid, ref, alt_1, alt_2]

            for field in info_fields:
                row.append(info_dict.get(field, "."))

            if sample_names and len(cols) > 9:
                fmt_keys = cols[8].split(":") if cols[8] != "." else []
                sample_data = cols[9:]

                if len(sample_data) < len(sample_names):
                    sample_data.extend(["."] * (len(sample_names) - len(sample_data)))

                for sample_value in sample_data[:len(sample_names)]:
                    values = sample_value.split(":") if sample_value != "." else []
                    fmt_map = {k: v for k, v in zip(fmt_keys, values)}

                    for field in format_fields:
                        row.append(fmt_map.get(field, "."))

            elif not sample_names:
                for _ in format_fields:
                    row.append(".")

            fout.write("\t".join(row) + "\n")
            n_written += 1

    print(f"Total VCF records:   {n_total}")
    print(f"Written TSV records: {n_written}")
    print(f"Output TSV:          {out_path}")


def main():
    parser = argparse.ArgumentParser(
        description="Convert VCF to TSV with normalized IDs, INFO fields, and FORMAT fields."
    )
    parser.add_argument("-i", "--input", required=True, help="Input VCF file")
    parser.add_argument("-o", "--output", required=True, help="Output TSV file")
    args = parser.parse_args()

    convert_vcf_to_tsv(args.input, args.output)


if __name__ == "__main__":
    main()