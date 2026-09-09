#!/usr/bin/env python3

import sys
import math


def parse_info(info_str):
    d = {}
    if info_str in ("", "."):
        return d
    for item in info_str.split(";"):
        if not item:
            continue
        if "=" in item:
            k, v = item.split("=", 1)
            d[k] = v
        else:
            d[item] = True
    return d


def info_to_string(info_dict):
    if not info_dict:
        return "."
    out = []
    for k, v in info_dict.items():
        if v is True:
            out.append(k)
        else:
            out.append(f"{k}={v}")
    return ";".join(out)


def parse_format_sample(fmt_str, sample_str):
    keys = fmt_str.split(":") if fmt_str not in ("", ".") else []
    vals = sample_str.split(":") if sample_str not in ("", ".") else []
    if len(vals) < len(keys):
        vals += [""] * (len(keys) - len(vals))
    return dict(zip(keys, vals))


def parse_gt(gt):
    if gt in ("", ".", "./.", ".|."):
        return [None, None], None

    sep = "|" if "|" in gt else "/"
    parts = gt.split(sep)

    alleles = []
    for p in parts[:2]:
        if p in ("", "."):
            alleles.append(None)
        else:
            try:
                alleles.append(int(p))
            except ValueError:
                alleles.append(None)

    while len(alleles) < 2:
        alleles.append(None)

    return alleles, sep


def safe_float(x):
    try:
        if x in ("", ".", None):
            return float("nan")
        return float(x)
    except Exception:
        return float("nan")


def fmt_num(x, digits=6):
    if x is None or not isinstance(x, (int, float)) or math.isnan(x):
        return "."
    x = float(x)
    if abs(x) > 0 and abs(x) < 10 ** (-digits):
        return f"{x:.{digits}e}"
    return f"{x:.{digits}f}".rstrip("0").rstrip(".")


def rotate_motif(s):
    return {s[i:] + s[:i] for i in range(len(s))} if s else set()


def reverse_complement(seq):
    comp = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    return seq.translate(comp)[::-1]


def motif_matches(query_motif, anno_motif):
    if not query_motif or not anno_motif:
        return 0

    query_motif = query_motif.upper()
    anno_motif = anno_motif.upper()

    if query_motif == anno_motif:
        return 2

    q_rots = rotate_motif(query_motif)
    a_rots = rotate_motif(anno_motif)

    if q_rots & a_rots:
        return 1

    rc_rots = rotate_motif(reverse_complement(query_motif))
    if rc_rots & a_rots:
        return 1

    return 0


def split_pipe_field(info_dict, key):
    raw = info_dict.get(key, ".")
    if raw in ("", ".", None):
        return []
    return raw.split("|")


def get_pipe_value(info_dict, key, idx):
    arr = split_pipe_field(info_dict, key)
    if idx is None or idx >= len(arr):
        return None
    val = arr[idx]
    if val in ("", ".", None):
        return None
    return val


def select_motif_index(info_dict, vcf_motif):
    motifs = split_pipe_field(info_dict, "stranno_motif")
    if not motifs:
        return None

    statuses = split_pipe_field(info_dict, "stranno_motif_status")

    # Medaka/STRdust may not have MOTIF/PERIOD from original caller.
    # If stranno gives only one motif, use it.
    if not vcf_motif and len(motifs) == 1:
        return 0

    best_idx = None
    best_score = -1

    for i, motif in enumerate(motifs):
        status = statuses[i] if i < len(statuses) else ""
        match_score = motif_matches(vcf_motif, motif)

        if match_score == 2:
            score = 100
        elif match_score == 1:
            score = 50
        else:
            score = 0

        if status == "E":
            score += 10
        elif status == "R":
            score += 5

        if score > best_score:
            best_score = score
            best_idx = i

    # If no motif match but stranno has motifs, keep the first motif instead of losing the record.
    if best_score <= 0:
        return 0

    return best_idx


def get_motif_specific_counts(info_dict, motif_idx):
    token = get_pipe_value(info_dict, "stranno_motif_counts", motif_idx)
    if token is None:
        return None

    counts = []
    for p in token.split(","):
        counts.append(safe_float(p))
    return counts


def get_motif_specific_stat(info_dict, key, motif_idx):
    token = get_pipe_value(info_dict, key, motif_idx)
    return safe_float(token)


def get_allele_seq(ref, alt_list, allele_index):
    if allele_index is None:
        return None
    if allele_index == 0:
        return ref

    alt_pos = allele_index - 1
    if 0 <= alt_pos < len(alt_list):
        return alt_list[alt_pos]

    return None


def get_repeat_count_from_sequence(ref, alt_list, allele_index, period):
    seq = get_allele_seq(ref, alt_list, allele_index)
    if seq is None or math.isnan(period) or period <= 0:
        return float("nan")
    return len(seq) / period


def zscore(v, mu, sd):
    v = safe_float(v)
    mu = safe_float(mu)
    sd = safe_float(sd)

    if any(math.isnan(x) for x in (v, mu, sd)):
        return float("nan")
    if sd == 0:
        return float("nan")

    return (v - mu) / sd


def norm_sf(z):
    z = safe_float(z)
    if math.isnan(z):
        return float("nan")
    return 0.5 * math.erfc(z / math.sqrt(2.0))


def p_adj_bh(values):
    vals = list(values)
    out = vals.copy()

    finite = [
        (i, v)
        for i, v in enumerate(vals)
        if isinstance(v, (int, float)) and math.isfinite(v)
    ]

    if not finite:
        return out

    finite_sorted = sorted(finite, key=lambda x: x[1])
    m = len(finite_sorted)

    adj = [float("nan")] * len(vals)
    prev = 1.0

    for rank_rev, (orig_idx, p) in enumerate(reversed(finite_sorted), start=1):
        rank = m - rank_rev + 1
        val = min(1.0, p * m / rank)
        prev = min(prev, val)
        adj[orig_idx] = prev

    for i, v in enumerate(vals):
        if isinstance(v, (int, float)) and math.isfinite(v):
            out[i] = adj[i]

    return out


def get_sample_col_idx(header_fields, sample_name):
    if sample_name in header_fields:
        return header_fields.index(sample_name)

    # If exact sample name does not match but VCF has exactly one sample column,
    # use it. This helps ATaRVa where header sample may be CD00014 instead of CD00014.URfixed.
    if len(header_fields) == 10:
        actual_sample = header_fields[9]
        print(
            f"WARNING: Sample '{sample_name}' not found. Using the only VCF sample column: '{actual_sample}'",
            file=sys.stderr
        )
        return 9

    sys.exit(
        f"ERROR: Sample '{sample_name}' not found in VCF header.\n"
        f"Available sample columns: {header_fields[9:]}"
    )


def get_vcf_motif_and_period(info):
    motif = info.get("MOTIF", "")

    if not motif:
        motifs = split_pipe_field(info, "stranno_motif")
        if len(motifs) == 1:
            motif = motifs[0]

    period = safe_float(info.get("PERIOD", "."))

    if math.isnan(period) or period <= 0:
        if motif:
            period = float(len(motif))

    return motif, period


def main():
    if len(sys.argv) != 4:
        sys.exit("Usage: python z_score.py <input.vcf> <sample_name> <output.vcf>")

    input_vcf = sys.argv[1]
    sample_name = sys.argv[2]
    output_vcf = sys.argv[3]

    new_headers = [
        '##INFO=<ID=TR_MATCHED_MOTIF_INDEX,Number=1,Type=Integer,Description="Selected index in stranno_motif used for count/stat lookup">',
        '##INFO=<ID=TR_MATCHED_MOTIF,Number=1,Type=String,Description="Selected motif from stranno_motif used for count/stat lookup">',
        '##INFO=<ID=TR_H1_REPCN,Number=1,Type=Float,Description="Allele 1 repeat count from matched stranno_motif_counts or sequence-length fallback">',
        '##INFO=<ID=TR_H2_REPCN,Number=1,Type=Float,Description="Allele 2 repeat count from matched stranno_motif_counts or sequence-length fallback">',

        '##INFO=<ID=TR_TENK10K_H1_Z,Number=1,Type=Float,Description="Allele 1 z-score using stranno_tenk10k median/stdev">',
        '##INFO=<ID=TR_TENK10K_H2_Z,Number=1,Type=Float,Description="Allele 2 z-score using stranno_tenk10k median/stdev">',
        '##INFO=<ID=TR_TENK10K_H1_P,Number=1,Type=Float,Description="Allele 1 upper-tail p-value">',
        '##INFO=<ID=TR_TENK10K_H2_P,Number=1,Type=Float,Description="Allele 2 upper-tail p-value">',
        '##INFO=<ID=TR_TENK10K_H1_PADJ,Number=1,Type=Float,Description="Allele 1 BH-adjusted p-value">',
        '##INFO=<ID=TR_TENK10K_H2_PADJ,Number=1,Type=Float,Description="Allele 2 BH-adjusted p-value">',

        '##INFO=<ID=TR_AOU1027_H1_Z,Number=1,Type=Float,Description="Allele 1 z-score using stranno_aou1027 median/stdev">',
        '##INFO=<ID=TR_AOU1027_H2_Z,Number=1,Type=Float,Description="Allele 2 z-score using stranno_aou1027 median/stdev">',
        '##INFO=<ID=TR_AOU1027_H1_P,Number=1,Type=Float,Description="Allele 1 upper-tail p-value">',
        '##INFO=<ID=TR_AOU1027_H2_P,Number=1,Type=Float,Description="Allele 2 upper-tail p-value">',
        '##INFO=<ID=TR_AOU1027_H1_PADJ,Number=1,Type=Float,Description="Allele 1 BH-adjusted p-value">',
        '##INFO=<ID=TR_AOU1027_H2_PADJ,Number=1,Type=Float,Description="Allele 2 BH-adjusted p-value">',

        '##INFO=<ID=TR_HPRC256_H1_Z,Number=1,Type=Float,Description="Allele 1 z-score using stranno_hprc256 median/stdev">',
        '##INFO=<ID=TR_HPRC256_H2_Z,Number=1,Type=Float,Description="Allele 2 z-score using stranno_hprc256 median/stdev">',
        '##INFO=<ID=TR_HPRC256_H1_P,Number=1,Type=Float,Description="Allele 1 upper-tail p-value">',
        '##INFO=<ID=TR_HPRC256_H2_P,Number=1,Type=Float,Description="Allele 2 upper-tail p-value">',
        '##INFO=<ID=TR_HPRC256_H1_PADJ,Number=1,Type=Float,Description="Allele 1 BH-adjusted p-value">',
        '##INFO=<ID=TR_HPRC256_H2_PADJ,Number=1,Type=Float,Description="Allele 2 BH-adjusted p-value">',
    ]

    meta_lines = []
    header_line = None
    records = []
    sample_col_idx = None

    with open(input_vcf) as fin:
        for line in fin:
            if line.startswith("##"):
                meta_lines.append(line)
                continue

            if line.startswith("#CHROM"):
                header_fields = line.rstrip("\n").split("\t")
                sample_col_idx = get_sample_col_idx(header_fields, sample_name)
                header_line = line
                continue

            fields = line.rstrip("\n").split("\t")
            if len(fields) < 10:
                continue

            chrom, pos, vid, ref, alt, qual, flt, info_str, fmt = fields[:9]
            sample_str = fields[sample_col_idx]

            fmt_dict = parse_format_sample(fmt, sample_str)
            gt = fmt_dict.get("GT", "")

            alleles, sep = parse_gt(gt)

            # Skip fully missing calls, but allow partial calls like STRdust 1|.
            if alleles[0] is None and alleles[1] is None:
                continue

            # Keep 0/0 records too if they have useful stranno statistics.
            # If you want to remove homozygous reference calls later, filter final TSV.
            info = parse_info(info_str)

            vcf_motif, period = get_vcf_motif_and_period(info)
            alt_list = [] if alt in ("", ".") else alt.split(",")

            motif_idx = None
            matched_motif = None

            a1_rep = get_repeat_count_from_sequence(ref, alt_list, alleles[0], period)
            a2_rep = get_repeat_count_from_sequence(ref, alt_list, alleles[1], period)

            z_values = {
                "tenk_h1": float("nan"),
                "tenk_h2": float("nan"),
                "aou_h1": float("nan"),
                "aou_h2": float("nan"),
                "hprc_h1": float("nan"),
                "hprc_h2": float("nan"),
            }

            has_stranno = "stranno_motif" in info and "stranno_motif_counts" in info

            if has_stranno:
                motif_idx = select_motif_index(info, vcf_motif)

                if motif_idx is not None:
                    matched_motif = get_pipe_value(info, "stranno_motif", motif_idx)
                    counts = get_motif_specific_counts(info, motif_idx)

                    if counts is not None:
                        if alleles[0] is not None and alleles[0] < len(counts):
                            a1_rep = counts[alleles[0]]
                        if alleles[1] is not None and alleles[1] < len(counts):
                            a2_rep = counts[alleles[1]]

                    stats = {
                        "tenk": (
                            get_motif_specific_stat(info, "stranno_tenk10k_median", motif_idx),
                            get_motif_specific_stat(info, "stranno_tenk10k_stdev", motif_idx),
                        ),
                        "aou": (
                            get_motif_specific_stat(info, "stranno_aou1027_median", motif_idx),
                            get_motif_specific_stat(info, "stranno_aou1027_stdev", motif_idx),
                        ),
                        "hprc": (
                            get_motif_specific_stat(info, "stranno_hprc256_median", motif_idx),
                            get_motif_specific_stat(info, "stranno_hprc256_stdev", motif_idx),
                        ),
                    }

                    for pop, (med, sd) in stats.items():
                        z_values[f"{pop}_h1"] = zscore(a1_rep, med, sd)
                        z_values[f"{pop}_h2"] = zscore(a2_rep, med, sd)

            rec = {
                "fields": fields,
                "info": info,
                "motif_idx": motif_idx,
                "matched_motif": matched_motif,
                "a1_rep": a1_rep,
                "a2_rep": a2_rep,
                **z_values,
            }

            records.append(rec)

    if header_line is None:
        sys.exit("ERROR: VCF header not found")

    for rec in records:
        for pop in ["tenk", "aou", "hprc"]:
            rec[f"p_{pop}_h1"] = norm_sf(rec[f"{pop}_h1"])
            rec[f"p_{pop}_h2"] = norm_sf(rec[f"{pop}_h2"])

    for pop in ["tenk", "aou", "hprc"]:
        padj_h1 = p_adj_bh([r[f"p_{pop}_h1"] for r in records])
        padj_h2 = p_adj_bh([r[f"p_{pop}_h2"] for r in records])

        for i, rec in enumerate(records):
            rec[f"padj_{pop}_h1"] = padj_h1[i]
            rec[f"padj_{pop}_h2"] = padj_h2[i]

    with open(output_vcf, "w") as fout:
        existing_header_ids = set()

        for line in meta_lines:
            fout.write(line)
            if line.startswith("##INFO=<ID="):
                try:
                    existing_header_ids.add(line.split("ID=", 1)[1].split(",", 1)[0])
                except Exception:
                    pass

        for h in new_headers:
            try:
                hid = h.split("ID=", 1)[1].split(",", 1)[0]
            except Exception:
                hid = None

            if hid and hid not in existing_header_ids:
                fout.write(h + "\n")

        fout.write(header_line)

        for rec in records:
            info = rec["info"]

            if rec["motif_idx"] is not None:
                info["TR_MATCHED_MOTIF_INDEX"] = str(rec["motif_idx"])

            if rec["matched_motif"] is not None:
                info["TR_MATCHED_MOTIF"] = rec["matched_motif"]

            info["TR_H1_REPCN"] = fmt_num(rec["a1_rep"])
            info["TR_H2_REPCN"] = fmt_num(rec["a2_rep"])

            info["TR_TENK10K_H1_Z"] = fmt_num(rec["tenk_h1"])
            info["TR_TENK10K_H2_Z"] = fmt_num(rec["tenk_h2"])
            info["TR_TENK10K_H1_P"] = fmt_num(rec["p_tenk_h1"])
            info["TR_TENK10K_H2_P"] = fmt_num(rec["p_tenk_h2"])
            info["TR_TENK10K_H1_PADJ"] = fmt_num(rec["padj_tenk_h1"])
            info["TR_TENK10K_H2_PADJ"] = fmt_num(rec["padj_tenk_h2"])

            info["TR_AOU1027_H1_Z"] = fmt_num(rec["aou_h1"])
            info["TR_AOU1027_H2_Z"] = fmt_num(rec["aou_h2"])
            info["TR_AOU1027_H1_P"] = fmt_num(rec["p_aou_h1"])
            info["TR_AOU1027_H2_P"] = fmt_num(rec["p_aou_h2"])
            info["TR_AOU1027_H1_PADJ"] = fmt_num(rec["padj_aou_h1"])
            info["TR_AOU1027_H2_PADJ"] = fmt_num(rec["padj_aou_h2"])

            info["TR_HPRC256_H1_Z"] = fmt_num(rec["hprc_h1"])
            info["TR_HPRC256_H2_Z"] = fmt_num(rec["hprc_h2"])
            info["TR_HPRC256_H1_P"] = fmt_num(rec["p_hprc_h1"])
            info["TR_HPRC256_H2_P"] = fmt_num(rec["p_hprc_h2"])
            info["TR_HPRC256_H1_PADJ"] = fmt_num(rec["padj_hprc_h1"])
            info["TR_HPRC256_H2_PADJ"] = fmt_num(rec["padj_hprc_h2"])

            rec["fields"][7] = info_to_string(info)
            fout.write("\t".join(rec["fields"]) + "\n")

    print(f"Input records retained: {len(records)}", file=sys.stderr)
    print(f"Done. Wrote annotated VCF: {output_vcf}", file=sys.stderr)


if __name__ == "__main__":
    main()