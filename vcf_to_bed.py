#!/usr/bin/env python3

import sys
import re


def parse_info(info_str):
    d = {}
    if info_str == "." or info_str.strip() == "":
        return d

    for item in info_str.split(";"):
        if not item:
            continue
        if "=" in item:
            k, v = item.split("=", 1)
            d[k] = v
        else:
            d[item] = "True"
    return d


def parse_coordinates_from_id(vid):
    """
    Supports Medaka-style IDs:
        chr1_57367043_57367118

    Also supports:
        chr1:57367043-57367118
    """

    m = re.match(r"^(.+?)_(\d+)_(\d+)$", vid)
    if m:
        return m.group(1), int(m.group(2)), int(m.group(3))

    m = re.match(r"^(.+?):(\d+)-(\d+)$", vid)
    if m:
        return m.group(1), int(m.group(2)), int(m.group(3))

    return None, None, None


def normalize_id(chrom, pos, original_id, info):
    """
    Make one stable ID used by both vcf_to_bed.py and vcf_to_tsv.py.

    Caller behavior:
      LongTR:
        VCF ID exists, e.g. 1-1435798-1435818-GGCGCGGAGC

      Medaka:
        VCF ID exists, e.g. chr1_57367043_57367118

      STRdust:
        VCF ID is usually '.', INFO has END only
        Use chr_POS_END

      ATaRVa:
        VCF ID is '.', real locus ID is in INFO/ID
        Use INFO/ID, e.g. HMNR7_VWA1
    """

    if original_id and original_id != ".":
        return original_id

    if "ID" in info and info["ID"] not in ["", "."]:
        return info["ID"]

    end_for_id = info.get("END", pos)
    return f"{chrom}_{pos}_{end_for_id}"


def infer_bed_interval(chrom, pos, vid, ref, info):
    """
    Returns:
        bed_chrom, bed_start, bed_end

    Priority:
      1. INFO START + END: LongTR and ATaRVa
      2. INFO END only: STRdust
      3. ID coordinates: Medaka
      4. POS + REF length fallback
    """

    # LongTR / ATaRVa
    if "START" in info and "END" in info:
        start = int(info["START"])
        end = int(info["END"])

        # Important:
        # ATaRVa START appears already 0-based in your examples.
        # LongTR START appears 1-based in your examples.
        #
        # We decide based on POS:
        #   if START == POS - 1, START is already BED-like
        #   if START == POS, convert to POS - 1
        pos_i = int(pos)

        if start == pos_i - 1:
            bed_start = start
        else:
            bed_start = start - 1

        return chrom, bed_start, end

    # STRdust
    if "END" in info:
        pos_i = int(pos)
        end = int(info["END"])
        return chrom, pos_i - 1, end

    # Medaka
    id_chrom, id_start, id_end = parse_coordinates_from_id(vid)
    if id_start is not None and id_end is not None:
        return chrom, id_start, id_end

    # Fallback
    pos_i = int(pos)
    bed_start = pos_i - 1
    bed_end = bed_start + len(ref)
    return chrom, bed_start, bed_end


def main():
    if len(sys.argv) != 3:
        sys.exit("Usage: python vcf_to_bed.py input.vcf output.bed")

    input_vcf = sys.argv[1]
    output_bed = sys.argv[2]

    n_total = 0
    n_written = 0
    n_skipped = 0

    with open(input_vcf) as fin, open(output_bed, "w") as fout:
        fout.write("#CHROM\tSTART\tEND\tSVTYPE\tID\n")

        for line in fin:
            if line.startswith("#"):
                continue

            n_total += 1

            fields = line.rstrip("\n").split("\t")
            if len(fields) < 8:
                n_skipped += 1
                continue

            chrom = fields[0]
            pos = fields[1]
            original_id = fields[2]
            ref = fields[3]
            info = parse_info(fields[7])

            vid = normalize_id(chrom, pos, original_id, info)

            try:
                bed_chrom, bed_start, bed_end = infer_bed_interval(
                    chrom=chrom,
                    pos=pos,
                    vid=vid if vid != original_id else original_id,
                    ref=ref,
                    info=info
                )
            except Exception as e:
                n_skipped += 1
                print(
                    f"WARNING: skipping {chrom}:{pos} ID={vid}: {e}",
                    file=sys.stderr
                )
                continue

            if bed_start < 0 or bed_end <= bed_start:
                n_skipped += 1
                print(
                    f"WARNING: invalid interval skipped: {chrom}:{pos} ID={vid} "
                    f"START={bed_start} END={bed_end}",
                    file=sys.stderr
                )
                continue

            svtype = info.get("SVTYPE", "DEL")

            fout.write(f"{bed_chrom}\t{bed_start}\t{bed_end}\t{svtype}\t{vid}\n")
            n_written += 1

    print(f"Total VCF records:   {n_total}", file=sys.stderr)
    print(f"Written BED records: {n_written}", file=sys.stderr)
    print(f"Skipped records:     {n_skipped}", file=sys.stderr)
    print(f"Output BED:          {output_bed}", file=sys.stderr)


if __name__ == "__main__":
    main()