#!/usr/bin/env python3
import sys

for line in sys.stdin:
    if line.startswith("##"):
        sys.stdout.write(line)
    elif line.startswith("#CHROM"):
        # inject header for TRID
        sys.stdout.write('##INFO=<ID=TRID,Number=1,Type=String,Description="Tandem repeat locus ID: chrom_POS_END_STR">\n')
        sys.stdout.write(line)
    else:
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 8:
            sys.stdout.write(line)
            continue

        chrom, pos, _id, ref, alt, qual, filt, info = fields[:8]
        pos_i = int(pos)

        # parse END from INFO
        end = None
        if info and info != ".":
            for kv in info.split(";"):
                if kv.startswith("END="):
                    try:
                        end = int(kv.split("=", 1)[1])
                    except ValueError:
                        pass

        # add TRID if END found
        if end is not None and "TRID=" not in info:
            trid = f"{chrom}_{pos}_{end}_STR"
            info = (info + ";" if info and info != "." else "") + f"TRID={trid}"

        fields[7] = info
        sys.stdout.write("\t".join(fields) + "\n")
