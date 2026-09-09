#!/usr/bin/env python3
import sys

for line in sys.stdin:
    if line.startswith("##"):
        sys.stdout.write(line)
    elif line.startswith("#CHROM"):
        # inject header definitions if not present
        sys.stdout.write('##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">\n')
        sys.stdout.write('##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Difference in length between REF and ALT alleles">\n')
        sys.stdout.write(line)
    else:
        fields = line.rstrip("\n").split("\t")
        chrom, pos, _id, ref, alt, qual, filt, info = fields[:8]
        pos = int(pos)

        # parse END from INFO (required for SVLEN)
        end = None
        for kv in info.split(";"):
            if kv.startswith("END="):
                try:
                    end = int(kv.split("=")[1])
                except ValueError:
                    pass

        # add SVTYPE=INS if missing
        if "SVTYPE=" not in info:
            if info == "." or info == "":
                info = "SVTYPE=DEL"
            else:
                info += ";SVTYPE=DEL"

        # add SVLEN if END is available
        if end is not None and "SVLEN=" not in info:
            svlen = end - pos
            info += ";SVLEN=%d" % svlen

        fields[7] = info
        sys.stdout.write("\t".join(fields) + "\n")
