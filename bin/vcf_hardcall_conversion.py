#!/usr/bin/env python
import gzip
from sys import argv

script, vcf, out = argv

with open(out, 'w') as outvcf:
    with gzip.open(vcf, 'rt') as vcf:
        headers = []
        for line in vcf:
            if line.startswith("#"):
                headers.append(line)
                outvcf.write(line)
            else:
                info = line.strip().split("\t")[:9]
                gts = line.strip().split("\t")[9:]
                hard = [xx.split(":")[0] for xx in gts]
                entry = info + hard
                outvcf.write("\t". join(entry)+ "\n")
