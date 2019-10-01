#!/usr/bin/env python
import pandas as pd
import numpy as np
from sys import argv
import gzip

script, impvcf, imphard, impmgf, samples = argv


gtkey = {"0|0": "0", "0|1":"1", "1|0":"1", "1|1":"2"}
snpinfo = ["snp", "ref", "alt"]
with open(samples, 'w') as sampfile:
	with gzip.open(imphard, 'wt') as hardfile:
		with gzip.open(impmgf, 'wt') as mgffile:
			with gzip.open(impvcf, 'rt') as vcf:
				headers = []
				for line in vcf:
					if line.startswith("#"):
						print(line)
						legend = line.strip().split("\t")[:9]
						iids = line.strip().split("\t")[9:]
						if iids != []:
							#iidline = ",".join(iids)
							#sampfile.write(iidline+"\n")
							[sampfile.write(str(xx) + "\n") for xx in iids]
					if not line.startswith("#"):
						info = line.strip().split("\t")[:9]
						gts = line.strip().split("\t")[9:]
						mgfinfo = info[2:5]
						dose = [xx.split(":")[1] for xx in gts]
						hard = [gtkey[xx.split(":")[0]] for xx in gts]
						doseentry = mgfinfo+dose
						hardentry = mgfinfo+hard
						doseout = ",".join(doseentry)
						hardout = ",".join(hardentry)
						mgffile.write(doseout+"\n")
						hardfile.write(hardout+"\n")
