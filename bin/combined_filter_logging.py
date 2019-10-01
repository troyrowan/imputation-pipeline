#!/usr/bin/env python
import glob
from sys import argv
script, runname = argv

filenames = [xx.strip("ref_alt/").strip(".log") for xx in glob.glob("ref_alt/" + runname + "/*log")]

for xx in filenames:
    snpfile = "snp_filtered/" + xx + ".log"
    indfile = "individual_filtered/" + xx + ".log"
    hwefile = "hwe_filtered/" + xx + ".log"
    with open(snpfile, "r") as snp:
        for line in snp:
            if "--bfile" in line:
                assay = line.split("/")[2].strip()
                print("Assay:\t" + assay + "\n")
                print("SNP_Filtering")
            if "variants loaded from .bim file" in line:
                start_snps = line.split()[0]
                print("StartingSNPs:\t", start_snps)
            if "removed due to missing genotype data" in line:
                rem_snps = line.split()[0]
                end_snps = int(start_snps) - int(rem_snps)
                proportion = int(rem_snps)/int(start_snps)
                percent = str(round(int(rem_snps)/int(start_snps)*100, 3))
                print("RemovedSNPs:\t" + rem_snps)
                print("EndSNPs:\t" + str(end_snps))
                print("PercentSNPRemoved:\t" + percent + "%")
                #print(proportion)
                if proportion > 0.1:
                    print("!!HOLY SHIT WE DROPPED TOO MUCH!!\n")
                else:
                    print("SNP Filtering okay\n")
    with open(indfile, "r") as ind:
        for line in ind:
            if "loaded from .fam" in line:
                start_inds = line.split()[0]
                print("Individual_Filtering")
                print("StartingInd:\t" + start_inds)
            if "removed due to missing genotype data" in line:
                rem_inds = line.split()[0]
                end_inds = int(start_inds)-int(rem_inds)
                proportion = int(rem_inds)/int(start_inds)
                percent = str(round(int(rem_inds)/int(start_inds)*100, 3))
                print("RemovedInd:\t" + rem_inds)
                print("EndInds:\t" + str(end_inds))
                print("PercentIndRemoved:\t" + percent + "%")
                if proportion > 0.1:
                    print("!!HOLY SHIT WE DROPPED TOO MUCH!!\n")
                else:
                    print("Individual Filtering okay\n")
    # with open(hwefile, "r") as hwe:
    #     for line in hwe:
    #         if "variants loaded from .bim file" in line:
    #             start_snps = line.split()[0]
    #             print("HWE_Filtering")
    #             print("StartingSNPs:\t", start_snps)
    #         if "removed due to Hardy-Weinberg exact test." in line:
    #             rem_snps = line.split()[1]
    #             end_snps = int(start_snps) - int(rem_snps)
    #             proportion = int(rem_snps)/int(start_snps)
    #             percent = str(round(int(rem_snps)/int(start_snps)*100, 3))
    #             print("RemovedSNPs:\t" + rem_snps)
    #             print("EndSNPs:\t" + str(end_snps))
    #             print("PercentSNPRemoved:\t" + percent + "%")
    #             #print(proportion)
    #             if proportion > 0.1:
    #                 print("!!HOLY SHIT WE DROPPED TOO MUCH!!\n\n")
    #             else:
    #                 print("HWE Filtering okay\n\n")
