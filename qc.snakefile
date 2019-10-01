rule filter_target:
	input:
		#targ = expand("merged_files/{run_name}/{run_name}.bed", run_name = config["run_name"])
		hwe = expand("hwe_filtered/{run_name}/{sample}.log", sample = config["sample"], run_name = config["run_name"])
	shell:
		"rm .snakemake/*_tracking/*"

#This map dictionary should be able to remain the same, and we can add new maps for whichever new assays become available in future datasets

#map dictionary can be found in the config files,
#It needs to be edited as new assays/maps become available
def mapdicter(WC): #Chooses which map file to use for associated ped file based on SNP number that is the first chunk of file name
	t = WC.sample
	num_sites = t.split('.')[0] #the file naming involves the number of sites in the genotype query, and so indicates which map to use.
	return config["mapdict"][num_sites]
def refdicter(WC): # To find ref allele
	t = WC.sample
	num_sites = t.split('.')[0] #the file naming involves the number of sites in the genotype query, and so indicates which map to use.
	return config["refdict"][num_sites]
def refalleledicter(WC):
	t = WC.sample
	num_sites = t.split('.')[0] #the file naming involves the number of sites in the genotype query, and so indicates which map to use.
	return config["refalleledict"][num_sites]

def sampdicter(wildcards):
	return samp_dict[wildcards.sample]

#This step sets the A1 allele to the actual reference allele.
#Also removes chromosome 0 variants and converts to a PLINK binary format
rule ref_alt:
	input:
		# bed = expand("{prefix}{{sample}}.bed", prefix = config["gt_path"]),
		# bim = expand("{prefix}{{sample}}.bim", prefix = config["gt_path"]),
		# fam = expand("{prefix}{{sample}}.fam", prefix = config["gt_path"]),
		ped = expand("{prefix}{{sample}}.ped", prefix = config["gt_path"]),
		map = mapdicter,
		ref = refdicter,
		allele = refalleledicter
	params:
		#inprefix = expand("{prefix}{{sample}}", prefix = config["gt_path"]),
		oprefix = "ref_alt/{run_name}/{sample}"
	benchmark:
		"benchmarks/ref_alt/{run_name}/{sample}.txt"
	log:
		"logs/ref_alt/{run_name}/{sample}.log"
	output:
		bed = "ref_alt/{run_name}/{sample}.bed",
		bim = "ref_alt/{run_name}/{sample}.bim",
		fam = "ref_alt/{run_name}/{sample}.fam"
	shell:
		#"(plink --bed {input.bed} --bim {input.bim} --fam {input.fam} --chr-set 32 --exclude maps/map_issues/190521_UMD_badsnps.txt --update-alleles {input.ref} --not-chr 0 --a1-allele {input.allele} --make-bed --real-ref-alleles --out {params.oprefix}) > {log}"
		"(plink --ped {input.ped} --map {input.map} --chr-set 32 --update-alleles {input.ref} --exclude maps/map_issues/190521_UMD_badsnps.txt --not-chr 0 --a1-allele {input.allele} --make-bed --real-ref-alleles --out {params.oprefix}) > {log}"

#Creative, easier way to identify duplicates (individuals genotyped multiple times on the same assay)
#This just merges files (from previous step) with themselves, and takes care of duplicates
rule no_duplicates:
	input:
		bed = "ref_alt/{run_name}/{sample}.bed",
		bim = "ref_alt/{run_name}/{sample}.bim",
		fam = "ref_alt/{run_name}/{sample}.fam"
	params:
		inprefix = "ref_alt/{run_name}/{sample}",
		oprefix = "no_duplicates/{run_name}/{sample}"
	benchmark:
		"benchmarks/no_duplicates/{run_name}/{sample}.txt"
	log:
		"logs/no_duplicates/{run_name}/{sample}.log"
	output:
		bed = "no_duplicates/{run_name}/{sample}.bed",
		bim = "no_duplicates/{run_name}/{sample}.bim",
		fam = "no_duplicates/{run_name}/{sample}.fam"
	shell:
 		"(plink --bfile {params.inprefix} --bmerge {params.inprefix} --chr-set 32 --real-ref-alleles --make-bed --out {params.oprefix}) > {log}"
#Calculates per-variant call rate
rule variant_stats:
	input:
		bed = "no_duplicates/{run_name}/{sample}.bed",
		bim = "no_duplicates/{run_name}/{sample}.bim",
		fam = "no_duplicates/{run_name}/{sample}.fam"
	threads : 4
	priority:100
	params:
		inprefix = "no_duplicates/{run_name}/{sample}",
		oprefix = "snp_stats/{run_name}/{sample}"
	benchmark:
		"benchmarks/variant_stats/{run_name}/{sample}.txt"
	log:
		"logs/variant_stats/{run_name}/{sample}.log"
	output:
		frq = "snp_stats/{run_name}/{sample}.frq", #This is input for plotting script
		log = "snp_stats/{run_name}/{sample}.log",
		#png = "snp_stats/figures/{run_name}/{sample}.snp_call_rate.png"
	shell:
		"(plink --bfile {params.inprefix} --chr-set 32 --memory 500 --real-ref-alleles --nonfounders --freq --out {params.oprefix}) > {log}"
		# "(plink --bfile {params.inprefix} --chr-set 32 --memory 500 --real-ref-alleles --nonfounders --freq --out {params.oprefix}; python bin/snp_call_rate_visualization.py {output.frq} {output.png}) > {log}"
#Python Script for Visualization (snp_call_rate_visualization.py) reads in .frq file generated by PLINK's --freq function
#Creates a histogram of the missing genotype rate for each SNP in the dataset.  This call rate is calculated for each SNP as (NCHROBS)/max(NCHROBS)
#Writes figure to a .png file
#Histogram uses 100 bins, and a y-axis max of 5000 SNPs, which look good for existing assays, but can be easily adjusted.

#After call rate calculations, this performs the by-variant filtering
rule filter_variants:
	input:
		bed = "no_duplicates/{run_name}/{sample}.bed",
		bim = "no_duplicates/{run_name}/{sample}.bim",
		fam = "no_duplicates/{run_name}/{sample}.fam",
		stats = "snp_stats/{run_name}/{sample}.frq",
		#png = "snp_stats/figures/{run_name}/{sample}.snp_call_rate.png"
	threads : 4
	priority:99
	params:
		inprefix = "no_duplicates/{run_name}/{sample}",
		oprefix="snp_filtered/{run_name}/{sample}",
		logprefix="filter_logs/{run_name}/{sample}"
	benchmark:
		"benchmarks/filter_variants/{run_name}/{sample}.txt"
	output:
		bed=temp("snp_filtered/{run_name}/{sample}.bed"),
		bim=temp("snp_filtered/{run_name}/{sample}.bim"),
		fam=temp("snp_filtered/{run_name}/{sample}.fam"),
		log="snp_filtered/{run_name}/{sample}.log"
	shell:
		"plink --bfile {params.inprefix} --threads 4 --chr-set 32 --memory 500 --real-ref-alleles --not-chr 0 --geno .1 --make-bed --out {params.oprefix}"



#PLINK statistic (missing) -- Gathers statistics on individual call rates (reports proportion of missing genotypes for each individual)
#Output: output file suffixes will be .imiss, .lmiss
rule individual_stats: #This step is performed on variant-filtered files
	input:
		bed = "snp_filtered/{run_name}/{sample}.bed",
		bim = "snp_filtered/{run_name}/{sample}.bim",
		fam = "snp_filtered/{run_name}/{sample}.fam",
	params:
		inprefix="snp_filtered/{run_name}/{sample}",
		oprefix="individual_stats/{run_name}/{sample}"
	benchmark:
		"benchmarks/individual_stats/{run_name}/{sample}.txt"
	output:
		imiss="individual_stats/{run_name}/{sample}.imiss",
		lmiss="individual_stats/{run_name}/{sample}.lmiss",
		log ="individual_stats/{run_name}/{sample}.log",
		#png="individual_stats/figures/{run_name}/{sample}.individual_call_rate.png"
	shell:
		"plink --bfile {params.inprefix} --chr-set 32 --memory 500 --real-ref-alleles --missing --out {params.oprefix}"
		# "plink --bfile {params.inprefix} --chr-set 32 --memory 500 --real-ref-alleles --missing --out {params.oprefix}; python bin/individual_call_rate_visualization.py {output.imiss} {output.png}"
#Python Script for Visualization (individual_call_rate_visualization.py) reads in .imiss file generated by PLINK's --missing function
#Creates a histogram of the missing genotyping rate for each individual in the dataset (F_MISS)
#Writes figure to a .png file
#Histogram uses 50 bins, and a y-axis max of 50 animals, which look good for existing assays, but can be easily adjusted.


#PLINK filter (mind) -- Filters individuals based on specified call rate.  Individuals missing more than given proportion of their genotypes will be filtered out of the dataset.
rule filter_individuals:
	input:
		bed = "snp_filtered/{run_name}/{sample}.bed",
		bim = "snp_filtered/{run_name}/{sample}.bim",
		fam = "snp_filtered/{run_name}/{sample}.fam",
		imiss="individual_stats/{run_name}/{sample}.imiss",
		#png="individual_stats/figures/{run_name}/{sample}.individual_call_rate.png"
	params:
		inprefix="snp_filtered/{run_name}/{sample}",
		oprefix="individual_filtered/{run_name}/{sample}"
	benchmark:
		"benchmarks/filter_individuals/{run_name}/{sample}.txt"
	output:
		bed="individual_filtered/{run_name}/{sample}.bed",
		bim="individual_filtered/{run_name}/{sample}.bim",
		fam="individual_filtered/{run_name}/{sample}.fam",
		log="individual_filtered/{run_name}/{sample}.log"

	shell: #Filter is set at 0.1 missing as a threshold for being dropped from the dataset
		"plink --bfile {params.inprefix} --chr-set 32 --memory 500 --real-ref-alleles --mind .1  --make-bed --out {params.oprefix}"#; python bin/individual_filtered_log_parsing.py {output.log} {params.csv}"

#Gathers statistics on individual Hardy Weinberg Equilibrium (reports HWE P value at each locus)
#Then visualizes with Python plotting script
rule hwe_stats:
	input:
		bed="individual_filtered/{run_name}/{sample}.bed",
		bim="individual_filtered/{run_name}/{sample}.bim",
		fam="individual_filtered/{run_name}/{sample}.fam"
	params:
		inprefix="individual_filtered/{run_name}/{sample}",
		oprefix="hwe_stats/{run_name}/{sample}"
	benchmark:
		"benchmarks/hwe_stats/{run_name}/{sample}.txt"
	output:
		hwe="hwe_stats/{run_name}/{sample}.hwe",
		log="hwe_stats/{run_name}/{sample}.log",
		#png="hwe_stats/figures/{run_name}/{sample}.hwe_pvalues.png"
	shell:
		"plink --bfile {params.inprefix} --chr-set 32 --memory 500 --real-ref-alleles --nonfounders --hardy --out {params.oprefix}"
		# "plink --bfile {params.inprefix} --chr-set 32 --memory 500 --real-ref-alleles --nonfounders --hardy --out {params.oprefix}; python bin/hwe_visualization.py {output.hwe} {output.png}"

#PLINK filter (hwe) -- Filters loci based on their HWE P values. Variants with HWE P values below a specified value will be removed from the dataset
# We've set this value at 1e-50, but this is context dependent on dataset. Composite reference data this will drop A TON of variants, even though they're perfectly fine
rule filter_hwe_variants:
	input:
		bed="individual_filtered/{run_name}/{sample}.bed",
		bim="individual_filtered/{run_name}/{sample}.bim",
		fam="individual_filtered/{run_name}/{sample}.fam",
		stats="hwe_stats/{run_name}/{sample}.hwe",
	params:
		inprefix="individual_filtered/{run_name}/{sample}",
		oprefix="hwe_filtered/{run_name}/{sample}"
	benchmark:
		"benchmarks/filter_hwe_variants/{run_name}/{sample}.txt"
	output:
		bed="hwe_filtered/{run_name}/{sample}.bed",
		bim="hwe_filtered/{run_name}/{sample}.bim",
		fam="hwe_filtered/{run_name}/{sample}.fam",
		log="hwe_filtered/{run_name}/{sample}.log"
	shell:
		"plink --bfile {params.inprefix} --chr-set 32 --memory 500 --real-ref-alleles --nonfounders --exclude subsetted_test/190903_UMD_GELtest/badsnps.txt --hwe 1e-50 --make-bed --out {params.oprefix}"

#We aren't doing any sort of sex-checking in our QC at this point. Will do this in the future as we go to impute X and Y chromosomes

#PLINK function (impute-sex) looks at sex provided by ped file, and at X chromosome heterozygosity (and y chromosome variants if provided), and determines whether an animal is male, female, or unknown. If sex from ped file is unknown, this will impute the sex if possible, and write that into the new bed file that it produces.
#Python Scripts: (missexed_animals_filter.py) reads the .sexcheck file produced by impute sex function.  If an individual's sex is changed from known M/F to the opposite sex, it's ID will be written to the {sample}.missexed_animals.txt file, and removed in subsequent step.  (missexed_animals_filter_logging.py) will count lines of filter output, and write to the assay's csv log file.
##Input File Location: hwe_filtered/
##Output File Location: sex_impute/
##Output: output file suffixes will be (.bed, .bim, .fam, .log, .sexcheck, .missexed_animals.txt)
# rule impute_sex:
# 	input:
# 		bed="hwe_filtered/{run_name}/{sample}.bed",
# 		bim="hwe_filtered/{run_name}/{sample}.bim",
# 		fam="hwe_filtered/{run_name}/{sample}.fam"
# 	params:
# 		logprefix="filter_logs/{run_name}/{sample}",
# 		inprefix="hwe_filtered/{run_name}/{sample}",
# 		oprefix="sex_impute/{run_name}/{sample}"
# 	benchmark:
# 		"benchmarks/impute_sex/{run_name}/{sample}.txt"
# 	output:
# 		bed=temp("sex_impute/{run_name}/{sample}.bed"),
# 		bim=temp("sex_impute/{run_name}/{sample}.bim"),
# 		fam=temp("sex_impute/{run_name}/{sample}.fam"),
# 		log="sex_impute/{run_name}/{sample}.log",
# 		sexcheck="sex_impute/{run_name}/{sample}.sexcheck",
# 		txt="sex_impute/{run_name}/{sample}.missexed_animals.txt"
# 	shell:
# 		"plink --bfile {params.inprefix} --chr-set 32 --memory 500 --real-ref-alleles --impute-sex ycount --make-bed --out {params.oprefix}; python bin/missexed_animals_filter.py {output.sexcheck} {output.txt}"# python bin/missexed_animals_filter_logging.py {output.txt} {output.log} {params.csv}" #{params.logprefix}.csv"
#
# #PLINK function (remove) takes "sex_impute/{run_name}/{sample}.missexed_animals.txt" and produces new bed file without listed IDs, removing missexed animals
# #Input File Locations: sex_impute/
# #Output File Locations: correct_sex/
# #Output File Types: (.bed, .bim, .fam, .log, .hh, .nosex)
#
# rule remove_missexed_animals:
# 	input:
# 		txt="sex_impute/{run_name}/{sample}.missexed_animals.txt",
# 		bed="sex_impute/{run_name}/{sample}.bed",
# 		bim="sex_impute/{run_name}/{sample}.bim",
# 		fam="sex_impute/{run_name}/{sample}.fam"
# 	params:
# 		inprefix="sex_impute/{run_name}/{sample}",
# 		oprefix="correct_sex/{run_name}/{sample}"
# 	benchmark:
# 		"benchmarks/remove_missexed_animals/{run_name}/{sample}.txt"
# 	output:
# 		bed="correct_sex/{run_name}/{sample}.bed",
# 		bim="correct_sex/{run_name}/{sample}.bim",
# 		fam="correct_sex/{run_name}/{sample}.fam",
# 		#hh=temp("correct_sex/{run_name}/{sample}.hh"),
# 		#nosex=temp("correct_sex/{run_name}/{sample}.nosex"),
# 		log="correct_sex/{run_name}/{sample}.log",
# 	shell:
# 		"plink --bfile {params.inprefix} --chr-set 32 --memory 500 --real-ref-alleles --remove {input.txt} --make-bed --out {params.oprefix}"

#This script performs filter logging for all assays after they're generated
#Parses log files and prints to sdout, point that to a filtering report.
rule filter_logging:
	input:
		snp = expand("snp_filtered/{{run_name}}/{sample}.log", sample = config["sample"]),
		ind = expand("individual_filtered/{{run_name}}/{sample}.log", sample = config["sample"]),
		hwe = expand("hwe_filtered/{{run_name}}/{sample}.log", sample = config["sample"])
		#sex = "sex_impute/{run_name}/{sample}.log",
		#mis = "sex_impute/{run_name}/{sample}.missexed_animals.txt",
		#cor = "correct_sex/{run_name}/{sample}.bed"
	params:
		prefix = "{run_name}"
	benchmark:
		"benchmarks/filter_logging/{run_name}/{run_name}.txt"
	log:
		"logs/filter_logging/{run_name}/{run_name}.log"
	output:
		log = "filter_logs/{run_name}/{run_name}_filtering_report.txt"
	shell:
		"python bin/combined_filter_logging.py {params.prefix} > {output.log}"
#Filtering report will be in the following format:
#For each assay:
	# Assay:	47843.180808.3095.B
	#
	# SNP_Filtering
	# StartingSNPs:	 45880
	# RemovedSNPs:	365
	# EndSNPs:	45515
	# PercentSNPRemoved:	0.796%
	# SNP Filtering okay
	#
	# Individual_Filtering
	# StartingInd:	3089
	# RemovedInd:	49
	# EndInds:	3040
	# PercentIndRemoved:	1.586%
	# Individual Filtering okay
	#
	# HWE_Filtering
	# StartingSNPs:	 45515
	# RemovedSNPs:	82
	# EndSNPs:	45433
	# PercentSNPRemoved:	0.18%
	# HWE Filtering okay
		#Will throw flags when too many SNPs or individuals are removed

#Prior to phasing/imputation, we merge all of the different assays into a single files
#The resulting files have ~850K variants for n individuals, with un-genotyped markers as missing
#Because we use reference-based phasing, these missing markers aren't filled in by Eagle (saving a step of removing Eagle-inferred genotypes prior to imputation)
rule merge_assays:
	input:
		expand("hwe_filtered/{{run_name}}/{sample}.bed", sample = config["sample"]),
		expand("hwe_filtered/{{run_name}}/{sample}.bim", sample = config["sample"]),
		expand("hwe_filtered/{{run_name}}/{sample}.fam", sample = config["sample"]),
		expand("hwe_filtered/{{run_name}}/{sample}.log", sample = config["sample"]),
		# expand("individual_filtered/{{run_name}}/{sample}.bed", sample = config["sample"]),
		# expand("individual_filtered/{{run_name}}/{sample}.bim", sample = config["sample"]),
		# expand("individual_filtered/{{run_name}}/{sample}.fam", sample = config["sample"]),
		# expand("individual_filtered/{{run_name}}/{sample}.log", sample = config["sample"]),
		"filter_logs/{run_name}/{run_name}_filtering_report.txt"
	params:
		oprefix="merged_files/{run_name}/{run_name}",
		pfiles = "hwe_filtered/{run_name}"
		#pfiles = "individual_filtered/{run_name}"
	output:
		mergefilelist= "individual_filtered/{run_name}/{run_name}.allfiles.txt",
		bim = "merged_files/{run_name}/{run_name}.bim",
		fam = "merged_files/{run_name}/{run_name}.fam",
		log = "merged_files/{run_name}/{run_name}.log",
		bed = "merged_files/{run_name}/{run_name}.bed"
	shell: #This list maker creates a list of all assays being imputed. This is given to PLINK's "--merge-list" command
		"python bin/file_list_maker.py {params.pfiles} {output.mergefilelist}; plink --merge-list {output.mergefilelist} --chr-set 32 --merge-equal-pos --real-ref-alleles --make-bed --out {params.oprefix}"
