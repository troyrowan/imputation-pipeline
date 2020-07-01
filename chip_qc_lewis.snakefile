import os
configfile: "source_functions/config/find_dups.config.yaml"
# Make log directories if they don't exist
for x in expand("log/slurm_out/{rules}", rules = config['rules']):
    os.makedirs(x, exist_ok = True)
for x in expand("log/psrecord/joint_genotyping/{rules}", rules = config['rules']):
    os.makedirs(x, exist_ok = True)

#Generates a single merged file with all filtering
rule filter_target:
	input:

		hwe = expand("imputation_runs/imputation_runs/{run_name}/merged_files/{run_name}.bed",
		run_name = config["run_name"])
	shell:
		"rm .snakemake/*_tracking/*"

#Map dictionary is copied into each new config file allowing run-by-run changes to map file
#This dict will need to be changed as new assays are added.
#If assay is absent, pipeline will fail in first step.
def mapchooser(WC): #Chooses which map file to use for associated ped file based on SNP number that is the first chunk of file name
	t = WC.sample
	num_sites = t.split('.')[0] #the file naming involves the number of sites in the genotype query, and so indicates which map to use.
	return config["mapdict"][num_sites]
def refchooser(WC): # To find ref allele
	t = WC.sample
	num_sites = t.split('.')[0] #the file naming involves the number of sites in the genotype query, and so indicates which map to use.
	return config["refdict"][num_sites]
def refallelechooser(WC):
	t = WC.sample
	num_sites = t.split('.')[0] #the file naming involves the number of sites in the genotype query, and so indicates which map to use.
	return config["refalleledict"][num_sites]

def sampchooser(wildcards):
	return samp_dict[wildcards.sample]

#This step sets the A1 allele to the actual reference allele.
#Also removes chromosome 0 variants and converts to a PLINK binary format
rule ref_alt:
	input:
		ped = expand("{prefix}{{sample}}.ped", prefix = config["gt_path"]),
		map = mapchooser,
		ref = refchooser,
		allele = refallelechooser
	params:
		threads = config["plink_threads"],
		mem = config["plink_mem"],
		oprefix = "imputation_runs/{run_name}/ref_alt/{sample}",
		psrecord = "log/{run_name}/psrecord/ref_alt/ref_alt.{sample}.log"
	output:
		bed = temp("imputation_runs/{run_name}/ref_alt/{sample}.bed"),
		bim = "imputation_runs/{run_name}/ref_alt/{sample}.bim",
		fam = "imputation_runs/{run_name}/ref_alt/{sample}.fam"
	shell:
		"""
		psrecord "plink --ped {input.ped} --map {input.map} --cow --update-alleles {input.ref} --not-chr 0 --a1-allele {input.allele} --threads {params.threads} --mem {params.mem} --make-bed --real-ref-alleles --out {params.oprefix}" --log {params.psrecord} --include-children --interval 5
		"""

#Creative, easier way to identify duplicates (individuals genotyped multiple times on the same assay)
#This just merges files (from previous step) with themselves, and takes care of duplicates
rule no_duplicates:
	input:
		bed = "imputation_runs/{run_name}/ref_alt/{sample}.bed",
		bim = "imputation_runs/{run_name}/ref_alt/{sample}.bim",
		fam = "imputation_runs/{run_name}/ref_alt/{sample}.fam"
	params:
		inprefix = "imputation_runs/{run_name}/ref_alt/{sample}",
		oprefix = "imputation_runs/{run_name}/no_duplicates/{sample}",
		threads = config["plink_threads"],
		mem = config["plink_mem"],
		psrecord = "log/{run_name}/psrecord/no_duplicates/no_duplicates.{sample}.log"
	output:
		bed = temp("imputation_runs/{run_name}/no_duplicates/{sample}.bed"),
		bim = "imputation_runs/{run_name}/no_duplicates/{sample}.bim",
		fam = "imputation_runs/{run_name}/no_duplicates/{sample}.fam"
	shell:
 		"""
		psrecord "plink --bfile {params.inprefix} --bmerge {params.inprefix} --cow --threads {params.threads} --mem {params.mem} --real-ref-alleles --make-bed --out {params.oprefix}" --log {params.psrecord} --include-children --interval 5
		"""
#Calculates per-variant call rate
rule variant_stats:
	input:
		bed = "imputation_runs/{run_name}/no_duplicates/{sample}.bed",
		bim = "imputation_runs/{run_name}/no_duplicates/{sample}.bim",
		fam = "imputation_runs/{run_name}/no_duplicates/{sample}.fam"
	params:
		inprefix = "imputation_runs/{run_name}/no_duplicates/{sample}",
		oprefix = "imputation_runs/{run_name}/snp_stats/{sample}",
		threads = config["plink_threads"],
		mem = config["plink_mem"],
		psrecord = "log/{run_name}/psrecord/variant_stats/variant_stats.{sample}.log"
	output:
		frq = "imputation_runs/{run_name}/snp_stats/{sample}.frq", #This is input for plotting script
		log = "imputation_runs/{run_name}/snp_stats/{sample}.log",
		#png = "snp_stats/figures/{sample}.snp_call_rate.png"
	shell:
		""""
		psrecord "plink --bfile {params.inprefix} --cow --threads {params.threads} --mem {params.mem} --real-ref-alleles --nonfounders --freq --out {params.oprefix}" --log {params.psrecord} --include-children --interval 5
		""""
		# "(plink --bfile {params.inprefix} --cow --memory 500 --real-ref-alleles --nonfounders --freq --out {params.oprefix}; python bin/snp_call_rate_visualization.py {output.frq} {output.png}) > {log}"
#Python Script for Visualization (snp_call_rate_visualization.py) reads in .frq file generated by PLINK's --freq function
#Creates a histogram of the missing genotype rate for each SNP in the dataset.  This call rate is calculated for each SNP as (NCHROBS)/max(NCHROBS)
#Writes figure to a .png file
#Histogram uses 100 bins, and a y-axis max of 5000 SNPs, which look good for existing assays, but can be easily adjusted.
#UPDATE: Python code for plotting not working at the moment
#After call rate calculations, this performs the by-variant filtering
rule filter_variants:
	input:
		bed = "imputation_runs/{run_name}/no_duplicates/{sample}.bed",
		bim = "imputation_runs/{run_name}/no_duplicates/{sample}.bim",
		fam = "imputation_runs/{run_name}/no_duplicates/{sample}.fam",
		stats = "imputation_runs/{run_name}/snp_stats/{sample}.frq",
		#png = "snp_stats/figures/{sample}.snp_call_rate.png"
	threads : 4
	priority:99
	params:
		inprefix = "imputation_runs/{run_name}/no_duplicates/{sample}",
		oprefix="imputation_runs/{run_name}/snp_filtered/{sample}",
		logprefix="imputation_runs/{run_name}/filter_logs/{sample}",
		threads = config["plink_threads"],
		mem = config["plink_mem"],
		filter = config["snp_callrate_filter"],
		psrecord = "log/{run_name}/psrecord/filter_variants/filter_variants.{sample}.log"
	output:
		bed=temp("imputation_runs/{run_name}/snp_filtered/{sample}.bed"),
		bim="imputation_runs/{run_name}/snp_filtered/{sample}.bim",
		fam="imputation_runs/{run_name}/snp_filtered/{sample}.fam",
		log="imputation_runs/{run_name}/snp_filtered/{sample}.log"
	shell:
		"""
		psrecord "plink --bfile {params.inprefix} --cow --threads {params.threads} --mem {params.mem} --real-ref-alleles --nonfounders --not-chr 0 --geno {params.filter} --make-bed --out {params.oprefix}" --log {params.psrecord} --include-children --interval 5
		"""



#PLINK statistic (missing) -- Gathers statistics on individual call rates (reports proportion of missing genotypes for each individual)
#Output: output file suffixes will be .imiss, .lmiss
rule individual_stats: #This step is performed on variant-filtered files
	input:
		bed = "imputation_runs/{run_name}/snp_filtered/{sample}.bed",
		bim = "imputation_runs/{run_name}/snp_filtered/{sample}.bim",
		fam = "imputation_runs/{run_name}/snp_filtered/{sample}.fam",
	params:
		inprefix="imputation_runs/{run_name}/snp_filtered/{sample}",
		oprefix="imputation_runs/{run_name}/individual_stats/{sample}",
		threads = config["plink_threads"],
		mem = config["plink_mem"],
		psrecord = "log/{run_name}/psrecord/individual_stats/individual_stats.{sample}.log"
	output:
		imiss="imputation_runs/{run_name}/individual_stats/{sample}.imiss",
		lmiss="imputation_runs/{run_name}/individual_stats/{sample}.lmiss",
		log ="imputation_runs/{run_name}/individual_stats/{sample}.log",
		#png="individual_stats/figures/{sample}.individual_call_rate.png"
	shell:
		"""
		psrecord "plink --bfile {params.inprefix} --cow --threads {params.threads} --mem {params.mem} --real-ref-alleles --missing --out {params.oprefix}" --log {params.psrecord} --include-children --interval 5
		"""



#PLINK filter (mind) -- Filters individuals based on specified call rate.  Individuals missing more than given proportion of their genotypes will be filtered out of the dataset.
rule filter_individuals:
	input:
		bed = "imputation_runs/{run_name}/snp_filtered/{sample}.bed",
		bim = "imputation_runs/{run_name}/snp_filtered/{sample}.bim",
		fam = "imputation_runs/{run_name}/snp_filtered/{sample}.fam",
		imiss="imputation_runs/{run_name}/individual_stats/{sample}.imiss",
		#png="individual_stats/figures/{sample}.individual_call_rate.png"
	params:
		inprefix="imputation_runs/{run_name}/snp_filtered/{sample}",
		oprefix="imputation_runs/{run_name}/individual_filtered/{sample}",
		threads = config["plink_threads"],
		mem = config["plink_mem"],
		filter = config["ind_callrate_filter"]
		psrecord = "log/{run_name}/psrecord/filter_individuals/filter_individuals.{sample}.log"
	output:
		bed=temp("imputation_runs/{run_name}/individual_filtered/{sample}.bed"),
		bim="imputation_runs/{run_name}/individual_filtered/{sample}.bim",
		fam="imputation_runs/{run_name}/individual_filtered/{sample}.fam",
		log="imputation_runs/{run_name}/individual_filtered/{sample}.log"

	shell: #Filter is set at 0.1 missing as a threshold for being dropped from the dataset
		"""
		psrecord "plink --bfile {params.inprefix} --cow --threads {params.threads} --mem {params.mem} --real-ref-alleles --mind {params.filter}  --make-bed --out {params.oprefix}" --log {params.psrecord} --include-children --interval 5
		"""
		#; python bin/individual_filtered_log_parsing.py {output.log} {params.csv}"

#Gathers statistics on individual Hardy Weinberg Equilibrium (reports HWE P value at each locus)
#Then visualizes with Python plotting script
rule hwe_stats:
	input:
		bed="imputation_runs/{run_name}/individual_filtered/{sample}.bed",
		bim="imputation_runs/{run_name}/individual_filtered/{sample}.bim",
		fam="imputation_runs/{run_name}/individual_filtered/{sample}.fam"
	params:
		inprefix="imputation_runs/{run_name}/individual_filtered/{sample}",
		oprefix="imputation_runs/{run_name}/hwe_stats/{sample}",
		threads = config["plink_threads"],
		mem = config["plink_mem"],
		psrecord = "log/{run_name}/psrecord/hwe_stats/hwe_stats.{sample}.log"
	output:
		hwe="imputation_runs/{run_name}/hwe_stats/{sample}.hwe",
		log="imputation_runs/{run_name}/hwe_stats/{sample}.log",
		#png="hwe_stats/figures/{sample}.hwe_pvalues.png"
	shell:
		"""
		psrecord "plink --bfile {params.inprefix} --cow --threads {params.threads} --mem {params.mem} --real-ref-alleles --nonfounders --hardy --out {params.oprefix}" --log {params.psrecord} --include-children --interval 5
		"""
		# "plink --bfile {params.inprefix} --cow --memory 500 --real-ref-alleles --nonfounders --hardy --out {params.oprefix}; python bin/hwe_visualization.py {output.hwe} {output.png}"

#PLINK filter (hwe) -- Filters loci based on their HWE P values. Variants with HWE P values below a specified value will be removed from the dataset
# We've set this value at 1e-50, but this is context dependent on dataset. Composite reference data this will drop A TON of variants, even though they're perfectly fine
rule filter_hwe_variants:
	input:
		bed="imputation_runs/{run_name}/individual_filtered/{sample}.bed",
		bim="imputation_runs/{run_name}/individual_filtered/{sample}.bim",
		fam="imputation_runs/{run_name}/individual_filtered/{sample}.fam",
		stats="imputation_runs/{run_name}/hwe_stats/{sample}.hwe",
	params:
		inprefix="imputation_runs/{run_name}/individual_filtered/{sample}",
		oprefix="imputation_runs/{run_name}/hwe_filtered/{sample}",
		threads = config["plink_threads"],
		mem = config["plink_mem"],
		filter = config["hwe_filter"],
		psrecord = "log/{run_name}/psrecord/filter_hwe_variants/filter_hwe_variants.{sample}.log"
	output:
		bed="imputation_runs/{run_name}/hwe_filtered/{sample}.bed",
		bim="imputation_runs/{run_name}/hwe_filtered/{sample}.bim",
		fam="imputation_runs/{run_name}/hwe_filtered/{sample}.fam",
		log="imputation_runs/{run_name}/hwe_filtered/{sample}.log"
	shell:
		"""
		psrecord "plink --bfile {params.inprefix} --cow --threads {params.threads} --mem {params.mem} --real-ref-alleles --nonfounders --hwe 1e-50 --make-bed --out {params.oprefix}" --log {params.psrecord} --include-children --interval 5
		"""

#We aren't doing any sort of sex-checking in our QC at this point. Will do this in the future as we go to impute X and Y chromosomes

#PLINK function (impute-sex) looks at sex provided by ped file, and at X chromosome heterozygosity (and y chromosome variants if provided), and determines whether an animal is male, female, or unknown. If sex from ped file is unknown, this will impute the sex if possible, and write that into the new bed file that it produces.
#Python Scripts: (missexed_animals_filter.py) reads the .sexcheck file produced by impute sex function.  If an individual's sex is changed from known M/F to the opposite sex, it's ID will be written to the {sample}.missexed_animals.txt file, and removed in subsequent step.  (missexed_animals_filter_logging.py) will count lines of filter output, and write to the assay's csv log file.
##Input File Location: hwe_filtered/
##Output File Location: sex_impute/
##Output: output file suffixes will be (.bed, .bim, .fam, .log, .sexcheck, .missexed_animals.txt)
# rule impute_sex:
# 	input:
# 		bed="hwe_filtered/{sample}.bed",
# 		bim="hwe_filtered/{sample}.bim",
# 		fam="hwe_filtered/{sample}.fam"
# 	params:
# 		logprefix="filter_logs/{sample}",
# 		inprefix="hwe_filtered/{sample}",
# 		oprefix="sex_impute/{sample}"
# 	benchmark:
# 		"imputation_runs/{run_name}/benchmarks/impute_sex/{sample}.txt"
# 	output:
# 		bed=temp("sex_impute/{sample}.bed"),
# 		bim=temp("sex_impute/{sample}.bim"),
# 		fam=temp("sex_impute/{sample}.fam"),
# 		log="sex_impute/{sample}.log",
# 		sexcheck="sex_impute/{sample}.sexcheck",
# 		txt="sex_impute/{sample}.missexed_animals.txt"
# 	shell:
# 		"plink --bfile {params.inprefix} --cow --memory 500 --real-ref-alleles --impute-sex ycount --make-bed --out {params.oprefix}; python bin/missexed_animals_filter.py {output.sexcheck} {output.txt}"# python bin/missexed_animals_filter_logging.py {output.txt} {output.log} {params.csv}" #{params.logprefix}.csv"
#
# #PLINK function (remove) takes "sex_impute/{sample}.missexed_animals.txt" and produces new bed file without listed IDs, removing missexed animals
# #Input File Locations: sex_impute/
# #Output File Locations: correct_sex/
# #Output File Types: (.bed, .bim, .fam, .log, .hh, .nosex)
#
# rule remove_missexed_animals:
# 	input:
# 		txt="sex_impute/{sample}.missexed_animals.txt",
# 		bed="sex_impute/{sample}.bed",
# 		bim="sex_impute/{sample}.bim",
# 		fam="sex_impute/{sample}.fam"
# 	params:
# 		inprefix="sex_impute/{sample}",
# 		oprefix="correct_sex/{sample}"
# 	benchmark:
# 		"imputation_runs/{run_name}/benchmarks/remove_missexed_animals/{sample}.txt"
# 	output:
# 		bed="correct_sex/{sample}.bed",
# 		bim="correct_sex/{sample}.bim",
# 		fam="correct_sex/{sample}.fam",
# 		#hh=temp("correct_sex/{sample}.hh"),
# 		#nosex=temp("correct_sex/{sample}.nosex"),
# 		log="correct_sex/{sample}.log",
# 	shell:
# 		"plink --bfile {params.inprefix} --cow --memory 500 --real-ref-alleles --remove {input.txt} --make-bed --out {params.oprefix}"

#This script performs filter logging for all assays after they're generated
#Parses log files and prints to sdout, point that to a filtering report.
rule filter_logging:
	input:
		snp = expand("{{run_name}}/snp_filtered/{sample}.log", sample = config["sample"]),
		ind = expand("{{run_name}}/individual_filtered/{sample}.log", sample = config["sample"]),
		hwe = expand("{{run_name}}/hwe_filtered/{sample}.log", sample = config["sample"])
	params:
		prefix = "{run_name}",
		psrecord = "log/{run_name}/psrecord/filter_logging/filter_logging.{sample}.log"
	output:
		log = "imputation_runs/{run_name}/filter_logs/{run_name}_filtering_report.txt"
	shell:
		"""
		psrecord "python bin/combined_filter_logging.py {params.prefix} > {output.log}" --log {params.psrecord} --include-children --interval 5
		"""
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
		expand("{{run_name}}/hwe_filtered/{sample}.bed", sample = config["sample"]),
		expand("{{run_name}}/hwe_filtered/{sample}.bim", sample = config["sample"]),
		expand("{{run_name}}/hwe_filtered/{sample}.fam", sample = config["sample"]),
		expand("{{run_name}}/hwe_filtered/{sample}.log", sample = config["sample"])
	params:
		oprefix="imputation_runs/{run_name}/merged_files/{run_name}",
		pfiles = "imputation_runs/{run_name}/hwe_filtered",
		threads = config["plink_threads"],
		mem = config["plink_mem"],
		psrecord = "log/{run_name}/psrecord/merge_assays/merge_assays.{sample}.log"
	output:
		mergefilelist= "imputation_runs/{run_name}/merged_files/{run_name}.allfiles.txt",
		bim = "imputation_runs/{run_name}/merged_files/{run_name}.bim",
		fam = "imputation_runs/{run_name}/merged_files/{run_name}.fam",
		log = "imputation_runs/{run_name}/merged_files/{run_name}.log",
		bed = "imputation_runs/{run_name}/merged_files/{run_name}.bed"
	shell: #This list maker creates a list of all assays being imputed. This is given to PLINK's "--merge-list" command
		"""
		psrecord "python bin/file_list_maker.py {params.pfiles} {output.mergefilelist}; plink --merge-list {output.mergefilelist} --cow --merge-equal-pos --real-ref-alleles --make-bed --out {params.oprefix}" --log {params.psrecord} --include-children --interval 5
		"""
