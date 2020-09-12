##include: "bigrefprep.snakefile"
#include: "190122_refcreation.snakefile"
include: "chip_qc_lewis.snakefile"
#include: "chip_ref_creation_lewis.snakefile"

#snakemake -s chip_imp_lewis.snakefile --jobs 1000 --rerun-incomplete --keep-going --latency-wait 30 --configfile chip_imp_lewis.config.yaml --cluster-config chip_imp_lewis.cluster.json --cluster "sbatch -p {cluster.p} -o {cluster.o} --account {cluster.account} -t {cluster.t} -c {cluster.c} --mem {cluster.mem} --account {cluster.account}" -np

import os
# Make log directories if they don't exist
for x in expand("log/{run_name}/slurm_out/{rule}", run_name = config['run_name'], rule = config['rules']):
	os.makedirs(x, exist_ok = True)
for x in expand("log/{run_name}/psrecord/{rule}", run_name = config['run_name'], rule = config['rules']):
	os.makedirs(x, exist_ok = True)

rule bigref_done: #Last file create is specified up here. Use expand to indicate how we want wild cards filled in.
	input:#Standard outputs for the pipeline are a dosage vcf file and a hardcall only vcf file. Have ability to make dosage input for GEMMA and other file types
		gen = expand("imputation_runs/{run_name}/imputed_genotypes/single_chrom/{run_name}.chr{chr}.reordered.vcf.gz",
		run_name = config["run_name"],
		chr = list(range(1,31)))
		# gen = expand("imputation_runs/{run_name}/assay_raw_vcf/{run_name}.chr{chr}.vcf.gz",
		# run_name = config["run_name"],
		# chr = list(range(1,31)))
		# mgf = expand("imputation_runs/{run_name}/imputed_genotypes/{run_name}.hardcall.vcf.gz",
		# run_name = config["run_name"])
		#chroms = expand("{run_name}/imputed_genotypes_single_chrom/{run_name}.chr{chr}.dose.mgf.gz", run_name = config["run_name"], chr = list(range(1,30)))
		#dose = expand("minimac_imputed/combined_imputed/mm4/{run_name}.chr{chr}.dose.vcf.gz", date = config["date"], run_name = config["run_name"], chr = list(range(1,30)))
		#mgf = expand("imputed_chromosomes/{run_name}.chr{chr}.dose.mgf.gz", run_name = config["run_name"], chr = "29")
	shell:
		"rm .snakemake/*_tracking/*"

#Input is the filtered merged file from PLINK filtering. Single set of bed/bim/fam files needed for this step
#Splits into chromosomes for parallelization
rule assay_chrsplit:
	input:
		bim = "imputation_runs/{run_name}/merged_files/{run_name}.bim",
		fam = "imputation_runs/{run_name}/merged_files/{run_name}.fam",
		log = "imputation_runs/{run_name}/merged_files/{run_name}.log",
		bed = "imputation_runs/{run_name}/merged_files/{run_name}.bed",
	params:
	 	inprefix = "imputation_runs/{run_name}/merged_files/{run_name}",
		vcf = "imputation_runs/{run_name}/assay_raw_vcf/{run_name}.chr{chr}.vcf",
		oprefix = "imputation_runs/{run_name}/assay_raw_vcf/{run_name}.chr{chr}",
		chr = "{chr}",
		threads=config["plink_threads"],
		mem=config["plink_mem"],
		psrecord = "log/{run_name}/psrecord/assay_chrsplit/assay_chrsplit.chr{chr}.log"
	output:
		vcf = temp("imputation_runs/{run_name}/assay_raw_vcf/{run_name}.chr{chr}.vcf.gz"), #These files are temporary, will be removed when no longer needed by new input rule
		tabix = temp("imputation_runs/{run_name}/assay_raw_vcf/{run_name}.chr{chr}.vcf.gz.tbi")
	shell: #PLINK is used to split merged files into 29 chromosomes and then recode as vcf, bgzipped post-recoding
		"""
		module load plink
		module load bcftools
		psrecord "plink --bfile {params.inprefix} --cow --real-ref-alleles --chr {params.chr} --nonfounders --threads {params.threads} --memory {params.mem} --recode vcf-iid --out {params.oprefix}; bgzip {params.vcf}; tabix {output.vcf}" --log {params.psrecord} --include-children --interval 5
		"""

rule bigref_phasing:
	input:
		refvcf = expand("reference_build/{ref_version}/vcf_per_assay/{ref}.chr{{chr}}.vcf.gz",
		ref = config["hdref"],
		ref_version = config["ref_version"]),
		vcf = "imputation_runs/{run_name}/assay_raw_vcf/{run_name}.chr{chr}.vcf.gz",
		tabix = "imputation_runs/{run_name}/assay_raw_vcf/{run_name}.chr{chr}.vcf.gz.tbi"
	params:
		imputemap = "bin/genetic_map_1cMperMb.txt",#This may be an issue if Eagle is installed somewhere else. One of the few hardcoded files in the pipeline
		out = "imputation_runs/{run_name}/eagle_merged/{run_name}.chr{chr}",
		chrom = "{chr}", #Chromosome wildcard for Eagle
		threads = config["eagle_threads"],
		psrecord = "log/{run_name}/psrecord/bigref_phasing/bigref_phasing.chr{chr}.log"
	output:
		vcf = "imputation_runs/{run_name}/eagle_merged/{run_name}.chr{chr}.vcf.gz", #Eagle outputs a bgzipped VCF file, we tabix in this step
		tabix = "imputation_runs/{run_name}/eagle_merged/{run_name}.chr{chr}.vcf.gz.tbi"
	shell: #This is for reference-based phasing in Eagle
	#Is allowRefAltSwap necessary here? What is that doing. Doublecheck.
		"""
		module load bcftools
		psrecord "~/Eagle_v2.4.1/eagle --vcfRef {input.refvcf} --vcfTarget {input.vcf} --geneticMapFile {params.imputemap} --allowRefAltSwap --chromX {params.chrom} --numThreads {params.threads} --outPrefix {params.out}; tabix {output.vcf}" --log {params.psrecord} --include-children --interval 5
		"""

rule imputation: #A single round of imputation for all target assays.
	input: #Reference files are cross-imputed HD/F250 assays for all individuals genotyped on either or both assays
		ref = expand("reference_build/{ref_version}/reference/combined/bigref.850k.chr{{chr}}.m3vcf.gz",
		ref_version = config["ref_version"]),
		haps = "imputation_runs/{run_name}/eagle_merged/{run_name}.chr{chr}.vcf.gz", #Phased VCF of target assays
		hapstbi = "imputation_runs/{run_name}/eagle_merged/{run_name}.chr{chr}.vcf.gz.tbi",
		#rmap = "phasing_maps/imputemap.chr{chr}.map" #Only needed when imputing with Minimac4
	params:
		oprefix = "imputation_runs/{run_name}/minimac_imputed/{run_name}.chr{chr}",
		chrom = "{chr}",
		threads = config["mm_threads"],
		psrecord = "log/{run_name}/psrecord/imputation/imputation.chr{chr}.log"
	output:
		gen = "imputation_runs/{run_name}/minimac_imputed/{run_name}.chr{chr}.dose.vcf.gz", #Outputs VCF file with both dosage and hardcall information
		#tbi = "imputation_runs/{run_name}/minimac_imputed/{run_name}.chr{chr}.dose.vcf.gz.tbi"
		#gen = "minimac_imputed/combined_imputed/mm4/{run_name}.chr{chr}.dose.vcf.gz"
		#vcf = "minimac_imputed/combined_imputed/{run_name}.chr{chr}.m3vcf.gz"
	shell: #Minimac3 appears to be working better, not sure what the hangup with Minimac4 is, but will explore in the near future
		"""
		module load bcftools
		psrecord "/storage/hpc/group/UMAG/SCRIPTS/Minimac4-1.0.2/release-build/minimac4 --refHaps {input.ref} --haps {input.haps} --allTypedSites --myChromosome {params.chrom} --cpu {params.threads} --prefix {params.oprefix};tabix {output.gen}" --log {params.psrecord} --include-children --interval 5

		"""

# #Minimac's dosage conversion does not work at this point. If it ever does, this'll be the rule that makes it work
# rule dosage_convert:
# 	input:
# 		gen = "minimac_imputed/{run_name}.chr{chr}.dose.vcf.gz",
# 		info = "minimac_imputed/{run_name}.chr{chr}.info"
# 	threads: 2
# 	params:
# 		oprefix = "minimac_imputed/combined_imputed/{run_name}.chr{chr}",
# 		chrom = "{chr}"
# 	log:
# 		"logs/dosage_convert/{run_name}.chr{chr}.log"
# 	benchmark:
# 		"benchmarks/dosage_convert/{run_name}.chr{chr}.benchmark.txt"
# 	output:
# 		dosage = "minimac_imputed/combined_imputed/{run_name}.chr{chr}.plink.dosage.gz",
# 		fam = "minimac_imputed/combined_imputed/{run_name}.chr{chr}.plink.fam",
# 		map = "minimac_imputed/combined_imputed/{run_name}.chr{chr}.plink.map"
# 	shell:
# 		"(DosageConvertor --vcfDose {input.gen} --info {input.info} --MyChromosome {params.chrom} --prefix {params.oprefix}) > {log}"

#VCF files need to be ordered for conversion to work correctly
#If this step has issues, it may be that the perl library isn't pointing to the correct vcftools install
#https://www.biostars.org/p/15163/
rule order_vcfs:
	input:
		vcf = "imputation_runs/{run_name}/minimac_imputed/{run_name}.chr{chr}.dose.vcf.gz", #each single-chromosome VCF file needs reordered
		template = "imputation_runs/{run_name}/minimac_imputed/{run_name}.chr29.dose.vcf.gz" #Chr29 is used as the "template" order
	params:
		vcf = "imputation_runs/{run_name}/imputed_genotypes/single_chrom/{run_name}.chr{chr}.reordered.vcf",
		psrecord = "log/{run_name}/psrecord/order_vcfs/order_vcfs.chr{chr}.log"
	output:
		vcf = temp("imputation_runs/{run_name}/imputed_genotypes/single_chrom/{run_name}.chr{chr}.reordered.vcf.gz"),
		tabix = temp("imputation_runs/{run_name}/imputed_genotypes/single_chrom/{run_name}.chr{chr}.reordered.vcf.gz.tbi")
	shell: #shuffle-cols does exactly what we need it to. Then bgzip and tabix output for concatenation with bcftools
		"""
		module load bcftools
		export PERL5LIB=bin/vcftools_0.1.13/perl/
		psrecord "bin/vcftools_0.1.13/perl/vcf-shuffle-cols -t {input.template} {input.vcf} > {params.vcf}; bgzip {params.vcf}; tabix {output.vcf}" --log {params.psrecord} --include-children --interval 5
		"""

# rule gwas_format: #Have changed this to do the longest step on a single chromosome basis, can't believe it took me this long...
# 	input:
# 		vcf = "imputation_runs/{run_name}/minimac_imputed/combined_reordered/{run_name}.chr{chr}.reordered.vcf.gz" #Reordered vcf file
#
# 	output:
# 		hard = "imputation_runs/{run_name}/imputed_genotypes_single_chrom/{run_name}.chr{chr}.hard.mgf.gz", #Hard calls for GEMMA in mgf format
# 		dose = "imputation_runs/{run_name}/imputed_genotypes_single_chrom/{run_name}.chr{chr}.dose.mgf.gz", #Dosage genotypes in mgf format
# 		sample = "imputation_runs/{run_name}/imputed_genotypes_single_chrom/{run_name}.chr{chr}.sample", #Sample file identifies the ordering of individuals in the genotype files
# 		#vcf = "minimac_imputed/850k_imputed/{run_name}.chr{chr}.dose.vcf.gz"
# 	shell: #My own custom script. Pretty damned slow, but it works. Breaking up into chromosomes should speed it WAY up. UPDATE it does
# 		"python bin/mm_to_mgf.py {input.vcf} {output.hard} {output.dose} {output.sample}"

rule hardcall_vcf:
	input:
		vcf = "imputation_runs/{run_name}/imputed_genotypes/single_chrom/{run_name}.chr{chr}.reordered.vcf.gz"
	params:
		vcf = "imputation_runs/{run_name}/imputed_genotypes/single_chrom/{run_name}.chr{chr}.hardcall.vcf",
		psrecord = "log/{run_name}/psrecord/hardcall_vcf/hardcall_vcf.chr{chr}.log"
	output:
		vcf = "imputation_runs/{run_name}/imputed_genotypes/single_chrom/{run_name}.chr{chr}.hardcall.vcf.gz"
	shell:#This script pulls out only the hardcall vcf information from the minimac imputation output, keeps vcf format then bgzips and tabix
		"""
		module load bcftools
		psrecord "python bin/vcf_hardcall_conversion.py {input.vcf} {params.vcf}; bgzip {params.vcf}; tabix {output.vcf}" --log {params.psrecord} --include-children --interval 5"""

rule concat_hardcall_vcf:#Puts individual chromosome files back into a single VCF
	input:
		concat = expand("imputation_runs/{{run_name}}/imputed_genotypes/single_chrom/{{run_name}}.chr{chr}.hardcall.vcf.gz", chr = list(range(1,32)))
	params:
		psrecord = "log/{run_name}/psrecord/concat_hardcall_vcf/concat_hardcall_vcf.log"
	output:
		vcf = "imputation_runs/{run_name}/imputed_genotypes/{run_name}.hardcall.vcf.gz",
		tbi = "imputation_runs/{run_name}/imputed_genotypes/{run_name}.hardcall.vcf.gz.tbi"
	shell:
		"""
		module load bcftools
		psrecord "bcftools concat {input.concat} -O z -o {output.vcf}; tabix {output.vcf}"  --log {params.psrecord} --include-children --interval 5
		"""
#Use BCFtools to combine individual sorted VCF files
rule concat_vcf:
	input:
		vcf = expand("imputation_runs/{{run_name}}/imputed_genotypes/single_chrom/{{run_name}}.chr{chr}.reordered.vcf.gz", chr = list(range(1,32))),
		tabix = expand("imputation_runs/{{run_name}}/imputed_genotypes/single_chrom/{{run_name}}.chr{chr}.reordered.vcf.gz.tbi", chr = list(range(1,32)))
	params:
		psrecord = "log/{run_name}/psrecord/concat_vcf/concat_vcf.log"
	output:
		vcf = "imputation_runs/{run_name}/imputed_genotypes/{run_name}.vcf.gz",
		tabix = "imputation_runs/{run_name}/imputed_genotypes/{run_name}.vcf.gz.tbi",
	shell: #List of all single-chromosome vcf files will be put onto the command line here for the concat command
		"""
		module load bcftools
		psrecord "bcftools concat {input.vcf} -O z -o {output.vcf}; tabix {output.vcf}" --log {params.psrecord} --include-children --interval 5"""

# rule concat_mgf: #Just a basic concatentation of the MGF files
# 	input:
# 		hard = expand("{{run_name}}/imputed_genotypes_single_chrom/{{run_name}}.chr{chr}.hard.mgf.gz", chr = list(range(1,30))),
# 		dose = expand("{{run_name}}/imputed_genotypes_single_chrom/{{run_name}}.chr{chr}.dose.mgf.gz", chr = list(range(1,30))),
# 		sample = "imputed_genotypes_single_chrom/{run_name}.chr29.sample" #All individual chromosome sample files are the same, this is the one that'll get copied
# 	params:
# 		hard = "imputation_runs/{run_name}/imputed_genotypes/{run_name}.hard.mgf", #These are the files that zcat will get pointed to before compression
# 		dose = "imputation_runs/{run_name}/imputed_genotypes/{run_name}.dose.mgf"
# 	output:
# 		hard = "imputation_runs/{run_name}/imputed_genotypes/{run_name}.hard.mgf.gz",
# 		dose = "imputation_runs/{run_name}/imputed_genotypes/{run_name}.dose.mgf.gz",
# 		sample = "imputation_runs/{run_name}/imputed_genotypes/{run_name}.sample",
# 	shell: #1. hardcall files    2. doseage files    3. Move sample file
# 		"zcat {input.hard} > {params.hard}; pigz {params.hard}; zcat {input.dose} > {params.dose}; pigz {params.dose}; mv {input.sample} > {output.sample}"
