##include: "bigrefprep.snakefile"
#include: "190122_refcreation.snakefile"
include: "qc.snakefile"
rule bigref_done: #Last file create is specified up here. Use expand to indicate how we want wild cards filled in.
	input:
		gen = expand("imputed_genotypes/{run_name}/{run_name}.vcf.gz", run_name = config["run_name"]),
		mgf = expand("imputed_genotypes/{run_name}/{run_name}.hardcall.vcf.gz", run_name = config["run_name"]),
		chroms = expand("imputed_genotypes_single_chrom/{run_name}/{run_name}.chr{chr}.dose.mgf.gz", run_name = config["run_name"], chr = list(range(1,30)))
		#dose = expand("minimac_imputed/{run_name}/combined_imputed/mm4/{run_name}.chr{chr}.dose.vcf.gz", date = config["date"], run_name = config["run_name"], chr = list(range(1,30)))
		#mgf = expand("imputed_chromosomes/{run_name}/{run_name}.chr{chr}.dose.mgf.gz", run_name = config["run_name"], chr = "29")
	shell:
		"rm .snakemake/*_tracking/*"

#Input is the filtered merged file from PLINK filtering. Single set of bed/bim/fam files needed for this step
#Splits into chromosomes for parallelization
rule assay_chrsplit:
	input:
		# bed = "subsetted_test/{run_name}/{run_name}.SNP50.bed",
		# bim = "subsetted_test/{run_name}/{run_name}.SNP50.bim",
		# fam = "subsetted_test/{run_name}/{run_name}.SNP50.fam",
		bed = "merged_files/{run_name}/{run_name}.bed",
		bim = "merged_files/{run_name}/{run_name}.bim",
 		fam = "merged_files/{run_name}/{run_name}.fam"
	params:
	#	inprefix = "subsetted_test/{run_name}/{run_name}.SNP50",
	 	inprefix = "merged_files/{run_name}/{run_name}",
		vcf = "assay_raw_vcf/{run_name}/{run_name}.chr{chr}.vcf",
		oprefix = "assay_raw_vcf/{run_name}/{run_name}.chr{chr}",
		chr = "{chr}"
	benchmark:
		"benchmarks/assay_chrsplit/{run_name}/{run_name}.chr{chr}.txt"
	log:
		"logs/assay_chrsplit/{run_name}.chr{chr}.log"
	output:
		vcf = temp("assay_raw_vcf/{run_name}/{run_name}.chr{chr}.vcf.gz"), #These files are temporary, will be removed when no longer needed by new input rule
		tabix = temp("assay_raw_vcf/{run_name}/{run_name}.chr{chr}.vcf.gz.tbi")
	shell: #PLINK is used to split merged files into 29 chromosomes and then recode as vcf, bgzipped post-recoding
		"(plink --bfile {params.inprefix} --keep-allele-order --chr {params.chr} --nonfounders --chr-set 29 --recode vcf --out {params.oprefix}; bgzip {params.vcf}; tabix {output.vcf})> {log}"

rule bigref_phasing:
	input: #This is using ALL HD individuals as the phasing reference (HD individuals are treated as true haplotypes)
		refvcf = expand("reference/{ref_version}/{sample}.chr{{chr}}.vcf.gz", ref_version = config["phasing_ref_version"], sample = config["hdref"]), #References what HD are called from config file, also which reference version to use
		#refvcf = expand("vcf_per_assay/{run_name}/{sample}.chr{{chr}}.vcf.gz", run_name = config["phasing_ref_version"], sample = config["hdref"]),
		vcf = "assay_raw_vcf/{run_name}/{run_name}.chr{chr}.vcf.gz", #Input file generated above
		tabix = "assay_raw_vcf/{run_name}/{run_name}.chr{chr}.vcf.gz.tbi"#Ensures tabix file was created for target genotypes
	params:
		imputemap = "/usr/local/bin/Eagle_v2.4/tables/genetic_map_1cMperMb.txt",
		out="eagle_merged/{run_name}/{run_name}.chr{chr}",
		chrx = "{chr}" #Chromosome wildcard for Eagle
	threads: 10 #Allow Eagle to use 10 cores per job, it will use this many, but not more. When we give 60 cores, 6 phasing jobs happen at a time
	priority: 30 #Is this necessary? or is this forcing everything to wait for phasing to end before imputation can start?
	benchmark:
		"benchmarks/bigref_phasing/{run_name}.chr{chr}.benchmark.txt"
	log:
		"logs/bigref_phasing/{run_name}.chr{chr}.log"
	output:
		vcf = "eagle_merged/{run_name}/{run_name}.chr{chr}.vcf.gz", #Eagle outputs a bgzipped VCF file, we tabix in this step
		tabix = "eagle_merged/{run_name}/{run_name}.chr{chr}.vcf.gz.tbi"
	shell: #This is for reference-based phasing in Eagle onlyself.
	#Is allowRefAltSwap necessary here? What is that doing. Doublecheckself.
		"(eagle --vcfRef {input.refvcf} --vcfTarget {input.vcf} --geneticMapFile {params.imputemap} --allowRefAltSwap --chromX {params.chrx} --numThreads 10 --outPrefix {params.out}; tabix {output.vcf})> {log}"

rule imputation: #A single round of imputation for all target assays.
	input: #Reference files are cross-imputed HD/F250 assays for all individuals genotyped on either or both assays
		ref = expand("reference/combined/{ref_version}/bigref.850k.chr{{chr}}.m3vcf.gz", ref_version = config["imp_ref_version"], sample = config["hdref"]),
		reftbi = expand("reference/combined/{ref_version}/bigref.850k.chr{{chr}}.m3vcf.gz", ref_version = config["imp_ref_version"], sample = config["hdref"]),
		haps = "eagle_merged/{run_name}/{run_name}.chr{chr}.vcf.gz", #Phased VCF of target assays
		hapstbi = "eagle_merged/{run_name}/{run_name}.chr{chr}.vcf.gz.tbi",
		rmap = "phasing_maps/imputemap.chr{chr}.map" #Only needed when imputing with Minimac4
	threads: 5 #Minimac rarely uses more than 5 cores, have set this as optimal allocation of resources.
	params:
		oprefix = "minimac_imputed/{run_name}/{run_name}.chr{chr}",
		chrom = "{chr}"
	log:
		"logs/imputation/{run_name}/{run_name}.chr{chr}.log"
	benchmark:
		"benchmarks/imputation/{run_name}/{run_name}.chr{chr}.benchmark.txt"
	output:
		gen = "minimac_imputed/{run_name}/{run_name}.chr{chr}.dose.vcf.gz", #Outputs VCF file with both dosage and hardcall information
		#gen = "minimac_imputed/{run_name}/combined_imputed/mm4/{run_name}.chr{chr}.dose.vcf.gz"
		#vcf = "minimac_imputed/{run_name}/combined_imputed/{run_name}.chr{chr}.m3vcf.gz"
	shell: #Minimac3 appears to be working better, not sure what the hangup with Minimac4 is, but will explore in the near future
		#"(minimac4 --refHaps {input.ref} --haps {input.haps} --mapFile {input.rmap} --allTypedSites --MyChromosome {params.chrom} --prefix {params.oprefix}) > {log}"
		"(Minimac3-omp --refHaps {input.ref} --haps {input.haps} --allTypedSites --myChromosome {params.chrom} --cpu 5 --prefix {params.oprefix}) > {log}"
#Issues occur when attempting to impute using minimac4, which should be MUCH faster than Minimac3. Fails as it goes to write imputed VCF files out
#Unsure of how to fix this, so for now, impute with Minimac3-omp (which allows multi-core imputation)

# #Minimac's dosage conversion does not work at this point. If it ever does, this'll be the rule that makes it work
# rule dosage_convert:
# 	input:
# 		gen = "minimac_imputed/{run_name}/{run_name}.chr{chr}.dose.vcf.gz",
# 		info = "minimac_imputed/{run_name}/{run_name}.chr{chr}.info"
# 	threads: 2
# 	params:
# 		oprefix = "minimac_imputed/{run_name}/combined_imputed/{run_name}.chr{chr}",
# 		chrom = "{chr}"
# 	log:
# 		"logs/dosage_convert/{run_name}/{run_name}.chr{chr}.log"
# 	benchmark:
# 		"benchmarks/dosage_convert/{run_name}/{run_name}.chr{chr}.benchmark.txt"
# 	output:
# 		dosage = "minimac_imputed/{run_name}/combined_imputed/{run_name}.chr{chr}.plink.dosage.gz",
# 		fam = "minimac_imputed/{run_name}/combined_imputed/{run_name}.chr{chr}.plink.fam",
# 		map = "minimac_imputed/{run_name}/combined_imputed/{run_name}.chr{chr}.plink.map"
# 	shell:
# 		"(DosageConvertor --vcfDose {input.gen} --info {input.info} --MyChromosome {params.chrom} --prefix {params.oprefix}) > {log}"

#VCF files need to be ordered for conversion to work correctly
rule order_vcfs:
	input:
		vcf = "minimac_imputed/{run_name}/{run_name}.chr{chr}.dose.vcf.gz", #each single-chromosome VCF file needs reordered
		template = "minimac_imputed/{run_name}/{run_name}.chr29.dose.vcf.gz" #Chr29 is used as the "template" order
	params:
		vcf = "minimac_imputed/{run_name}/combined_reordered/{run_name}.chr{chr}.reordered.vcf"
	log:
		"logs/order_vcfs/{run_name}/{run_name}.chr{chr}.log"
	benchmark:
		"benchmarks/order_vcfs/{run_name}/{run_name}.chr{chr}.benchmark.txt"
	output:
		vcf = temp("minimac_imputed/{run_name}/combined_reordered/{run_name}.chr{chr}.reordered.vcf.gz"),
		tabix = temp("minimac_imputed/{run_name}/combined_reordered/{run_name}.chr{chr}.reordered.vcf.gz.tbi")
	shell: #shuffle-cols does exactly what we need it to. Then bgzip and tabix output for concatenation with bcftools
		"(vcf-shuffle-cols -t {input.template} {input.vcf} > {params.vcf}; bgzip {params.vcf}; tabix {output.vcf}) > {log}"

rule gwas_format: #Have changed this to do the longest step on a single chromosome basis, can't believe it took me this long...
	input:
		vcf = "minimac_imputed/{run_name}/combined_reordered/{run_name}.chr{chr}.reordered.vcf.gz" #Reordered vcf file
	log:
		"logs/gwas_format/{run_name}/{run_name}.log"
	benchmark:
		"benchmarks/gwas_format/{run_name}/{run_name}.benchmark.txt"
	output:
		hard = "imputed_genotypes_single_chrom/{run_name}/{run_name}.chr{chr}.hard.mgf.gz", #Hard calls for GEMMA in mgf format
		dose = "imputed_genotypes_single_chrom/{run_name}/{run_name}.chr{chr}.dose.mgf.gz", #Dosage genotypes in mgf format
		sample = "imputed_genotypes_single_chrom/{run_name}/{run_name}.chr{chr}.sample", #Sample file identifies the ordering of individuals in the genotype files
		#vcf = "minimac_imputed/{run_name}/850k_imputed/{run_name}.chr{chr}.dose.vcf.gz"
	shell: #My own custom script. Pretty damned slow, but it works. Breaking up into chromosomes should speed it WAY up. UPDATE it does
		"(python bin/mm_to_mgf.py {input.vcf} {output.hard} {output.dose} {output.sample}) > {log}"

rule hardcall_vcf:
	input:
		vcf = "minimac_imputed/{run_name}/combined_reordered/{run_name}.chr{chr}.reordered.vcf.gz"
	params:
		vcf = "concat_vcf/{run_name}/{run_name}.chr{chr}.hardcall.vcf"
	log:
		"logs/hardcall_vcf/{run_name}/{run_name}.chr{chr}.log"
	benchmark:
		"benchmarks/hardcall_vcf/{run_name}/{run_name}.chr{chr}.benchmark.txt"
	output:
		vcf = "concat_vcf/{run_name}/{run_name}.chr{chr}.hardcall.vcf.gz"
	shell:
		"(python bin/vcf_hardcall_conversion.py {input.vcf} {params.vcf}; bgzip {params.vcf}; tabix {output.vcf}) > {log}"

rule concat_hardcall_vcf:
	input:
		concat = expand("concat_vcf/{{run_name}}/{{run_name}}.chr{chr}.hardcall.vcf.gz", chr = list(range(1,30)))
	log:
		"logs/concat_hardcall_vcf/{run_name}/{run_name}.log"
	benchmark:
		"benchmarks/concat_hardcall_vcf/{run_name}/{run_name}.benchmark.txt"
	output:
		vcf = "imputed_genotypes/{run_name}/{run_name}.hardcall.vcf.gz"
	shell:
		"(bcftools concat {input.concat} -O z -o {output.vcf})"
#Use BCFtools to combine individual sorted VCF files
rule concat_vcf:
	input:
		vcf = expand("minimac_imputed/{{run_name}}/combined_reordered/{{run_name}}.chr{chr}.reordered.vcf.gz", chr = list(range(1,30))),
		tabix = expand("minimac_imputed/{{run_name}}/combined_reordered/{{run_name}}.chr{chr}.reordered.vcf.gz.tbi", chr = list(range(1,30))),
	log:
		"logs/concat_vcf/{run_name}/{run_name}.log"
	benchmark:
		"benchmarks/concat_vcf/{run_name}/{run_name}.benchmark.txt"
	output:
		vcf = "imputed_genotypes/{run_name}/{run_name}.vcf.gz",
		tabix = "imputed_genotypes/{run_name}/{run_name}.vcf.gz.tbi",
	shell: #List of all single-chromosome vcf files will be put onto the command line here for the concat command
		"(bcftools concat {input.vcf} -O z -o {output.vcf}; tabix {output.vcf}) > {log}"

rule concat_mgf: #Just a basic concatentation of the MGF files
	input:
		hard = expand("imputed_genotypes_single_chrom/{{run_name}}/{{run_name}}.chr{chr}.hard.mgf.gz", chr = list(range(1,30))),
		dose = expand("imputed_genotypes_single_chrom/{{run_name}}/{{run_name}}.chr{chr}.dose.mgf.gz", chr = list(range(1,30))),
		sample = "imputed_genotypes_single_chrom/{run_name}/{run_name}.chr29.sample" #All individual chromosome sample files are the same, this is the one that'll get copied
	params:
		hard = "imputed_genotypes/{run_name}/{run_name}.hard.mgf", #These are the files that zcat will get pointed to before compression
		dose = "imputed_genotypes/{run_name}/{run_name}.dose.mgf"
	log:
		"logs/concat_mgf/{run_name}/{run_name}.log"
	benchmark:
		"benchmarks/concat_mgf/{run_name}/{run_name}.benchmark.txt"
	output:
		hard = "imputed_genotypes/{run_name}/{run_name}.hard.mgf.gz",
		dose = "imputed_genotypes/{run_name}/{run_name}.dose.mgf.gz",
		sample = "imputed_genotypes/{run_name}/{run_name}.sample",
	shell: #1. hardcall files    2. doseage files    3. Move sample file
		"(zcat {input.hard} > {params.hard}; pigz {params.hard}; zcat {input.dose} > {params.dose}; pigz {params.dose}; mv {input.sample} > {output.sample}) > {log}"
