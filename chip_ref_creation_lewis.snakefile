import os
# Make log directories if they don't exist
for x in expand("log/{run_name}/slurm_out/{rule}", run_name = config['run_name'], rule = config['rules']):
	os.makedirs(x, exist_ok = True)
for x in expand("log/{run_name}/psrecord/{rule}", run_name = config['run_name'], rule = config['rules']):
	os.makedirs(x, exist_ok = True)

rule targ:
	input:
		#refs=expand("reference/{run_name}/{sample}.chr{chr}.vcf.gz.tbi", sample=config["sample"], date=config["date"], chr=list(range(1,30)))
		ref=expand("imputation_runs/{run_name}/reference/combined/bigref.850k.chr{chr}.vcf.gz", run_name=config["run_name"], chr=list(range(1,30)))
	shell:
		"rm .snakemake/*_tracking/*"
include: "qc.snakefile"

#This rule will create single-chromosome vcf files of the reference
rule recode_vcf:
	input:
		bim="imputation_runs/{run_name}/merged_files/{run_name}.bim",
		fam="imputation_runs/{run_name}/merged_files/{run_name}.fam",
		log="imputation_runs/{run_name}/merged_files/{run_name}.log",
		bed="imputation_runs/{run_name}/merged_files/{run_name}.bed"
	params:
		chrom="{chr}",
		inprefix="imputation_runs/{run_name}/merged_files/{run_name}",
		oprefix="imputation_runs/{run_name}/merged_files/{run_name}.chr{chr}",
		vcf="imputation_runs/{run_name}/merged_files/{run_name}.chr{chr}.vcf",
		threads=config["plink_threads"],
		mem=config["plink_mem"]
	output:
		vcf=temp("imputation_runs/{run_name}/merged_files/{run_name}.chr{chr}.vcf.gz"),
		tabix=temp("imputation_runs/{run_name}/merged_files/{run_name}.chr{chr}.vcf.gz.tbi"),
	shell:
		#This step is an issue in within-breed situatoins. Not sure why the --maf works for bigref and not for within breed, but leads to there being SNPs that are actually fixed, but are guessed to have an 01 genotype in eagle (and thus what is actually T . can't handle knowing that there is a 01)
		"plink --bfile {params.inprefix} --nonfounders --chr {params.chrom} --cow --memory {params.mem} --threads {params.threads} --real-ref-alleles --recode vcf --maf 0.0000001 --out {params.oprefix}; bgzip {params.vcf}; tabix {output.vcf}"
		#"(plink --bfile {params.inprefix} --nonfounders --chr {params.chrom} --chr-set 33 --memory 500 --maf 0.0001 --real-ref-alleles --make-bed --out {params.oprefix}) > {log}"

rule eagle_merged:
	input:
		vcf="imputation_runs/{run_name}/merged_files/{run_name}.chr{chr}.vcf.gz",
		map="genetic_map_1cMperMb.txt" #Map should be changed
	params:
		iprefix="imputation_runs/{run_name}/merged_files/{run_name}.chr{chr}.merged",
		out="imputation_runs/{run_name}/eagle_merged/bigref.chr{chr}.phased",
		threads=config["eagle_threads"]
	output:
		vcf="imputation_runs/{run_name}/eagle_merged/bigref.chr{chr}.phased.vcf.gz",
		tabix="imputation_runs/{run_name}/eagle_merged/bigref.chr{chr}.phased.vcf.gz.tbi"
	shell:
		"eagle --vcf {input.vcf} --geneticMapFile {input.map} --numThreads {params.threads} --chromX 32 --outPrefix {params.out}; tabix {output.vcf}"

rule make_phasing_vcf_extract_lists:
	input:
		bim=expand("imputation_runs/{{run_name}}/hwe_filtered/{sample}.bim", sample=config["ref_assays"]),
		fam=expand("imputation_runs/{{run_name}}/hwe_filtered/{sample}.fam", sample=config["ref_assays"])
	output:
		keep_ids="imputation_runs/{run_name}/vcf_per_assay/{sample}.keepvcf",
		keep_snps="imputation_runs/{run_name}/vcf_per_assay/{sample}.vcfregion"
	shell:
		"python bin/vcf_extraction_maker.py {input.bim} {input.fam} {output.keep_snps} {output.keep_ids}"

# rule fix_alleles:
# 	input:
# 		vcfgz="eagle_merged/{run_name}/bigref.chr{chr}.phased.vcf.gz",
# 		index="eagle_merged/{run_name}/bigref.chr{chr}.phased.vcf.gz.tbi"
# 	params:
# 		vcf="eagle_merged/{run_name}/bigref.chr{chr}.phased.fixed.vcf"
# 	output:
# 		vcfgz="eagle_merged/{run_name}/bigref.chr{chr}.phased.fixed.vcf.gz",
# 		index="eagle_merged/{run_name}/bigref.chr{chr}.phased.fixed.vcf.gz.tbi"
# 	shell:
# 		"zcat"
# 		#"zcat {input.vcfgz} | awk "{{OFS="\\t"}};{{gsub(".","A",\$5)}}1; bgzip {params.vcf}; tabix {output.vcf}"

rule refvcf_per_assay: #filter the vcfs on a per assay basis
	input:
		# vcfgz="imputation_runs/{run_name}/eagle_merged/bigref.chr{chr}.phased.fixed.vcf.gz",
		# index="imputation_runs/{run_name}/eagle_merged/bigref.chr{chr}.phased.fixed.vcf.gz.tbi",
		vcfgz="imputation_runs/{run_name}/eagle_merged/bigref.chr{chr}.phased.vcf.gz",
		index="imputation_runs/{run_name}/eagle_merged/bigref.chr{chr}.phased.vcf.gz.tbi",
		keep_ids="imputation_runs/{run_name}/vcf_per_assay/{sample}.keepvcf", #This iteration of the pipeline is only for testing Minimac's imputation accuracy
		#keep_ids="extract_lists/{run_name}/{sample}.chr{chr}.keepvcf",
		keep_maps="imputation_runs/{run_name}/vcf_per_assay/{sample}.chr{chr}.vcfregion"
	params:
		vcf="imputation_runs/{run_name}/vcf_per_assay/{sample}.chr{chr}.vcf"
	output:
		vcf="imputation_runs/{run_name}/vcf_per_assay/{sample}.chr{chr}.vcf.gz",
		tbi="imputation_runs/{run_name}/vcf_per_assay/{sample}.chr{chr}.vcf.gz.tbi"
	shell:
		"bcftools view {input.vcfgz} -R {input.keep_maps} -S {input.keep_ids} --force-samples -O z -o {output.vcf}; tabix {output.vcf}"

rule mm4_convert:
	input:
		#vcf=expand("phased_refs/{run_name}/{sample}.chr{{chr}}.vcf.gz", ref_version=config["ref_version"], sample=config["refs"]),
		#tbi=expand("phased_refs/{run_name}/{sample}.chr{{chr}}.vcf.gz.tbi", ref_version=config["ref_version"], sample=config["refs"])
		#vcf="phased_refs/{run_name}/{sample}.chr{chr}.vcf.gz",
		#tbi="phased_refs/{run_name}/{sample}.chr{chr}vcf.gz.tbi"
		vcf="imputation_runs/{run_name}/vcf_per_assay/{sample}.chr{chr}.vcf.gz",
		tbi="imputation_runs/{run_name}/vcf_per_assay/{sample}.chr{chr}.vcf.gz.tbi"
	params:
		oprefix="imputation_runs/{run_name}/reference/{sample}.chr{chr}",
		chrom="{chr}"
	output:
		hd="imputation_runs/{run_name}/reference/{sample}.chr{chr}.m3vcf.gz"
	shell:
		"Minimac3-omp --refHaps {input.vcf} --processReference --cpu 5 --myChromosome {params.chrom} --prefix {params.oprefix}"

rule refcreation_hd:
	input:
		f250=expand("imputation_runs/{run_name}/imputation_runs/{run_name}/reference/{sample}.chr{{chr}}.m3vcf.gz",
		run_name=config["run_name"],
		sample=config["f250ref"]),
		hd=expand("imputation_runs/{run_name}/vcf_per_assay/{sample}.chr{{chr}}.vcf.gz",
		run_name=config["run_name"],
		sample=config["hdref"]),
		hdtbi=expand("imputation_runs/{run_name}/vcf_per_assay/{sample}.chr{{chr}}.vcf.gz.tbi",
		run_name=config["run_name"],
		sample=config["hdref"])
	params:
		hdimputedprefix="imputation_runs/{run_name}/reference/crossimp/hd.850k.chr{chr}",
		f250imputedprefix="imputation_runs/{run_name}/reference/crossimp/f250.850k.chr{chr}",
		chrom ="{chr}"
	output:
		hdimputed="imputation_runs/{run_name}/reference/crossimp/hd.850k.chr{chr}.dose.vcf.gz",
		sorted="imputation_runs/{run_name}/reference/crossimp/hd.850k.chr{chr}.vcf.gz",
		hdimputedtbi="imputation_runs/{run_name}/reference/crossimp/hd.850k.chr{chr}.vcf.gz.tbi"
	shell:
		"minimac4 --refHaps {input.f250} --haps {input.hd} --myChromosome {params.chrom} --format GT --cpu 5 --allTypedSites --prefix {params.hdimputedprefix}; bcftools sort {output.hdimputed} -O z -o {output.sorted};tabix {output.sorted}"

rule refcreation_f250:
	input:
		hd=expand("imputation_runs/{run_name}/reference/{sample}.chr{{chr}}.m3vcf.gz",
		run_name=config["run_name"],
		sample=config["hdref"]),
		f250tbi=expand("imputation_runs/{run_name}/vcf_per_assay/{sample}.chr{{chr}}.vcf.gz.tbi",
		run_name=config["run_name"],
		sample=config["f250ref"]),
		f250=expand("imputation_runs/{run_name}/vcf_per_assay/{sample}.chr{{chr}}.vcf.gz",
		run_name=config["run_name"],
		sample=config["f250ref"]),
	params:
		hdimputedprefix="imputation_runs/{run_name}/reference/crossimp/hd.850k.chr{chr}",
		f250imputedprefix="imputation_runs/{run_name}/reference/crossimp/f250.850k.chr{chr}",
		chrom ="{chr}"
	output:
		f250imputed="imputation_runs/{run_name}/reference/crossimp/f250.850k.chr{chr}.dose.vcf.gz",
		sorted="imputation_runs/{run_name}/reference/crossimp/f250.850k.chr{chr}.vcf.gz",
		f250imputedtabix="imputation_runs/{run_name}/reference/crossimp/f250.850k.chr{chr}.vcf.gz.tbi"

	shell:
		"minimac4 --refHaps {input.hd} --haps {input.f250} --myChromosome {params.chrom} --format GT --cpu 5 --allTypedSites --prefix {params.f250imputedprefix}; bcftools sort {output.f250imputed} -O z -o {output.sorted}; tabix {output.sorted}"

rule combinerefs:
	input:
		hdimputed="imputation_runs/{run_name}/reference/crossimp/hd.850k.chr{chr}.vcf.gz",
		f250imputed="imputation_runs/{run_name}/reference/crossimp/f250.850k.chr{chr}.vcf.gz"
	params:
		chrom="{chr}",
		oprefix="imputation_runs/{run_name}/reference/combined/bigref.850k.chr{chr}"
	output:
		ref="imputation_runs/{run_name}/reference/combined/bigref.850k.chr{chr}.vcf.gz",
		reftbi="imputation_runs/{run_name}/reference/combined/bigref.850k.chr{chr}.vcf.gz.tbi"
	shell:
		"vcf-merge {input.hdimputed} {input.f250imputed} | bgzip -c > {output.ref}; tabix {output.ref}"

rule reformat_refs:
	input:
		ref="imputation_runs/{run_name}/reference/combined/bigref.850k.chr{chr}.vcf.gz",
		reftbi="imputation_runs/{run_name}/reference/combined/bigref.850k.chr{chr}.vcf.gz.tbi"
	params:

	output:
		refm3vcf="imputation_runs/{run_name}/reference/combined/bigref.850k.chr{chr}.m3vcf.gz"
	shell:
		"Minimac3-omp --refHaps {output.ref} --processReference --myChromosome {params.chrom} --prefix {params.oprefix}"
