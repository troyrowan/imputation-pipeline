#include: "chip_ref_creation_qc_lewis.snakefile"

import os
# Make log directories if they don't exist
for x in expand("log/{run_name}/slurm_out/{rule}", run_name = config['run_name'], rule = config['rules']):
	os.makedirs(x, exist_ok = True)
for x in expand("log/{run_name}/psrecord/{rule}", run_name = config['run_name'], rule = config['rules']):
	os.makedirs(x, exist_ok = True)

#snakemake -s chip_ref_creation_lewis.snakefile --jobs 1000 --rerun-incomplete --keep-going --latency-wait 30 --configfile chip_ref_creation_lewis.config.yaml --cluster-config chip_ref_creation_lewis.cluster.json --cluster "sbatch -p {cluster.p} -o {cluster.o} --account {cluster.account} -t {cluster.t} -c {cluster.c} --mem {cluster.mem} --account {cluster.account}" -np

rule targ:
	input:
		targ=lambda wildcards: expand("reference_build/{run_name}/reference/stats/bigref.850k.chr{chr}.stats.txt",
		run_name=config["run_name"],
		chr=list(range(1,31)))
		# targ=lambda wildcards: expand("imputation_runs/{run_name}/merged_files/{run_name}.chr{chr}.vcf.gz",
		# run_name=config["run_name"],
		# chr=list(range(1,31)))

		# targ=lambda wildcards: expand("reference_build/{run_name}/vcf_per_assay/{sample}.chr{chr}.vcf.gz",
		# run_name=config["run_name"],
		# chr=list(range(1,2)),
		# sample=config["ref_assays"])
		# targ=lambda wildcards: expand("reference_build/{run_name}/reference/{sample}.chr{chr}.m3vcf.gz",
		# run_name=config["run_name"],
		# chr=list(range(1,2)),
		# sample=config["ref_assays"])

		#"reference_build/{run_name}/reference/combined/bigref.850k.chr{chr}.m3vcf.gz"
		#refs=expand("reference/{run_name}/{sample}.chr{chr}.vcf.gz.tbi", sample=config["sample"], date=config["date"], chr=list(range(1,30)))
		# ref=expand("refer{run_name}/reference/combined/bigref.850k.chr{chr}.vcf.gz",
		# run_name=config["run_name"],
		# chr=list(range(1,32)))
	shell:
		"rm .snakemake/*_tracking/*"

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
		mem=config["plink_mem"],
		psrecord = "log/{run_name}/psrecord/recode_vcf/recode_vcf.chr{chr}.log"
	output:
		vcf=temp("imputation_runs/{run_name}/merged_files/{run_name}.chr{chr}.vcf.gz"),
		tabix=temp("imputation_runs/{run_name}/merged_files/{run_name}.chr{chr}.vcf.gz.tbi"),
	shell:
		#This step is an issue in within-breed situatoins. Not sure why the --maf works for bigref and not for within breed, but leads to there being SNPs that are actually fixed, but are guessed to have an 01 genotype in eagle (and thus what is actually T . can't handle knowing that there is a 01)
		"""
		module load plink
		module load bcftools
		psrecord "plink --bfile {params.inprefix} --nonfounders --chr {params.chrom} --cow --memory {params.mem} --threads {params.threads} --real-ref-alleles --recode vcf --maf 0.000000001 --out {params.oprefix}; bgzip {params.vcf}; tabix {output.vcf}" --log {params.psrecord} --include-children --interval 5
		"""
		#"(plink --bfile {params.inprefix} --nonfounders --chr {params.chrom} --chr-set 33 --memory 500 --maf 0.0001 --real-ref-alleles --make-bed --out {params.oprefix}) > {log}"
#Memory requirements from Eagle 2.4.1 documentation: When phasing without a reference panel, Eagle’s memory use scales linearly with the number of samples (N) and the number of SNPs (M). For our tests on N=150K UK Biobank samples, Eagle required ≈1 GB RAM per 1,000 SNPs.
rule eagle_merged:
	input:
		vcf="imputation_runs/{run_name}/merged_files/{run_name}.chr{chr}.vcf.gz",
		map="bin/genetic_map_1cMperMb.txt" #Map should be changed
	params:
		out="reference_build/{run_name}/eagle_merged/bigref.chr{chr}.phased",
		threads=config["eagle_threads"],
		psrecord = "log/{run_name}/psrecord/eagle_merged/eagle_merged.chr{chr}.log"
	output:
		vcf="reference_build/{run_name}/eagle_merged/bigref.chr{chr}.phased.vcf.gz",
		tabix="reference_build/{run_name}/eagle_merged/bigref.chr{chr}.phased.vcf.gz.tbi"
	shell:
		"""
		module load bcftools
		psrecord " ~/Eagle_v2.4.1/eagle --vcf {input.vcf} --geneticMapFile {input.map} --numThreads {params.threads} --chromX 32 --outPrefix {params.out}; tabix {output.vcf}" --log {params.psrecord} --include-children --interval 5
		"""

# rule make_phasing_vcf_extract_lists:
# 	input:
# 		bim=expand("imputation_runs/{{run_name}}/hwe_filtered/{sample}.bim", sample=config["ref_assays"]),
# 		fam=expand("imputation_runs/{{run_name}}/hwe_filtered/{sample}.fam", sample=config["ref_assays"])
# 	params:
# 		psrecord = "log/{run_name}/psrecord/make_phasing_vcf_extract_lists/make_phasing_vcf_extract_lists.{sample}.log",
# 		bim="imputation_runs/{run_name}/hwe_filtered/{sample}.bim",
# 		fam="imputation_runs/{run_name}/hwe_filtered/{sample}.fam"
# 	output:
# 		keep_ids="reference_build/{run_name}/vcf_per_assay/{sample}.keepvcf",
# 		keep_snps="reference_build/{run_name}/vcf_per_assay/{sample}.vcfregion"
# 	shell:
# 		"""
# 		psrecord "python bin/vcf_extraction_maker.py {params.bim} {params.fam} {output.keep_snps} {output.keep_ids}" --log {params.psrecord} --include-children --interval 5
# 		"""

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
		# vcfgz="reference_build/{run_name}/eagle_merged/bigref.chr{chr}.phased.fixed.vcf.gz",
		# index="reference_build/{run_name}/eagle_merged/bigref.chr{chr}.phased.fixed.vcf.gz.tbi",
		vcfgz="reference_build/{run_name}/eagle_merged/bigref.chr{chr}.phased.vcf.gz",
		index="reference_build/{run_name}/eagle_merged/bigref.chr{chr}.phased.vcf.gz.tbi",
		keep_ids="imputation_runs/{run_name}/vcf_per_assay/{sample}.keepvcf",
		keep_maps="imputation_runs/{run_name}/vcf_per_assay/{sample}.vcfregion"
	params:
		vcf="reference_build/{run_name}/vcf_per_assay/{sample}.chr{chr}.vcf",
		psrecord = "log/{run_name}/psrecord/refvcf_per_assay/refvcf_per_assay.{sample}.chr{chr}.log"
	output:
		vcf="reference_build/{run_name}/vcf_per_assay/{sample}.chr{chr}.vcf.gz",
		tbi="reference_build/{run_name}/vcf_per_assay/{sample}.chr{chr}.vcf.gz.tbi"
	shell:
		"""
		module load bcftools
		psrecord "bcftools view {input.vcfgz} -R {input.keep_maps} -S {input.keep_ids} --force-samples -O z -o {output.vcf}; tabix {output.vcf}" --log {params.psrecord} --include-children --interval 5
		"""

rule mm4_convert:
	input:
		#vcf=expand("phased_refs/{run_name}/{sample}.chr{{chr}}.vcf.gz", ref_version=config["ref_version"], sample=config["refs"]),
		#tbi=expand("phased_refs/{run_name}/{sample}.chr{{chr}}.vcf.gz.tbi", ref_version=config["ref_version"], sample=config["refs"])
		#vcf="phased_refs/{run_name}/{sample}.chr{chr}.vcf.gz",
		#tbi="phased_refs/{run_name}/{sample}.chr{chr}vcf.gz.tbi"
		vcf="reference_build/{run_name}/vcf_per_assay/{sample}.chr{chr}.vcf.gz",
		tbi="reference_build/{run_name}/vcf_per_assay/{sample}.chr{chr}.vcf.gz.tbi"
	params:
		oprefix="reference_build/{run_name}/temp_reference/{sample}.chr{chr}",
		chrom="{chr}",
		psrecord = "log/{run_name}/psrecord/mm4_convert/mm4_convert.{sample}.chr{chr}.log"
	output:
		hd="reference_build/{run_name}/temp_reference/{sample}.chr{chr}.m3vcf.gz"
	shell:
		"""
		psrecord "/home/tnr343/Minimac3/bin/Minimac3-omp --refHaps {input.vcf} --processReference --cpu 5 --myChromosome {params.chrom} --prefix {params.oprefix}" --log {params.psrecord} --include-children --interval 5
		"""

rule refcreation_hd:
	input:
		f250=expand("reference_build/{run_name}/temp_reference/{sample}.chr{{chr}}.m3vcf.gz",
		run_name=config["run_name"],
		sample=config["f250ref"]),
		hd=expand("reference_build/{run_name}/vcf_per_assay/{sample}.chr{{chr}}.vcf.gz",
		run_name=config["run_name"],
		sample=config["hdref"]),
		hdtbi=expand("reference_build/{run_name}/vcf_per_assay/{sample}.chr{{chr}}.vcf.gz.tbi",
		run_name=config["run_name"],
		sample=config["hdref"])
	params:
		hdimputedprefix="reference_build/{run_name}/reference/crossimp/hd.850k.chr{chr}",
		f250imputedprefix="reference_build/{run_name}/reference/crossimp/f250.850k.chr{chr}",
		chrom ="{chr}",
		threads=config["mm_threads"],
		psrecord = "log/{run_name}/psrecord/refcreation_hd/refcreation_hd.chr{chr}.log"
	output:
		hdimputed="reference_build/{run_name}/reference/crossimp/hd.850k.chr{chr}.dose.vcf.gz",
		sorted="reference_build/{run_name}/reference/crossimp/hd.850k.chr{chr}.vcf.gz",
		hdimputedtbi="reference_build/{run_name}/reference/crossimp/hd.850k.chr{chr}.vcf.gz.tbi"
	shell:
		"""
		module load bcftools
		psrecord "/home/tnr343/Minimac4/release-build/minimac4 --refHaps {input.f250} --haps {input.hd} --myChromosome {params.chrom} --format GT --cpu {params.threads} --allTypedSites --prefix {params.hdimputedprefix}; bcftools sort {output.hdimputed} -O z -o {output.sorted};tabix {output.sorted}" --log {params.psrecord} --include-children --interval 5
		"""
#hd="reference_build/{run_name}/reference/{sample}.chr{chr}.m3vcf.gz"
rule refcreation_f250:
	input:
		hd=expand("reference_build/{run_name}/temp_reference/{sample}.chr{{chr}}.m3vcf.gz",
		run_name=config["run_name"],
		sample=config["hdref"]),
		f250tbi=expand("reference_build/{run_name}/vcf_per_assay/{sample}.chr{{chr}}.vcf.gz.tbi",
		run_name=config["run_name"],
		sample=config["f250ref"]),
		f250=expand("reference_build/{run_name}/vcf_per_assay/{sample}.chr{{chr}}.vcf.gz",
		run_name=config["run_name"],
		sample=config["f250ref"])
	params:
		hdimputedprefix="reference_build/{run_name}/reference/crossimp/hd.850k.chr{chr}",
		f250imputedprefix="reference_build/{run_name}/reference/crossimp/f250.850k.chr{chr}",
		chrom ="{chr}",
		threads=config["mm_threads"],
		psrecord = "log/{run_name}/psrecord/refcreation_f250/refcreation_f250.chr{chr}.log"
	output:
		f250imputed="reference_build/{run_name}/reference/crossimp/f250.850k.chr{chr}.dose.vcf.gz",
		sorted="reference_build/{run_name}/reference/crossimp/f250.850k.chr{chr}.vcf.gz",
		f250imputedtabix="reference_build/{run_name}/reference/crossimp/f250.850k.chr{chr}.vcf.gz.tbi"
	shell:
		"""
		module load bcftools
		psrecord "/home/tnr343/Minimac4/release-build/minimac4 --refHaps {input.hd} --haps {input.f250} --myChromosome {params.chrom} --format GT --cpu {params.threads} --allTypedSites --prefix {params.f250imputedprefix}; bcftools sort {output.f250imputed} -O z -o {output.sorted}; tabix {output.sorted}" --log {params.psrecord} --include-children --interval 5
		"""

rule combinerefs:
	input:
		hdimputed="reference_build/{run_name}/reference/crossimp/hd.850k.chr{chr}.vcf.gz",
		f250imputed="reference_build/{run_name}/reference/crossimp/f250.850k.chr{chr}.vcf.gz"
	params:
		chrom="{chr}",
		oprefix="reference_build/{run_name}/reference/combined/bigref.850k.chr{chr}",
		psrecord = "log/{run_name}/psrecord/combinerefs/combinerefs.bigref.850k.chr{chr}.log"
	output:
		ref="reference_build/{run_name}/reference/combined/bigref.850k.chr{chr}.vcf.gz",
		reftbi="reference_build/{run_name}/reference/combined/bigref.850k.chr{chr}.vcf.gz.tbi"
	shell:
		"""
		module load bcftools
		psrecord "bcftools merge {input.hdimputed} {input.f250imputed} --force-samples -O z -o {output.ref}; tabix {output.ref}" --log {params.psrecord} --include-children --interval 5
		"""
#Above rule will
rule reformat_refs:
	input:
		ref="reference_build/{run_name}/reference/combined/bigref.850k.chr{chr}.vcf.gz",
		reftbi="reference_build/{run_name}/reference/combined/bigref.850k.chr{chr}.vcf.gz.tbi"
	params:
		psrecord = "log/{run_name}/psrecord/reformat_refs/reformat_refs.chr{chr}.log",
		chrom = "{chr}",
		oprefix="reference_build/{run_name}/reference/combined/bigref.850k.chr{chr}"
	output:
		refm3vcf="reference_build/{run_name}/reference/combined/bigref.850k.chr{chr}.m3vcf.gz"
	shell:
		"""
		psrecord "/home/tnr343/Minimac3/bin/Minimac3-omp --refHaps {input.ref} --processReference --myChromosome {params.chrom} --prefix {params.oprefix}" --log {params.psrecord} --include-children --interval 5
		"""

rule chip_ref_stats:
	input:
		ref="reference_build/{run_name}/reference/combined/bigref.850k.chr{chr}.vcf.gz",
		reftbi="reference_build/{run_name}/reference/combined/bigref.850k.chr{chr}.vcf.gz.tbi"
	#params:
		#psrecord = "log/{run_name}/psrecord/chip_ref_stats/chip_ref_stats.chr{chr}.log"
	output:
		stats="reference_build/{run_name}/reference/stats/bigref.850k.chr{chr}.stats.txt"
	shell:
		"""
		module load bcftools
		bcftools stats {input.ref} > {output.stats}
		"""
