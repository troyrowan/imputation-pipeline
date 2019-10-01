rule targ:
	input:
		#refs = expand("reference/{run_name}/{sample}.chr{chr}.vcf.gz.tbi", sample = config["sample"], date = config["date"], chr = list(range(1,30)))
		ref = expand("reference/combined/{run_name}/bigref.850k.chr{chr}.vcf.gz", run_name = config["run_name"], chr = list(range(1,30)))
	shell:
		"rm .snakemake/*_tracking/*"
#include: "qc.snakefile"
# rule testset_chrsplit:
# 	input:
# 		bed = "subsetted_test/{run_name}/{sample}.bed"
# 	params:
# 		inprefix = "subsetted_test/{run_name}/{sample}",
# 		oprefix = "assay_chrsplit/{run_name}/{sample}.chr{chr}",
# 		chr = "{chr}"
# 	benchmark:
# 		"benchmarks/assay_chrsplit/{run_name}/{sample}.chr{chr}.txt"
# 	log:
# 		"logs/eagle_split_chromosomes/{sample}.chr{chr}.log"
# 	output:
# 		bed = "assay_chrsplit/{run_name}/{sample}.chr{chr}.bed",
# 		bim = "assay_chrsplit/{run_name}/{sample}.chr{chr}.bim",
# 		fam = "assay_chrsplit/{run_name}/{sample}.chr{chr}.fam",
# 		log = "assay_chrsplit/{run_name}/{sample}.chr{chr}.log"
# 	shell:
# 		"(plink --bfile {params.inprefix}  --real-ref-alleles --chr {params.chr} --make-bed  --nonfounders --chr-set 29 --out {params.oprefix})> {log}"

# rule make_merge_list: #makes a merge list of raw genotypes, doesn't need to have a run because its just going to be a list related thing.
# 	input:
# 		filelist = expand("/{{run_name}}/{sample}.chr{{chr}}.bed", sample = config["sample"])
# 	log:
# 		"logs/make_merge_list/{run_name}/bigref.txt"
# 	output:
# 		"merged_files/{run_name}/bigref.chr{chr}.txt"
# 	shell:
# 		"(python ./bin/merge_file_maker.py {input.filelist} {output}) > {log}"
# # #
rule recode_vcf:
	input:
		bim = "merged_files/{run_name}/{run_name}.bim",
		fam = "merged_files/{run_name}/{run_name}.fam",
		log = "merged_files/{run_name}/{run_name}.log",
		bed = "merged_files/{run_name}/{run_name}.bed"
	params:
		chrom = "{chr}",
		inprefix="merged_files/{run_name}/{run_name}",
		oprefix = "merged_files/{run_name}/{run_name}.chr{chr}",
		vcf = "merged_files/{run_name}/{run_name}.chr{chr}.vcf"
	benchmark:
		"benchmarks/merged_chrsplit/{run_name}/{run_name}.chr{chr}.txt"
	log:
		"logs/merged_chrsplit/{run_name}/{run_name}.chr{chr}.log"
	output:
		vcf = temp("merged_files/{run_name}/{run_name}.chr{chr}.vcf.gz"),
		tabix = temp("merged_files/{run_name}/{run_name}.chr{chr}.vcf.gz.tbi"),
		#oprefix = "merged_files/{run_name}/{run_name}.chr{chr}.merged.bed"
	shell:
		#This step is an issue in within-breed situatoins. Not sure why the --maf works for bigref and not for within breed, but leads to there being SNPs that are actually fixed, but are guessed to have an 01 genotype in eagle (and thus what is actually T . can't handle knowing that there is a 01)
		"(plink --bfile {params.inprefix} --nonfounders --chr {params.chrom} --chr-set 33 --memory 500 --real-ref-alleles --recode vcf --out {params.oprefix}; bgzip {params.vcf}; tabix {output.vcf}) > {log}"
		#"(plink --bfile {params.inprefix} --nonfounders --chr {params.chrom} --chr-set 33 --memory 500 --maf 0.0001 --real-ref-alleles --make-bed --out {params.oprefix}) > {log}"
rule eagle_merged:
	input:
		vcf = "merged_files/{run_name}/{run_name}.chr{chr}.vcf.gz",
		#plink = "merged_files/{run_name}/{run_name}.chr{chr}.merged.bed",
		map = "/usr/local/bin/Eagle_v2.4/tables/genetic_map_1cMperMb.txt"
	params:
		iprefix = "merged_files/{run_name}/{run_name}.chr{chr}.merged",
		out="eagle_merged/{run_name}/bigref.chr{chr}.phased"
	threads: 10
	priority: 30
	benchmark:
		"benchmarks/eagle_merged/{run_name}/bigref.chr{chr}.benchmark.txt"
	log:
		"logs/eagle_merged/{run_name}/bigref.chr{chr}.log"
	output:
		vcf = "eagle_merged/{run_name}/bigref.chr{chr}.phased.vcf.gz",
		tabix = "eagle_merged/{run_name}/bigref.chr{chr}.phased.vcf.gz.tbi"
	shell:
		"(eagle --vcf {input.vcf} --geneticMapFile {input.map} --numThreads 10 --chromX 32 --outPrefix {params.out}; tabix {output.vcf})> {log}"

rule make_phasing_vcf_extract_lists: # This rule should be much easier when now that we're not trying to make lists talk back and forth.  Directly reference
	input:
		#bim = expand("assay_chrsplit/{{run_name}}/{sample}.chr{{chr}}.bim", sample = REFASSAYS),
		#fam = expand("assay_chrsplit/{{run_name}}/{sample}.chr{{chr}}.fam", sample = REFASSAYS)
		bim = "subsetted_test/{run_name}/{sample}.bim",
		fam = "subsetted_test/{run_name}/{sample}.fam"
	benchmark:
		"benchmarks/make_phasing_vcf_extract_lists/{run_name}/bigref.chr{chr}.benchmark.txt"
	log:
		"logs/make_phasing_vcf_extract_lists/{run_name}/bigref.chr{chr}.log"
	output:
		keep_ids = "phasing_extract_lists/{run_name}/{sample}.chr{chr}.keepvcf",
		keep_snps = "phasing_extract_lists/{run_name}/{sample}.chr{chr}.vcfregion"
	shell:
		"python bin/vcf_extraction_maker.py {input.bim} {input.fam} {output.keep_snps} {output.keep_ids}"

rule fix_alleles:
	input:
		vcfgz="eagle_merged/{run_name}/bigref.chr{chr}.phased.vcf.gz",
		index="eagle_merged/{run_name}/bigref.chr{chr}.phased.vcf.gz.tbi"
	params:
		vcf="eagle_merged/{run_name}/bigref.chr{chr}.phased.fixed.vcf"
	output:
		vcfgz="eagle_merged/{run_name}/bigref.chr{chr}.phased.fixed.vcf.gz",
		index="eagle_merged/{run_name}/bigref.chr{chr}.phased.fixed.vcf.gz.tbi"
	shell:
		"zcat"
		#"zcat {input.vcfgz} | awk "{{OFS="\\t"}};{{gsub(".","A",\$5)}}1; bgzip {params.vcf}; tabix {output.vcf}"

rule refvcf_per_assay: #filter the vcfs on a per assay basis
	input:
		vcfgz="eagle_merged/{run_name}/bigref.chr{chr}.phased.fixed.vcf.gz",
		index="eagle_merged/{run_name}/bigref.chr{chr}.phased.fixed.vcf.gz.tbi",
		keep_ids = "phasing_extract_lists/{run_name}/{sample}.chr{chr}.keepvcf", #This iteration of the pipeline is only for testing Minimac's imputation accuracy
		#keep_ids = "extract_lists/{run_name}/{sample}.chr{chr}.keepvcf",
		keep_maps = "phasing_extract_lists/{run_name}/{sample}.chr{chr}.vcfregion"
	params:
		vcf = "vcf_per_assay/{run_name}/{sample}.chr{chr}.vcf"
	benchmark:
		"benchmarks/vcf_per_assay/{run_name}/{sample}.chr{chr}.benchmark.txt"
	log:
		"logs/vcf_per_assay/{run_name}/{sample}.chr{chr}.log"
	output:
		vcf = "vcf_per_assay/{run_name}/{sample}.chr{chr}.vcf.gz",
		tbi = "vcf_per_assay/{run_name}/{sample}.chr{chr}.vcf.gz.tbi"
	shell:
		"(bcftools view {input.vcfgz} -R {input.keep_maps} -S {input.keep_ids} --force-samples -O z -o {output.vcf}; tabix {output.vcf}) > {log}"

rule mm4_convert:
	input:
		#vcf = expand("phased_refs/{run_name}/{sample}.chr{{chr}}.vcf.gz", ref_version = config["ref_version"], sample = config["refs"]),
		#tbi = expand("phased_refs/{run_name}/{sample}.chr{{chr}}.vcf.gz.tbi", ref_version = config["ref_version"], sample = config["refs"])
		#vcf = "phased_refs/{run_name}/{sample}.chr{chr}.vcf.gz",
		#tbi = "phased_refs/{run_name}/{sample}.chr{chr}vcf.gz.tbi"
		vcf = "vcf_per_assay/{run_name}/{sample}.chr{chr}.vcf.gz",
		tbi = "vcf_per_assay/{run_name}/{sample}.chr{chr}.vcf.gz.tbi"
	params:
		oprefix = "reference/{run_name}/{sample}.chr{chr}",
		chrom = "{chr}"
	threads: 5
	benchmark:
		"benchmarks/mm4_convert/{run_name}/combined.chr{chr}.benchmark.txt"
	log:
		"logs/mm4_convert/{run_name}/combined.chr{chr}.log"
	output:
		hd = "reference/{run_name}/{sample}.chr{chr}.m3vcf.gz"
	shell:
		"(Minimac3-omp --refHaps {input.vcf} --processReference --cpu 5 --myChromosome {params.chrom} --prefix {params.oprefix})>{log}"

rule refcreation_hd:
	input:
		f250 = expand("reference/{run_name}/{sample}.chr{{chr}}.m3vcf.gz", run_name = config["run_name"], sample = config["f250ref"]),
		hd = expand("vcf_per_assay/{run_name}/{sample}.chr{{chr}}.vcf.gz", run_name = config["run_name"], sample = config["hdref"]),
		hdtbi = expand("vcf_per_assay/{run_name}/{sample}.chr{{chr}}.vcf.gz.tbi", run_name = config["run_name"], sample = config["hdref"])
	params:
		hdimputedprefix = "reference/crossimp/{run_name}/hd.850k.chr{chr}",
		f250imputedprefix = "reference/crossimp/{run_name}/f250.850k.chr{chr}",
		chrom  = "{chr}"
	benchmark:
		"benchmarks/refcreation/{run_name}/combined.chr{chr}.benchmark.txt"
	log:
		"logs/refcreation/{run_name}/combined.chr{chr}.log"
	threads: 5
	output:
		hdimputed = "reference/crossimp/{run_name}/hd.850k.chr{chr}.dose.vcf.gz",
		sorted = "reference/crossimp/{run_name}/hd.850k.chr{chr}.vcf.gz",
		hdimputedtbi = "reference/crossimp/{run_name}/hd.850k.chr{chr}.vcf.gz.tbi"
	shell:
		"(minimac4 --refHaps {input.f250} --haps {input.hd} --myChromosome {params.chrom} --format GT --cpu 5 --allTypedSites --prefix {params.hdimputedprefix}; bcftools sort {output.hdimputed} -O z -o {output.sorted};tabix {output.sorted}) > {log}"

rule refcreation_f250:
	input:
		hd = expand("reference/{run_name}/{sample}.chr{{chr}}.m3vcf.gz", run_name = config["run_name"], sample = config["hdref"]),
		f250tbi = expand("vcf_per_assay/{run_name}/{sample}.chr{{chr}}.vcf.gz.tbi", run_name = config["run_name"], sample = config["f250ref"]),
		f250 = expand("vcf_per_assay/{run_name}/{sample}.chr{{chr}}.vcf.gz", run_name = config["run_name"], sample = config["f250ref"]),
	params:
		hdimputedprefix = "reference/crossimp/{run_name}/hd.850k.chr{chr}",
		f250imputedprefix = "reference/crossimp/{run_name}/f250.850k.chr{chr}",
		chrom  = "{chr}"
	benchmark:
		"benchmarks/refcreation/{run_name}/combined.chr{chr}.benchmark.txt"
	log:
		"logs/refcreation/{run_name}/combined.chr{chr}.log"
	threads: 5
	output:
		f250imputed = "reference/crossimp/{run_name}/f250.850k.chr{chr}.dose.vcf.gz",
		sorted = "reference/crossimp/{run_name}/f250.850k.chr{chr}.vcf.gz",
		f250imputedtabix = "reference/crossimp/{run_name}/f250.850k.chr{chr}.vcf.gz.tbi"

	shell:
		"(minimac4 --refHaps {input.hd} --haps {input.f250} --myChromosome {params.chrom} --format GT --cpu 5 --allTypedSites --prefix {params.f250imputedprefix}; bcftools sort {output.f250imputed} -O z -o {output.sorted}; tabix {output.sorted}) > {log}"

rule combinerefs:
	input:
		hdimputed = "reference/crossimp/{run_name}/hd.850k.chr{chr}.vcf.gz",
		f250imputed = "reference/crossimp/{run_name}/f250.850k.chr{chr}.vcf.gz"
	params:
		chrom = "{chr}",
		oprefix = "reference/combined/{run_name}/bigref.850k.chr{chr}"
	benchmark:
		"benchmarks/combinerefs/{run_name}/combined.chr{chr}.benchmark.txt"
	log:
		"logs/combinerefs/{run_name}/combined.chr{chr}.log"
	output:
		ref = "reference/combined/{run_name}/bigref.850k.chr{chr}.vcf.gz",
		reftbi = "reference/combined/{run_name}/bigref.850k.chr{chr}.vcf.gz.tbi",
		refm3vcf = "reference/combined/{run_name}/bigref.850k.chr{chr}.m3vcf.gz"
	shell:
		"(vcf-merge {input.hdimputed} {input.f250imputed} | bgzip -c > {output.ref}; tabix {output.ref}; Minimac3-omp --refHaps {output.ref} --processReference --myChromosome {params.chrom} --prefix {params.oprefix}) > {log}"
