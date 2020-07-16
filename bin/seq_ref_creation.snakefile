import os
# Make log directories if they don't exist
for x in expand("log/{run_name}/slurm_out/{rule}", run_name = config['run_name'], rule = config['rules']):
	os.makedirs(x, exist_ok = True)
for x in expand("log/{run_name}/psrecord/{rule}", run_name = config['run_name'], rule = config['rules']):
	os.makedirs(x, exist_ok = True)

rule targets:
	input:
		#phased = expand("eagle_phased/{run_name}.{filter}.chr{chr}.vcf.gz", run_name = "190805_seqref", filter = "T99", chr = list(range(1,29)))
		info = expand("imputation_runs/{run_name}/{run_name}.chr{chr}.info",
		run_name = config["run_name"],
		chr = list(range(1,30)) + ["X", "Y"])
		# filter = expand("reference_build/{run_name}/reference/{run_name}.{filter}.chr{chr}.m3vcf.gz",
		# run_name = config["run_name"],
		# filter = config["tranche"]+"_"+config["allele_count"],
		# chr = list(range(1,30)) + ["X", "Y"])
	# shell:
	# 	"rm .snakemake/*_tracking/*"

rule variant_stats:
	input:
		vcf = expand("{dir}Chr{{chr}}-Run8-TAUIND-raw-toDistribute.vcf.gz",
		dir = config["1kbulls_dir"])
	params:
		psrecord = "log/{run_name}/psrecord/variant_stats/variant_stats.chr{chr}.log"
	output:
		info = "imputation_runs/{run_name}/{run_name}.chr{chr}.info"
	shell:
		"""
		module load bcftools
		psrecord "bcftools query -f '%CHROM %POS %REF %ALT %AC %AF %AN %DP %ExcessHet %FS %InbreedingCoeff %MLEAC %MLEAF %MQ %QD %SOR %VQSLOD %culprit \n' {input.vcf} > {output.info}" --log {params.psrecord} --interval 5"""


rule tranche_extract:
	input:
		vcf = expand("{ref_dir}/Chr{{chr}}-Run7-TAU-tranche90-toDistribute.vcf.gz", ref_dir = config["1kbulls_dir"]),
		tbi = expand("{ref_dir}/Chr{{chr}}-Run7-TAU-tranche90-toDistribute.vcf.gz.tbi", ref_dir = config["1kbulls_dir"])
	params:
		vcf = "reference_build/{run_name}/tranche_extract/{run_name}.{filter}.chr{chr}.vcf",
		threads = config["bcftools_threads"],
		AC = config["allele_count"],
		psrecord = "log/{run_name}/psrecord/tranche_extract/tranche_extract.chr{chr}.log"
	output:
		vcf = "reference_build/{run_name}/tranche_extract/{run_name}.{filter}.chr{chr}.vcf.gz",
		tbi = "reference_build/{run_name}/tranche_extract/{run_name}.{filter}.chr{chr}.vcf.gz.tbi"
	shell:
		#"(/cluster/spack/opt/spack/linux-centos7-x86_64/gcc-4.8.5/bcftools-1.8-eexz77zeqrygzxeluq2fjenrwmeirovk/bin/bcftools view -i 'VQSLOD>=-2.8631' --max-alleles 2 -v snps {input.vcf} --threads 10 -O z -o {output.vcf}) > {log}"
		"""
		module load bcftools
		psrecord "bcftools view -i 'VQSLOD>=-2.8631 & AC>=20' --max-alleles 2 -v snps {input.vcf} --threads {params.threads} -O z -o {output.vcf}; tabix {output.vcf}" --log {params.psrecord} --include-children --interval 5
		"""

#BTA1 is 5.85% genome length, so that'll set the boundries for our memory/runtime when phasing
#Expecting BTA1 to have 5,537,259 SNPs (from T99 based on how many came out of BTA28)
#Based on Eagle Manual (150K individuals take 1GB per 1,000 SNPs) we should need 150GB for BTA1 (less for others)
# Can we split
rule phasing:
	input:
		vcf = "reference_build/{run_name}/tranche_extract/{run_name}.{filter}.chr{chr}.vcf.gz",
		tbi = "reference_build/{run_name}/tranche_extract/{run_name}.{filter}.chr{chr}.vcf.gz.tbi",
		map = "bin/genetic_map_1cMperMb.txt"
	params:
		oprefix = "reference_build/{run_name}/eagle_phased/{run_name}.{filter}.chr{chr}",
		threads = config["eagle_threads"],
		psrecord = "log/{run_name}/psrecord/phasing/phasing.chr{chr}.log"
	output:
		vcf = "reference_build/{run_name}/eagle_phased/{run_name}.{filter}.chr{chr}.vcf.gz",
		tbi = "reference_build/{run_name}/eagle_phased/{run_name}.{filter}.chr{chr}.vcf.gz.tbi"
	shell:
		"""
		module load bcftools
		psrecord "~/Eagle_v2.4.1/eagle --vcf {input.vcf} --geneticMapFile {input.map} --numThreads {params.threads} --chromX 32 --vcfOutFormat=z --outPrefix {params.oprefix}; tabix {output.vcf}" --log {params.psrecord} --include-children --interval 5
		"""

rule convert_mm:
	input:
		vcf = "reference_build/{run_name}/eagle_phased/{run_name}.{filter}.chr{chr}.vcf.gz",
		tbi = "reference_build/{run_name}/eagle_phased/{run_name}.{filter}.chr{chr}.vcf.gz.tbi"
	params:
		oprefix = "reference_build/{run_name}/reference/{run_name}.{filter}.chr{chr}",
		chrom = "{chr}",
		threads = config["mm_threads"],
		psrecord = "log/{run_name}/psrecord/phasing/phasing.chr{chr}.log"
	output:
		"reference_build/{run_name}/reference/{run_name}.{filter}.chr{chr}.m3vcf.gz"
	shell:
		"""
		psrecord "/storage/hpc/group/UMAG/SCRIPTS/Minimac3/bin/Minimac3-omp --refHaps {input.vcf} --processReference --cpu {params.threads} --lowMemory --myChromosome {params.chrom} --prefix {params.oprefix}" --log {params.psrecord} --include-children --interval 5
		"""
