rule imputed_seq:
	input:
		# vcf = expand("imputed/{run_name}/{run_name}.chr{chr}.dose.vcf.gz",
		# run_name = config["run_name"],
		# chr = list(range(1,30)))
		vcf = expand("imputation_runs/{run_name}/seq_imputed/{run_name}.chr{chr}.filtered.vcf.gz",
		run_name = config["run_name"],
		chr = list(range(1,30)))

include: "chip_imp_lewis.snakefile"
#To run this on Lewis, copy and paste the following:
#snakemake -s seq_imputation.snakefile --configfile config/190402_RAN.config.json --cluster-config 190930_imputation_config.json --cluster "sbatch -p {cluster.p} -o {cluster.o} --account {cluster.account} -t {cluster.t} -c {cluster.c} -n {cluster.n} --mem {cluster.mem}" --rerun-incomplete -np

rule seq_imputation:
	input:
		gen = "imputation_runs/{run_name}/imputed_genotypes/single_chrom/{run_name}.chr{chr}.reordered.vcf.gz",
		ref = expand("reference_build/{ref}/reference/{ref}.{filter}.chr{{chr}}.m3vcf.gz",
		ref = config["seq_ref_version"],
		filter = config["filter"])
	params:
		oprefix = "imputation_runs/{run_name}/seq_imputed/{run_name}.chr{chr}",
		chrom = "{chr}",
		psrecord = "log/{run_name}/psrecord/seq_imputation/seq_imputation.{run_name}.chr{chr}.log",
		threads = config["mm_threads"]
	output:
		vcf = "imputation_runs/{run_name}/seq_imputed/{run_name}.chr{chr}.dose.vcf.gz",
		info = "imputation_runs/{run_name}/seq_imputed/{run_name}.chr{chr}.info"
	shell:
		"""
		psrecord "/storage/hpc/group/UMAG/SCRIPTS/Minimac4-1.0.2/release-build/minimac4 --refHaps {input.ref} --haps {input.gen} --cpus {params.threads} --myChromosome {params.chrom} --prefix {params.oprefix}" --log {params.psrecord} --interval 60 --include-children
		"""

rule extract_high_quality:
	input:
		vcf = "imputation_runs/{run_name}/seq_imputed/{run_name}.chr{chr}.dose.vcf.gz",
		info = "imputation_runs/{run_name}/seq_imputed/{run_name}.chr{chr}.info"
	params:
		MAF = config["MAF"],
		R2 = config["R2"]
	output:
		vcf = "imputation_runs/{run_name}/seq_imputed/{run_name}.chr{chr}.filtered.vcf.gz",
		tbi = "imputation_runs/{run_name}/seq_imputed/{run_name}.chr{chr}.filtered.vcf.gz.tbi",
	shell:
		"""
		module load bcftools
		bcftools view -i 'MAF>={params.MAF} & R2>={params.R2}' {input.vcf} -O z -o {output.vcf}
		tabix {output.vcf}
		"""

rule convert_sequence_mach:
	input:
		vcf = "imputation_runs/{run_name}/seq_imputed/{run_name}.chr{chr}.dose.vcf.gz",
		info = "imputation_runs/{run_name}/seq_imputed/{run_name}.chr{chr}.info"
	params:
		chrom = "{chr}",
		oprefix = "imputation_runs/{run_name}/seq_imputed/{run_name}.chr{chr}",
		psrecord = "log/{run_name}/psrecord/convert_seq_mach/convert_seq_mach.{run_name}.chr{chr}.log"
	output:
		mach = "imputation_runs/{run_name}/seq_imputed/{run_name}.chr{chr}.mach.dose.gz",
		info = "imputation_runs/{run_name}/seq_imputed/{run_name}.chr{chr}.mach.info"
	shell:
		"""
		psrecord "/storage/hpc/group/UMAG/SCRIPTS/DosageConvertor/build/DosageConvertor --vcfDose {input.vcf} --info {input.info} --myChromosome {params.chrom} --type mach --format 1 --prefix {params.oprefix}" --log {params.psrecord} --interval 60 --include-children
		"""
