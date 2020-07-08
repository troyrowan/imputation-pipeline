rule imputed_seq:
	input:
		# vcf = expand("imputed/{run_name}/{run_name}.chr{chr}.dose.vcf.gz",
		# run_name = config["run_name"],
		# chr = list(range(1,30)))
		vcf = expand("imputed/{run_name}/{run_name}.chr{chr}.mach.dose.gz",
		run_name = config["run_name"],
		chr = list(range(1,30)))

#To run this on Lewis, copy and paste the following:
#snakemake -s seq_imputation.snakefile --configfile config/190402_RAN.config.json --cluster-config 190930_imputation_config.json --cluster "sbatch -p {cluster.p} -o {cluster.o} --account {cluster.account} -t {cluster.t} -c {cluster.c} -n {cluster.n} --mem {cluster.mem}" --rerun-incomplete -np

rule imputation:
	input:
		gen = "imputation_runs/{run_name}/imputed_genotypes/single_chrom/{run_name}.chr{chr}.reordered.vcf.gz",
		ref = expand("reference_build/{run_name}/reference/{run_name}.{filter}.chr{{chr}}.m3vcf.gz",
		run_name = config["run_name"],
		filter = config["tranche"]+"_"+config["allele_count"])
	params:
		oprefix = "imputation_runs/{run_name}/seq_imputed/{run_name}/{run_name}.chr{chr}",
		chrom = "{chr}",
		threads = config["mm_threads"]
	output:
		vcf = "imputation_runs/{run_name}/seq_imputed/{run_name}/{run_name}.chr{chr}.dose.vcf.gz"
	shell:
		"/storage/hpc/group/UMAG/SCRIPTS/Minimac4-1.0.2/release-build/minimac4 --refHaps {input.ref} --haps {input.gen} --cpus {params.threads} --myChromosome {params.chrom} --prefix {params.oprefix}"


rule convert_mach:
	input:
		vcf = "imputation_runs/{run_name}/seq_imputed/{run_name}/{run_name}.chr{chr}.dose.vcf.gz",
		info = "imputation_runs/{run_name}/seq_imputed/{run_name}/{run_name}.chr{chr}.info"
	params:
		chrom = "{chr}",
		oprefix = "imputation_runs/{run_name}/imputed/{run_name}/{run_name}.chr{chr}"
	output:
		mach = "imputation_runs/{run_name}/seq_imputed/{run_name}/{run_name}.chr{chr}.mach.dose.gz",
		info = "imputation_runs/{run_name}/seq_imputed/{run_name}/{run_name}.chr{chr}.mach.info"
	shell:
		"/storage/hpc/group/UMAG/SCRIPTS/DosageConvertor/build/DosageConvertor --vcfDose {input.vcf} --info {input.info} --myChromosome {params.chrom} --type mach --format 1 --prefix {params.oprefix}"
