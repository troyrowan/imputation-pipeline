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
		gen = "chip_imputed/{run_name}/{run_name}.chr{chr}.hardcall.vcf.gz",
		ref = expand("/storage/htc/schnabellab/results/9913/wgs/1kbulls_ars1.2/reference_haps/sequence/{ref_run}.{filter}.chr{{chr}}.m3vcf.gz",
		ref_run = config["seqref"],
		filter = config["seqref_filter"])
	threads: 20
	# benchmark:
	# 	"benchmarks/{run_name}/imputation/{run_name}.chr{chr}.benchmark.txt"
	params:
		oprefix = "imputed/{run_name}/{run_name}.chr{chr}",
		chrom = "{chr}"
	output:
		vcf = "imputed/{run_name}/{run_name}.chr{chr}.dose.vcf.gz"
	shell:
		"/storage/hpc/group/UMAG/SCRIPTS/Minimac4-1.0.2/release-build/minimac4 --refHaps {input.ref} --haps {input.gen} --cpus {threads} --myChromosome {params.chrom} --prefix {params.oprefix}"


rule convert_mach:
	input:
		vcf = "imputed/{run_name}/{run_name}.chr{chr}.dose.vcf.gz",
		info = "imputed/{run_name}/{run_name}.chr{chr}.info"
	params:
		chrom = "{chr}",
		oprefix = "imputed/{run_name}/{run_name}.chr{chr}"
	output:
		mach = "imputed/{run_name}/{run_name}.chr{chr}.mach.dose.gz",
		info = "imputed/{run_name}/{run_name}.chr{chr}.mach.info"
	shell:
		"/storage/hpc/group/UMAG/SCRIPTS/DosageConvertor/build/DosageConvertor --vcfDose {input.vcf} --info {input.info} --myChromosome {params.chrom} --type mach --format 1 --prefix {params.oprefix}"
