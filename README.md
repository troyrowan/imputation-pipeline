# UMAG Imputation Pipeline
## Setup and basic pipeline instructions

#### Cloning repository and setting up environments 
Navigate to directory where you want to setup imputation pipeline and clone repository \
`git clone git@github.com:troyrowan/imputation-pipeline.git`

You'll need to activate a conda environment on Lewis that includes the relevant packages to run the pipeline.
This string of commands should load all of the packages that you'll need to run the pipeline.
`module load miniconda3; source activate /storage/hpc/data/tnr343/miniconda/envs/snake`
NOTE: We may want to change the path to this environment to be somewhere shared

Activate environment (This needs to be done every time prior to running Snakemake)\
`source activate snake`

### Setting up config files for imputation
Create copy of master config file for your imputation run\
`cp config/master_config.json config/my_imprun_config.json`

Add run-specific information into config file: 
* Only the following items should need changed
* Updates to reference should be reflected in updated master config (which you can get from a `git pull`)
* `sample`
    + List of `.ped` file prefixes in the format: `nsnps.date.nindividuals.manifest`
* `gt_path`
    + Absolute path to where raw ped files reside
* `run_name`
    + Whatever you choose your imputation run name to be
* Changes to filtering parameters can also be altered here, but traditional filtering is propogated through when master config is copied
    
Example top of config file is given below:

```python
rules: ["ref_alt", "ref_build_start", "no_duplicates", "variant_stats", "filter_variants", "individual_stats", "filter_individuals", "hwe_stats", "filter_hwe_variants", "filter_monomorphic", "filter_logging", "merge_assays", "assay_chrsplit", "bigref_phasing", "imputation", "order_vcfs", "hardcall_vcf", "concat_hardcall_vcf", "concat_vcf", "seq_imputation", "convert_seq_mach", "tabix", "convert_plink", "concat_plink"]
gt_path: "/storage/hpc/group/UMAG/WORKING/tnr343/imputation/chip_imputation/imputation-pipeline/imputation_runs/200921_HairShed/raw_genotypes/"
sample: ["139977.200921.42.C", "139977.200921.74.B"]
ref_assays: ["F250_all_ref", "HD_all_ref"]
ref_version: "200812_chip_refbuild"
seq_ref_version: "200708_Run8_refbuild"
hdref: "HD_all_ref"
f250ref: "F250_all_ref"
run_name: "200921_HairShed"
filter: "T90_20"
ind_callrate_filter: 0.1
snp_callrate_filter: 0.1
hwe_filter: 1e-75
mac_filter: 1
plink_threads: 4
plink_mem: 2000
eagle_threads: 16
mm_threads: 8
mapdict: {"26504":"/storage/hpc/group/UMAG/PLINK_FILES/snp_number/9913_ARS1.2_26504_GGPLDV3_snp_number_200806.map", ...}
refdict: {"26504":"/storage/hpc/group/UMAG/PLINK_FILES/snp_number/9913_ARS1.2_26504_GGPLDV3_snp_number_200806.REF", ...}
refalleledict: {"26504":"/storage/hpc/group/UMAG/PLINK_FILES/snp_number/9913_ARS1.2_26504_GGPLDV3_snp_number_200806.REF_ALLELE", ...}
}
```

**NOTE**
As new assays become available, the dictionaries `"mapdict"`, `"refdict"`, and `"refalleledict"` that refer to corresponding `.map`, `.REF`, and `.REF_ALLELE` 
files need to be updated as well. Where the key is `"nsnps"` and value is the absolute file path to each of the three files

### Running imputation pipeline (850K SNPs)
Upon setting up config file, execute a "dry-run" of Snakemake to ensure that it sees all necessary files, etc.\
`snakemake -s bin/chip_imp_lewis.snakefile --configfile config/my_imprun_config.json -np`

Snakemake will cycle through rules and you can check to make sure shell commands/outputs look correct. Errors will be thrown if there are issues in config file, etc.

Then run imputation pipeline with the following command, specifying the appropriate cluster config file which can be edited if available resources don't suffice:\
`snakemake -s bin/chip_imp_lewis.snakefile --configfile config/my_imprun_config.json --cluster-config code/cluster/cluster/chip_imp_lewis.cluster.json --cluster "sbatch -p {cluster.p} -o {cluster.o} --account {cluster.account} -t {cluster.t} -c {cluster.c} --mem {cluster.mem}" --jobs 30 -p &> imputation_runs/my_imprun/my_imprun_snakemakerun.log`

### Endpoint files
The pipeline as is outputs the following files in this directory (`imputation_runs/my_imprun/imputed_genotypes`):
* Filter logging `my_imprun_filtering_report.txt`
    + A single file that has filter logging for each filtering step (SNP call rate, individual call rate, and HWE p-value) for each assay. 
* Concatenated hardcall/dosage VCF file (in `imputation_runs/my_imprun/imputed_genotypes/` subdirectory)
    + This is the version of VCF file that comes out of Minimac where genotyeps are represnted like: `1|1:1.999`
    + Here the first half of each SNP for each individual is a phased hard call and after the : is the additive dosage genotype
* Binary PLINK files (*.bed, *.bim, *.fam)


### Running imputation pipeline (Sequence)
Chip-level QC and imputation (shown previously) is linked to this pipeline, so won't need to be run separately \
Config files should be the same as for running 850K imputation \
Upon setting up config file, execute a "dry-run" of Snakemake to ensure that it sees all necessary files, etc. \
`snakemake -s bin/seq_imputation.snakefile --configfile config/my_imprun_config.json -np`

Snakemake will cycle through rules and you can check to make sure shell commands/outputs look correct. Errors will be thrown if there are issues in config file, etc.

Then run imputation pipeline with the following command, specifying the appropriate cluster config file which can be edited if available resources don't suffice:\
`snakemake -s bin/seq_imputation.snakefile --configfile config/my_imprun_config.json --cluster-config code/cluster/cluster/seq_imputation.cluster.json --cluster "sbatch -p {cluster.p} -o {cluster.o} --account {cluster.account} -t {cluster.t} -c {cluster.c} --mem {cluster.mem} --qos {cluster.qos}" --jobs 30 -p &> imputation_runs/my_imprun/my_imprun_snakemakerun.log`\

*NOTE*: The one thing to keep an eye on is the QOS that is needed for imputation. Most chromosomes (in 100K animal datasets) will run in under 2 days. However larger chromosomes may need to be run on the `biolong` QOS. This can all be changed in the cluster config file. 


### Endpoint files
The pipeline as is outputs the following files in this directory (`imputation_runs/my_imprun/seq_imputed`): 
* Single-chromosome hardcall/dosage VCF files (in `imputation_runs/my_imprun/seq_imputed/` subdirectory)
    + This is the version of VCF file that comes out of Minimac where genotyeps are represnted like: `1|1:1.999`
    + Here the first half of each SNP for each individual is a phased hard call and after the : is the additive dosage genotype
    + Minimac info file with imputation accuracy, allele frequencies, etc. 
* The 850K genotypes and endpoint files mentioned above will also exist in the same `imputation_runs/my_imprun/imputed_genotypes` directory. 
