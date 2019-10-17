# UMAG 850K Imputation Pipeline
## Setup and basic pipeline instructions

#### Cloning repository and setting up environments 
Navigate to directory where you want to setup imputation pipeline and clone repository \
`git clone git@github.com:troyrowan/imputation-pipeline.git`

Set up imputation conda environment from .yml file in directory\
The first weird bit of the command deals with some weird dependency issues that I was having before\
Solution came from: https://stackoverflow.com/questions/55661167/how-to-fix-condavalueerror-invalid-environment-name-in-conda-terminal-when-i \
`CONDA_RESTORE_FREE_CHANNEL=1 conda env create -f imputation_env.yml`

Activate environment (This needs to be done every time prior to running Snakemake)\
`source activate Imp3`

### Setting up config files for imputation
Create copy of master config file for your imputation run\
`cp config/master_config.json config/my_imprun_config.json`

Add run-specific information into config file: 
* Only first three items should need changed
* Updates to reference should be reflected in updated master config (which you can get from a `git pull`)
* `sample`
    + List of `.ped` file prefixes in the format: `nsnps.date.nindividuals.manifest`
* `gt_path`
    + Absolute path to where raw ped files reside
* `run_name`
    + Whatever you choose your imputation run name to be
    
Example top of config file is given below:
```python
{
  "sample":["56789.191002.500.A", "123456.191002.250.B"],
  "gt_path":"/CIFS/MUG01_N/deckerje/tnr343/191002_genodump/",
  "run_name":"191002_IMPRUN",
  "ref_dir":"/data/REF_HAPS/850K/"
  "hdref":"777962.180624.8853.A",
  "f250ref":"227234.180624.28421.A",
  "imp_ref_version":"180919_bigref",
  "phasing_ref_version":"180919_bigref"
  "mapdict" : {"139977":"/CIFS/MUG01_N/schnabelr/PLINK_FILES/9913_ARS1.2_139977_GGPHDV3_snpnumber_180613.map", ...},
  "refdict" : {"139977":"/CIFS/MUG01_N/schnabelr/PLINK_FILES/9913_ARS1.2_139977_GGPHDV3_snpnumber_180613.REF", ...},
  "refalleledict" : {"139977":"/CIFS/MUG01_N/schnabelr/PLINK_FILES/9913_ARS1.2_139977_GGPHDV3_snpnumber_180613.REF_ALLELE", ...}
}
```
**NOTE**
As new assays become available, the dictionaries `"mapdict"`, `"refdict"`, and `"refalleledict"` that refer to corresponding `.map`, `.REF`, and `.REF_ALLELE` 
files need to be updated as well. Where the key is `"nsnps"` and value is the absolute file path to each of the three files

### Running imputation pipeline
Upon setting up config file, execute a "dry-run" of Snakemake to ensure that it sees all necessary files, etc.\
`snakemake -s imputation.snakefile --configfile config/my_imprun_config.json -np`

Snakemake will cycle through rules and you can check to make sure shell commands/outputs look correct

Then run imputation pipeline with the following command, specifying an appropriate number of cores and pointing output information to 
a log file:\
`snakemake -s imputation.snakefile --configfile config/my_imprun_config.json --cores 60 -p &> my_imprun_snakemakerun.log`

### Endpoint files
The pipeline as is outputs the following files:
* Filter logging counts
    + A single file that has filter logging for each filtering step (SNP call rate, individual call rate, and HWE p-value) for each assay. 
* Concatenated hardcall/dosage VCF file (in `run_name/imputed_genotypes/` subdirectory)
    + This is the version of VCF file that comes out of Minimac where genotyeps are represnted like: `1|1:1.999`
    + Here the first half of each SNP for each individual is a phased hard call and after the : is the additive dosage genotype
* Hardcall only VCF file (in `run_name/imputed_genotypes/` subdirectory)
    + I use a script to manually extract only hardcall genotypes from above
    + This can be used in other softwares to convert to nearly any other genotype file format
