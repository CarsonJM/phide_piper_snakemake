# Phide Piper instructions

## 1. Make and activate snakemake conda environment using following command

`mamba create -n phide_piper`

`mamba activate phide_piper`

`mamba install -c conda-forge -c bioconda snakemake snakefmt snakedeploy git biopython plotnine matplotlib -y`

## 2. Make a new directory where you want the analysis to take place

`mkdir <insert directory name here>`

## 3. Change to the location of the cloned directory

`cd <insert directory name here>`

## 4. Run the following command to deploy the phide_piper workflow in the specified directory

`snakedeploy deploy-workflow https://github.com/CarsonJM/phide_piper_dev.git . --branch master`

*You should see a 'config' and 'workflow' directories now*

## 5. Modify the config/config.yaml file so that:

- the samples_df is the location of the samples.tsv file
- the resource directory is where you want to store large downloaded databases
- the results directory is where you want to store your results
- other workflow parameters are set as desired

## 6. Modify the sample.tsv file so that it matches the desired sample names and paths

## 7. Run a dry run of the workflow to verify everything is set up correctly

`snakemake --configfile <path to config file> --dry-run`

## 8. To run using SLURM, use the following command:

`snakemake --profile config/simple/ --configfile <path to config file>`

## 9. After a successful run, use following to generate a report from the workflow:

`snakemake --profile config/simple/ --configfile <path to config file> --report <path to desired report location>`