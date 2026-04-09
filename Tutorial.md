# Tutorial ProCESS-DLP 

This is an easy to use tutorial to reproduce single-cell DNA sequencing data. 

## 0. Set-up
In order to perform a DLP-like ProCESS simulation, you need to install a specific ProCESS version from a specific branch.

```r
devtools::install_github("caravagnalab/ProCESS", ref="on-the-fly_mutations")
```

Or either using the docker image with singularity:

```bash
singularity pull docker://giorgiagandolfi97/process_on_the_fly:v2
```

## 1. Simulate tissue and sampling

Run the Rscript in `scripts/0_simulation/ProCESS_dlp_sim7.R`. This scripts will generate a three clone CRC tumour. Each clone is
characterized by a specific growht rate and a fixed SNV/INDEL rate. CN rate is changing across clones. A total of 900 cells are sampled.
You can run the script by submitting the `scripts/0_simulation/run_process_sim.sh` to the SLURM jobs scheduler.

## 2. Simulate sequencing

Run the `sh` script that you find in `scripts/1_sequencing/run_sequencing.sh`.This script will perform sequencing of 
sampled cells in parallel as a job array. 

## 3. Preprocess SAM files and generate FASTQs

Run the `sh` script that you find in `scripts/2_preprocessing/run_preprocessing.sh`. This script will:

1. merge all SAM files for each chromosome
2. create FASTQs from SAM files

The output folder will contain paired end FASTQs files for each sampled cell.
