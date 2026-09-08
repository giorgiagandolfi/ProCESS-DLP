#!/bin/bash
#SBATCH --partition=EPYC
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=40gb
#SBATCH --time=2:00:00
#SBATCH --output=SPN01_DLP.out
#SBATCH --error=SPN01_DLP.err
#SBATCH --job-name=SPN01_DLP

module load singularity
image="/orfeo/cephfs/scratch/cdslab/ggandolfi/process_on_the_fly_muts/process_on_the_fly_v1.sif"
# change with your path to the simulate_tissue.R and simulate_mutation.R scripts
base="/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-DLP/scripts/SCOUT/SPN01"
#

singularity exec --bind /orfeo:/orfeo --no-home $image Rscript $base/simulate_mutation.R