#!/bin/bash
#SBATCH --partition=GENOA
#SBATCH --nodes=1
#SBATCH --cpus-per-task=5
#SBATCH --mem=40gb
#SBATCH --time=2:00:00
#SBATCH --output=sequencing_%A_%a.out
#SBATCH --error=sequencing_%A_%a.err
#SBATCH --job-name=SPN01_DLP_sequencing
#SBATCH --array=1-299%20


module load singularity
image="/orfeo/cephfs/scratch/cdslab/ggandolfi/process_on_the_fly_muts/process_on_the_fly_v1.sif"
# change with your path to the simulate_tissue.R and simulate_mutation.R scripts
base="/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-DLP/scripts/SCOUT/SPN01"
#


echo "singularity exec --bind /orfeo:/orfeo --no-home $image Rscript $base/simulate_sequencing.R $SLURM_ARRAY_TASK_ID"
singularity exec --bind /orfeo:/orfeo --no-home $image Rscript $base/simulate_sequencing.R $SLURM_ARRAY_TASK_ID
