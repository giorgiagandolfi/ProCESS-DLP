#!/bin/bash
#SBATCH --partition=GENOA
#SBATCH --nodes=1
#SBATCH --cpus-per-task=5
#SBATCH --mem=40gb
#SBATCH --time=2:00:00
#SBATCH --output=preprocessing_%A_%a.out
#SBATCH --error=preprocessing_%A_%a.err
#SBATCH --job-name=SPN01_DLP_preprocessing
#SBATCH --array=1-1000%20


module load singularity
image="/orfeo/cephfs/scratch/cdslab/ggandolfi/process_on_the_fly_muts/process_on_the_fly_v1.sif"
# change with your path to the simulate_tissue.R and simulate_mutation.R scripts
base="/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-DLP/scripts/SCOUT/SPN01"
#

cell_ids=($(ls /orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-DLP/SCOUT/SPN01/process | grep seq))
echo "singularity exec --bind /orfeo:/orfeo --no-home $image Rscript $base/preprocess_SAM.R ${cell_ids[$SLURM_ARRAY_TASK_ID-1]}"
singularity exec --bind /orfeo:/orfeo --no-home $image Rscript $base/preprocess_SAM.R ${cell_ids[$SLURM_ARRAY_TASK_ID-1]}
