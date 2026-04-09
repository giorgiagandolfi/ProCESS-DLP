#!/bin/bash
#SBATCH --partition=EPYC
#SBATCH --nodes=1
#SBATCH --cpus-per-task=9
#SBATCH --mem=100gb
#SBATCH --time=4:00:00
#SBATCH --output=Preprocess_DLP.out
#SBATCH --error=Preprocess_DLP.err
#SBATCH --job-name=Preprocess_DLP
#SBATCH --array=1-2%2 ### change this with the number of your sampled cells 


module load singularity
## change this with the path to the singularity image
image="/orfeo/cephfs/scratch/cdslab/ggandolfi/process_on_the_fly_muts/process_on_the_fly_v2.sif"
# change with your path to the github clones repo
base="/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-DLP/"
phylo_path="/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-DLP/dlp_simulation_7/phylo_forest_7.sff"
outdir="/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-DLP/sequencing_dlp7"



singularity exec --bind /orfeo:/orfeo --no-home $image Rscript "$base/scripts/2_preprocessing/preprocess_SAM.R" $phylo_path $outdir $SLURM_ARRAY_TASK_ID
