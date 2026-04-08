#!/bin/bash
#SBATCH --partition=EPYC
#SBATCH --nodes=1
#SBATCH --cpus-per-task=9
#SBATCH --mem=100gb
#SBATCH --time=4:00:00
#SBATCH --output=Preprocess_DLP.out
#SBATCH --error=Preprocess_DLP.err
#SBATCH --job-name=Preprocess_DLP

module load singularity
image="/orfeo/cephfs/scratch/cdslab/ggandolfi/process_on_the_fly_muts/process_on_the_fly_v2.sif"
# change with your path to the simulate_tissue.R and simulate_mutation.R scripts
base="/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-DLP/scripts/2_preprocessing"
phylo_path="/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-DLP/dlp_simulation_6/phylo_forest_6.sff"
outdir="/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-DLP/sequencing_dlp6"



singularity exec --bind /orfeo:/orfeo --no-home $image Rscript $base/preprocess_SAM.R $phylo_path $outdir
