#!/bin/bash
#SBATCH --partition=EPYC
#SBATCH --nodes=1
#SBATCH --cpus-per-task=5
#SBATCH --mem=40gb
#SBATCH --time=2:00:00
#SBATCH --output=process_sim_%J.out
#SBATCH --error=process_sim_%J.err
#SBATCH --job-name=process_sim_dlp



module load singularity
image="/orfeo/cephfs/scratch/cdslab/ggandolfi/process_on_the_fly_muts/process_on_the_fly_v2.sif"
# change with your path to the simulate_tissue.R and simulate_mutation.R scripts
base="/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-DLP/scripts/0_simultation"
outdir="/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-DLP/sequencing_dlp7"

#


echo "singularity exec --bind /orfeo:/orfeo --no-home $image Rscript $base/ProCESS_dlp_sim7.R"
singularity exec --bind /orfeo:/orfeo --no-home $image Rscript $base/ProCESS_dlp_sim7.R
