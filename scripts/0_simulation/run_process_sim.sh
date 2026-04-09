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
## change this with the path to the singularity image
image="/orfeo/cephfs/scratch/cdslab/ggandolfi/process_on_the_fly_muts/process_on_the_fly_v2.sif"
# change with your path to the github clones repo
base="/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-DLP/"
outdir="/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-DLP/dlp_simulation_7" ## to change

#


echo "singularity exec --bind /orfeo:/orfeo --no-home $image Rscript "$base/scripts/0_simulation/ProCESS_dlp_sim7.R" $outdir"
singularity exec --bind /orfeo:/orfeo --no-home $image Rscript "$base/scripts/0_simulation/ProCESS_dlp_sim7.R" $outdir
