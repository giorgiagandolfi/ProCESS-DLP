#!/bin/bash
#SBATCH --partition=EPYC
#SBATCH --nodes=1
#SBATCH --cpus-per-task=5
#SBATCH --mem=40gb
#SBATCH --time=2:00:00
#SBATCH --output=sequencing_%A_%a.out
#SBATCH --error=sequencing_%A_%a.err
#SBATCH --job-name=DLP_sequencing
#SBATCH --array=1-2 ### change this with the number of your sampled cells 



module load singularity
## change this with the path to the singularity image
image="/orfeo/cephfs/scratch/cdslab/ggandolfi/process_on_the_fly_muts/process_on_the_fly_v2.sif"
# change with your path to the github clones repo
base="/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-DLP/"
## change with the path of your phylogenetic forest
phylo_path="/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-DLP/dlp_simulation_7/phylo_forest_7.sff"
## change with the output folder path
outdir="/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-DLP/sequencing_dlp7"
## coverage for single cells data
sc_coverage="0.2"


echo "singularity exec --bind /orfeo:/orfeo --no-home $image Rscript "$base/scripts/1_sequencing/simulate_sequencing.R" $SLURM_ARRAY_TASK_ID $phylo_path $outdir $sc_coverage"
singularity exec --bind /orfeo:/orfeo --no-home $image Rscript "$base/scripts/1_sequencing/simulate_sequencing.R" $SLURM_ARRAY_TASK_ID $phylo_path $outdir $sc_coverage
