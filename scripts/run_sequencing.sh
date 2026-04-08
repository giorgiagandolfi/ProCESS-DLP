#!/bin/bash
#SBATCH --partition=EPYC
#SBATCH --nodes=1
#SBATCH --cpus-per-task=5
#SBATCH --mem=40gb
#SBATCH --time=2:00:00
#SBATCH --output=sequencing_%A_%a.out
#SBATCH --error=sequencing_%A_%a.err
#SBATCH --job-name=DLP6_sequencing
#SBATCH --array=1-900%40 


module load singularity
image="/orfeo/cephfs/scratch/cdslab/ggandolfi/process_on_the_fly_muts/process_on_the_fly_v2.sif"
# change with your path to the simulate_tissue.R and simulate_mutation.R scripts
base="/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-DLP/scripts/"
phylo_path="/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-DLP/dlp_simulation_6/phylo_forest_6.sff"
outdir="/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-DLP/sequencing_dlp6"

#


echo "singularity exec --bind /orfeo:/orfeo --no-home $image Rscript $base/simulate_sequencing.R $SLURM_ARRAY_TASK_ID $phylo_path $outdir"
singularity exec --bind /orfeo:/orfeo --no-home $image Rscript $base/simulate_sequencing.R $SLURM_ARRAY_TASK_ID $phylo_path $outdir
