#!/bin/bash
#SBATCH --partition=EPYC
#SBATCH --nodes=1
#SBATCH --cpus-per-task=5
#SBATCH --mem=40gb
#SBATCH --time=2:00:00
#SBATCH --output=sequencing_normal_%A_%a.out
#SBATCH --error=sequencing_normal_%A_%a.err
#SBATCH --job-name=BULK_sequencing_normal
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
### Add here the chunk coverage
### SUGGESTION: the chunk coverage should be around 5x
### change the size of the job array: total coverage (eg 100x)/5x = 20
coverage="1"
ref_path="/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-DLP/dlp_simulation_7/"

echo "singularity exec --bind /orfeo:/orfeo --no-home $image Rscript "$base/scripts/1_sequencing/simulate_sequencing_normal.R" $SLURM_ARRAY_TASK_ID $phylo_path $outdir $coverage $ref_path"
singularity exec --bind /orfeo:/orfeo --no-home $image Rscript "$base/scripts/1_sequencing/simulate_sequencing_normal.R" $SLURM_ARRAY_TASK_ID $phylo_path $outdir $coverage $ref_path
