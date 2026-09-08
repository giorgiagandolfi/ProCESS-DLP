library(dplyr)
library(ProCESS)
library(parallel)

args <- commandArgs(trailingOnly = TRUE)
cell_id <- args[1]
step=1
cat(paste0(step, ". Mergin SAM file by chromosomes..."))
step <- step + 1
merge_sams <- function(output_local_dir, BAM_file,
                       chromosomes, num_of_cores) {
  
  SAM_files <- ""
  for (i in 1:length(chromosomes)) {
    chr_SAM_file <- file.path(output_local_dir,
                              paste0(chromosomes[i], ".sam"))
    SAM_files <- paste(SAM_files, chr_SAM_file)
  }
  
  cmd <- paste("samtools merge -fc -@", num_of_cores,
               "-o", BAM_file, SAM_files)
  print(cmd)
  invisible(system(cmd, intern = TRUE))
}

merged_sam_dir <- paste0("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-DLP/SCOUT/SPN01/process/",cell_id,"/MERGED_SAM")
dir.create(merged_sam_dir)
output_local_dir <- paste0("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-DLP/SCOUT/SPN01/process/",cell_id,"/SAM")
chromosomes <-paste0("chr_",c(seq_along(1:22),"X","Y"))

merge_sams(output_local_dir = output_local_dir,BAM_file = file.path(merged_sam_dir,paste0(cell_id,".bam")),
           num_of_cores = 8,chromosomes = chromosomes)
# ##################################################################################
# ##################################################################################
#

cat(paste0(step, ". Generating FASTQ files for each sample..."))
step <- step + 1

BAM_files <- list.files(merged_sam_dir, pattern = "\\.bam$")
# 
generate_fastq <- function(orig_file, fastq_dir) {
  base_orig_file <- tools::file_path_sans_ext(basename(orig_file))
  
  file_prefix <- file.path(fastq_dir, base_orig_file)
  R1 <- paste0(file_prefix, ".R1.fastq.gz")
  R2 <- paste0(file_prefix, ".R2.fastq.gz")
  unpaired <- paste0(file_prefix, ".unpaired.fastq.gz")
  singleton <- paste0(file_prefix, ".singleton.fastq.gz")
  
  prefix <- strsplit(file_prefix,'/') %>% unlist()
  prefix <- prefix[length(prefix)]
  cmd <- paste("samtools fastq -@ 8 -c 9 -N -1", R1, "-2", R2, "-0", unpaired,
               "-s", singleton, orig_file)
  print(cmd)
  invisible(system(cmd, intern = TRUE))
}
# 

fastq_dir <- paste0("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-DLP/SCOUT/SPN01/process/FASTQ")
dir.create(fastq_dir)

for (c in BAM_files){
  curr_BAM_file <- file.path(merged_sam_dir, c)
  generate_fastq(curr_BAM_file, fastq_dir)
}
