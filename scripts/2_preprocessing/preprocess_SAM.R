library(dplyr)
library(ProCESS)
library(parallel)
args <- commandArgs(trailingOnly = TRUE)
phylo_path <- args[1]
outdir <- args[2]
cell_idx <- as.numeric(args[3])

phylo_forest <- load_phylogenetic_forest(phylo_path)
cell_names <- phylo_forest$get_samples_info()$name
cell_name <- cell_names[cell_idx]

outdir_cell <- paste0(outdir,"/sequencing_",cell_name)
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
# merged_sam_dir <- paste0(outdir_cell,"/MERGED_SAM")
# dir.create(merged_sam_dir,recursive = T)

chromosomes <-paste0("chr_",c(seq_along(1:22),"X","Y"))

merge_sams(output_local_dir = outdir_cell,BAM_file = file.path(outdir_cell,paste0(cell_name,".bam")),
                     num_of_cores = 2,chromosomes = chromosomes)
# ##################################################################################
# ##################################################################################
#
# cat(paste0(step, ". Splitting BAM file by sample..."))
# step <- step + 1
# 
# split_bam_by_samples <- function(output_local_dir, BAM_file, remove_local_bam, num_of_cores) {
#   name <- strsplit(BAM_file, '/') %>% unlist()
#   name <- name[length(name)]
#   cmd <- paste0("samtools split -f ", file.path(output_local_dir, "%*_%!.bam "), BAM_file, " -@ ", num_of_cores)
#   print(cmd)
#   invisible(system(cmd, intern = TRUE))
# 
#   if (remove_local_bam) {
#     file.remove(BAM_file)
#   }
# }
# splitted_bam_dir <- paste0(outdir_cell,"/SPLITTED_BAM")
# dir.create(splitted_bam_dir)
# 
# 
# split_bam_by_samples(output_local_dir=splitted_bam_dir,
#                                BAM_file=file.path(merged_sam_dir,"merged.bam"),
#                                remove_local_bam=FALSE, num_of_cores=8)

# 
# ##################################################################################
# ##################################################################################
# 
cat(paste0(step, ". Generating FASTQ files for each sample..."))
step <- step + 1

BAM_files <- list.files(outdir_cell, pattern = "\\.bam$")
# 
generate_fastq <- function(orig_file, fastq_dir,remove_local_bam) {
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
  if (remove_local_bam) {
    file.remove(orig_file)
  }
}
# 

# fastq_dir <- paste0(outdir_cell,"/FASTQ")
# dir.create(fastq_dir,recursive = T)

for (c in BAM_files){
  curr_BAM_file <- file.path(outdir_cell,c)
  generate_fastq(curr_BAM_file, outdir_cell,remove_local_bam=T)
}
