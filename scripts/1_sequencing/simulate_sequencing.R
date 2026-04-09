rm(list=ls())
library(ProCESS)
library(dplyr)
library(ggplot2)
source("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-DLP/scripts/utils/DLP.R")
source("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-DLP/scripts/utils/utils_plot.R")
set.seed(12345)

setwd("/orfeo/cephfs/scratch/cdslab/ggandolfi/process_on_the_fly_muts/reference")

args <- commandArgs(trailingOnly = TRUE)
cell_idx <- as.numeric(args[1])
phylo_path <- args[2]
outdir <- args[3]
sc_coverage <- as.numeric(args[4])
print(cell_idx)

phylo_forest <- load_phylogenetic_forest(phylo_path)
cell_names <- phylo_forest$get_samples_info()$name
cell_name <- cell_names[cell_idx]
print(cell_name)
subforest <-  phylo_forest$get_subforest_for(cell_name)



# phylo_forest <- load_phylogenetic_forest(paste0(outdir,"phylo_forest.sff"))
# cov <- DLP.coverage(phylo_forest, 100, sample_prefix = "SPN01_")


no_error_seq <- ErrorlessIlluminaSequencer()

dir.create(path = paste0(outdir,"/sequencing_",cell_name),recursive = T)
			 
chromosomes <- subforest$get_absolute_chromosome_positions()$chr

start_time <- Sys.time()
results <- parallel::mclapply(
  chromosomes,
  function(chr) {
    simulate_seq(
      subforest,
      coverage = sc_coverage,
      write_SAM = TRUE,
      with_normal_sample = FALSE,chromosomes = c(chr),
      read_size = 150,
      sequencer = no_error_seq,
      insert_size_mean = 350,
      insert_size_stddev = 10,
      missed_SID_statistics=T, germline_statistics=F,
      wide_format=FALSE,
      output_dir = paste0(outdir,"/sequencing_",cell_name),
      update_SAM = TRUE
    )
  },
  mc.cores = 4
)
end_time <- Sys.time()
elapsed_time <- end_time - start_time
elapsed_time <- as.numeric(elapsed_time, units = "mins")
print(elapsed_time)

saveRDS(object = results,file = paste0(outdir,"/sequencing_",cell_name,"/seq_res.rds"))
