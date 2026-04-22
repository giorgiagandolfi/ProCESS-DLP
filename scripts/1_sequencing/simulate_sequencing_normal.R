rm(list=ls())
library(ProCESS)
library(dplyr)
library(ggplot2)
source("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-DLP/scripts/utils/DLP.R")
source("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-DLP/scripts/utils/utils_plot.R")
set.seed(12345)

#setwd("/orfeo/cephfs/scratch/cdslab/ggandolfi/process_on_the_fly_muts/reference")

args <- commandArgs(trailingOnly = TRUE)
chunk_id <- args[1]
chunk_name <- paste0("chunk_0",chunk_id)
phylo_path <- args[2]
outdir <- args[3]
bulk_coverage <- as.numeric(args[4])
ref_path <- args[5]

phylo_forest <- load_phylogenetic_forest(phylo_path)
setwd(ref_path)


no_error_seq <- ErrorlessIlluminaSequencer()

dir.create(path = paste0(outdir,"/normal_sample/sequencing_",chunk_name),recursive = T)

chromosomes <- phylo_forest$get_absolute_chromosome_positions()$chr

start_time <- Sys.time()
results <- parallel::mclapply(
  chromosomes,
  function(chr) {
    simulate_normal_seq(
      phylo_forest,
      coverage = bulk_coverage,
      write_SAM = TRUE,
      chromosomes = c(chr),
      read_size = 150,
      sequencer = no_error_seq,
      insert_size_mean = 350,
      insert_size_stddev = 10,
      missed_SID_statistics=T,
      wide_format=FALSE,
      output_dir = paste0(outdir,"/normal_sample/sequencing_",chunk_name),
      update_SAM = TRUE
    )
  },
  mc.cores = 4
)
end_time <- Sys.time()
elapsed_time <- end_time - start_time
elapsed_time <- as.numeric(elapsed_time, units = "mins")
print(elapsed_time)

saveRDS(object = results,file = paste0(outdir,"/normal_sample/sequencing_",chunk_name,"/seq_res.rds"))