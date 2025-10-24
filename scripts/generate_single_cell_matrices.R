library(ProCESS)
library(dplyr)
setwd("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-DLP/dlp_simulation_1")
phylo_forest <- load_phylogenetic_forest("phylo_forest_1.sff")
sample_forest <- load_sample_forest("sample_forest_1.sff")

chromosomes <- phylo_forest$get_absolute_chromosome_positions()$chr
seq_res <- readRDS("seq_res_2.rds") ### seq res
somatic_dp <- lapply(seq_along(chromosomes), function(x){
  seq_res[[x]]$mutations %>% 
    filter(classes!="germinal") %>% 
    mutate(mutationID=paste0("chr",chr,":",from,":",ref,":",alt)) %>% 
    select(mutationID,sample,DP)
}) %>% bind_rows()

somatic_nv <- lapply(seq_along(chromosomes), function(x){
  seq_res[[x]]$mutations %>% 
    filter(classes!="germinal") %>% 
    mutate(mutationID=paste0("chr",chr,":",from,":",ref,":",alt)) %>% 
    select(mutationID,sample,NV)
}) %>% bind_rows()


saveRDS(object = somatic_nv,file = "/orfeo/cephfs/scratch/cdslab/shared/process_dlp_data/somatic_nv.rds")
saveRDS(object = somatic_dp,file = "/orfeo/cephfs/scratch/cdslab/shared/process_dlp_data/somatic_dp.rds")
