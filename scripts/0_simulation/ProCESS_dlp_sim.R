rm(list=ls())
library(ProCESS)
library(dplyr)
library(ggplot2)
source("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-DLP/scripts/utils/DLP.R")
source("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-DLP/scripts/utils/utils_plot.R")
set.seed(12345)
sim <- TissueSimulation(name = "DLP")
sim$add_mutant(name = "A", growth_rates = 0.1, death_rates = 0)
sim$place_cell("A", 500, 500)
sim$run_up_to_size(species = 'A', 5000)

# sampling tissue
n_w <- n_h <- 30
ncells <- 50
bbox <- sim$search_sample(c("A" = ncells), n_w, n_h)
# collect a DLP+ sample from the tissue simulation `sim`
DLP.sample(sim, bbox$lower_corner, bbox$upper_corner, sample_prefix="test")
sim$get_samples_info()$tumour_cells %>% table()

# adding second mutant
sim$add_mutant(name = "B", growth_rates = 0.5, death_rates = 0.0)
#sim$update_rates(species = "A", rates = c(growth = 0.01, death = 0.1))
sim$mutate_progeny(sim$choose_cell_in("A"), "B")
sim$run_up_to_size(species = 'B', 3000)
#'
# find a tissue rectangle containing 5 cells of type A and B at least
bbox <- sim$search_sample(c("B" = ncells), n_w, n_h)

#'
DLP.sample(sim, bbox$lower_corner, bbox$upper_corner, sample_prefix="test")
sim$get_samples_info()$tumour_cells %>% table()


forest <- sim$get_sample_forest()
forest$save("/orfeo/cephfs/scratch/cdslab/ggandolfi/DLP/sample_forest_500cells.sff")
# 
# color_palette <- Polychrome::createPalette(nrow(sim$get_samples_info()), seedcolors=c("#000000", "#FFFFFF"))
# plot_forest(forest,highlight_sample = T,color_map = color_palette)
# my_plot_forest(forest = forest,highlight_sample = T,color_sample=color_palette)+
#   theme(legend.position = "none")
# 


setwd("/orfeo/cephfs/scratch/cdslab/shared/SCOUT/")
m_engine <- MutationEngine(setup_code = "GRCh38",tumour_type = "COADREAD", context_sampling = 20)
# COSMIC_account = list("email"="giorgia.gandolfi@phd.units.it","password"="2*db!XQ4sgQ!dbg"))


mu_SNV = 1e-8
mu_CNA = 0 #2e-10
mu_INDELs = 1e-9
##112707518-112846239 
CNA_Clone2 = ProCESS::CNA(type = "D", "5",
                          from = 107707518, len = 2e7,allele = 0)

## Drivers for the tumors
m_engine$add_mutant(mutant_name = "A",
                    passenger_rates = c(SNV = mu_SNV, CNA = 0,indel=mu_INDELs),drivers = list(list("APC R1450*", allele = 1)))
m_engine$add_mutant(mutant_name = "B",passenger_rates = c(SNV = mu_SNV, CNA = mu_CNA,indel=mu_INDELs),drivers = list(CNA_Clone2))
m_engine$add_exposure(time = 0,coefficients = c(SBS1 = 0.15,SBS5 = 0.40,
                                                SBS18 = 0.15,SBS17b = 0.20,ID1 = 0.40,ID2 = 0.40,ID18=0.2,SBS88 = 0.10))
phylo_forest <- m_engine$place_mutations(forest, num_of_preneoplatic_SNVs=800, num_of_preneoplatic_indels=200)
phylo_forest$save("/orfeo/cephfs/scratch/cdslab/ggandolfi/DLP/phylo_forest_500cells.sff")
cov <- DLP.coverage(phylo_forest, 50, sample_prefix = "test")
cov


#basic_seq <- BasicIlluminaSequencer(1e-3) ## only for testing purpose
no_error_seq <- ErrorlessIlluminaSequencer()

chromosomes <- phylo_forest$get_absolute_chromosome_positions()$chr
# chromosomes <- c("5","22")
#ref_path <- file.path(output_local_dir, "reference.fasta")



results <- parallel::mclapply(
  chromosomes,
  function(chr) {
    simulate_seq(
      phylo_forest,
      coverage = cov,
      write_SAM = FALSE,
      with_normal_sample = FALSE,chromosomes = c(chr),
      read_size = 150,
      sequencer = no_error_seq,
      insert_size_mean = 350,
      insert_size_stddev = 10,
      missed_SID_statistics=F, germline_statistics=F,
      wide_format=FALSE,
      # include_non_sequenced_mutations = TRUE,
      output_dir = "/orfeo/cephfs/scratch/cdslab/ggandolfi/DLP/sequencing_test_500cells",
      update_SAM = TRUE
    )
  },
  mc.cores = 4
)

saveRDS(object = results,file = "/orfeo/cephfs/scratch/cdslab/ggandolfi/DLP/seq_res_500cells.rds")
# t = readRDS("seq_res_new.rds")
# t[[1]]$mutations %>% 
#   select(c("chr","chr_pos","ref","alt"),contains("test_446.452")) %>% 
#   arrange(desc(test_446.452.occurrences)) %>% head()
# 
# phylo_forest <- load_phylogenetic_forest("/orfeo/cephfs/scratch/cdslab/ggandolfi/DLP/phylo_forest.sff")
