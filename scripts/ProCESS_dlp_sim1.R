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
n_w <- n_h <- 20
ncells <- 40

# adding second mutant
sim$add_mutant(name = "B", growth_rates = 0.5, death_rates = 0.0)
#sim$update_rates(species = "A", rates = c(growth = 0.01, death = 0.1))
sim$mutate_progeny(sim$choose_cell_in("A"), "B")
sim$run_up_to_size(species = 'B', 3000)
#'
# find a tissue rectangle containing 5 cells of type A and B at least
bbox <- sim$search_sample(c("B" = ncells/2,"A"=ncells/2), n_w, n_h)

#'
DLP.sample(sim, bbox$lower_corner, bbox$upper_corner, sample_prefix="DLP1_")
sim$get_samples_info()$tumour_cells %>% table()


forest <- sim$get_sample_forest()
#forest$save("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-DLP/sample_forest_1.sff")
## 
color_palette <- Polychrome::createPalette(nrow(sim$get_samples_info()), seedcolors=c("#000000", "#FFFFFF"))
plot_forest(forest,highlight_sample = T,color_map = color_palette)
my_plot_forest(forest = forest,highlight_sample = T,color_sample=color_palette)+
  theme(legend.position = "none")


#
#setwd("/orfeo/cephfs/scratch/cdslab/ggandolfi/process_on_the_fly_muts/reference")
#m_engine <- MutationEngine(setup_code = "GRCh38",tumour_type = "COADREAD", context_sampling = 20)
#
#
#mu_SNV = 1e-8
#mu_CNA = 0
#mu_INDELs = 1e-9
#
#CNA_Clone2 = ProCESS::CNA(type = "D", "5",
#                          from = 107707518, len = 2e7,allele = 0)
#
### Drivers for the tumors
#m_engine$add_mutant(mutant_name = "A",
#                    passenger_rates = c(SNV = mu_SNV, CNA = 0,indel=mu_INDELs),drivers = list(list("APC R1450*", allele = 1)))
#m_engine$add_mutant(mutant_name = "B",passenger_rates = c(SNV = mu_SNV, CNA = mu_CNA,indel=mu_INDELs),drivers = list(CNA_Clone2))
#m_engine$add_exposure(time = 0,coefficients = c(SBS1 = 0.15,SBS5 = 0.40,
#                                                SBS18 = 0.15,SBS17b = 0.20,ID1 = 0.40,ID2 = 0.40,ID18=0.2,SBS88 = 0.10))
#phylo_forest <- m_engine$place_mutations(forest, num_of_preneoplatic_SNVs=800, num_of_preneoplatic_indels=200)
#phylo_forest$save("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-DLP/phylo_forest_1.sff")
#phylo_forest <- load_phylogenetic_forest("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-DLP/phylo_forest_1.sff")
phylo_forest <- load_phylogenetic_forest("/orfeo/cephfs/scratch/cdslab/shared/process_dlp_data/phylo_forest_1.sff")
cov <- DLP.coverage(phylo_forest, 20, sample_prefix = "DLP1_")
cov

no_error_seq <- ErrorlessIlluminaSequencer()

chromosomes <- phylo_forest$get_absolute_chromosome_positions()$chr

results <- parallel::mclapply(
  chromosomes,
  function(chr) {
    simulate_seq(
      phylo_forest,
      coverage = cov,
      write_SAM = F,#TRUE
      with_normal_sample = FALSE,chromosomes = c(chr),
      read_size = 150,
      sequencer = no_error_seq,
      insert_size_mean = 350,
      insert_size_stddev = 10,
      missed_SID_statistics=F, germline_statistics=F,
      wide_format=FALSE,
      #output_dir = "/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-DLP/sequencing_dlp1",
      update_SAM = FALSE #TRUE
    )
  },
  mc.cores = 4
)

saveRDS(object = results,file = "/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-DLP/seq_res_1.rds")

