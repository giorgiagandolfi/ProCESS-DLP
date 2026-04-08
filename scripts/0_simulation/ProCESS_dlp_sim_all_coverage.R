rm(list=ls())
library(ProCESS)
library(dplyr)
library(ggplot2)
source("/orfeo/cephfs/scratch/cdslab/ggandolfi/DLP/DLP.R")
source("/orfeo/cephfs/scratch/cdslab/ggandolfi/DLP/utils_plot.R")
set.seed(12345)
sim <- TissueSimulation()
sim$add_mutant(name = "A", growth_rates = 0.1, death_rates = 0)
sim$place_cell("A", 500, 500)
sim$run_up_to_size(species = 'A', 10000)

# sampling tissue
n_w <- n_h <- 10
ncells <- 0.8 * n_w * n_h
bbox <- sim$search_sample(c("A" = ncells), n_w, n_h)
# collect a DLP+ sample from the tissue simulation `sim`
sim$sample_cells("S1", bbox$lower_corner, bbox$upper_corner)
plot_tissue(sim)
# adding second mutant
sim$add_mutant(name = "B", growth_rates = 0.5, death_rates = 0.0)
#sim$update_rates(species = "A", rates = c(growth = 0.01, death = 0.1))
sim$mutate_progeny(sim$choose_cell_in("A"), "B")
sim$run_up_to_size(species = 'B', 8000)
#'
# find a tissue rectangle containing 5 cells of type A and B at least
bbox <- sim$search_sample(c("B" = ncells), n_w, n_h)
sim$sample_cells("S2", bbox$lower_corner, bbox$upper_corner)

#'



forest <- sim$get_sample_forest()
color_palette <- c("S1"="plum2","S2"="plum3")
my_plot_forest(forest = forest,highlight_sample = T,color_sample=color_palette)+
  theme(legend.position = "none")





setwd("/orfeo/cephfs/scratch/cdslab/shared/SCOUT/")
m_engine <- MutationEngine(setup_code = "GRCh38",tumour_type = "COADREAD", context_sampling = 20)
# COSMIC_account = list("email"="giorgia.gandolfi@phd.units.it","password"="2*db!XQ4sgQ!dbg"))


mu_SNV = 1e-8
mu_CNA = 2e-10
mu_INDELs = 1e-9
##112707518-112846239 
CNA_Clone2 = ProCESS::CNA(type = "D", "5",
                          chr_pos = 107707518, len = 2e7,allele = 0)

## Drivers for the tumors
m_engine$add_mutant(mutant_name = "A",
                    passenger_rates = c(SNV = mu_SNV, CNA = 0,indel=mu_INDELs),drivers = list(list("APC R1450*", allele = 1)))
m_engine$add_mutant(mutant_name = "B",passenger_rates = c(SNV = mu_SNV, CNA = mu_CNA,indel=mu_INDELs),drivers = list(CNA_Clone2))
m_engine$add_exposure(time = 0,coefficients = c(SBS1 = 0.15,SBS5 = 0.40,
                                                SBS18 = 0.15,SBS17b = 0.20,ID1 = 0.40,ID2 = 0.40,ID18=0.2,SBS88 = 0.10))
phylo_forest <- m_engine$place_mutations(forest, num_of_preneoplatic_SNVs=800, num_of_preneoplatic_indels=200)

seq_results <- simulate_seq(phylo_forest, coverage = 10, write_SAM = T,
                            with_normal_sample = FALSE,chromosomes = c("22"),output_dir = "/orfeo/cephfs/scratch/cdslab/ggandolfi/DLP/sequencing_test",update_SAM = T)

seq_results$mutations %>% 
  filter(classes!="germinal") %>% 
  filter(S1.occurrences!=0) %>% 
  ggplot(aes(x=S1.coverage))+geom_histogram()






