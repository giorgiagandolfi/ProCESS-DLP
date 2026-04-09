rm(list=ls())
library(ProCESS)
library(dplyr)
library(ggplot2)
source("../utils/DLP.R")
source("../utils/utils_plot.R")
set.seed(12345)
args <- commandArgs(trailingOnly = TRUE)
outdir <- args[1]


dir.create(outdir)
setwd(outdir)
sim_id <- "DPL7"
sim <- TissueSimulation(name = sim_id,save_snapshots = F)
sim$add_mutant(name = "A", growth_rates = 0.1, death_rates = 0)
sim$place_cell("A", 500, 500)
sim$run_up_to_size(species = 'A', 3000)

# sampling tissue
n_w <- n_h <- 30
ncells <- 50

# adding second mutant
sim$add_mutant(name = "B", growth_rates = 0.5, death_rates = 0.0)
#sim$update_rates(species = "A", rates = c(growth = 0.01, death = 0.1))
sim$mutate_progeny(sim$choose_cell_in("A"), "B")
sim$run_up_to_size(species = 'B', 2000)

# adding third mutant
sim$add_mutant(name = "C", growth_rates = 0.8, death_rates = 0.0)
sim$mutate_progeny(sim$choose_cell_in("A"), "C")
sim$run_up_to_size(species = 'C', 2000)
# find a tissue rectangle containing 5 cells of type A and B at least
bbox1 <- sim$search_sample(c("B" = ncells/3,"A"=ncells/3,"C"=ncells/3), n_w, n_h)
DLP.sample(sim, bbox1$lower_corner, bbox1$upper_corner, sample_prefix=paste0(sim_id,"_A"))

#bbox1 <- sim$search_sample(c("A"=ncells), n_w, n_h)
#DLP.sample(sim, bbox1$lower_corner, bbox1$upper_corner, sample_prefix="DLP2_A")
sim$get_samples_info()$tumour_cells %>% table()

forest <- sim$get_sample_forest()
forest$save(paste0(outdir,"sample_forest.sff"))
## Only for plotting purpose
# color_palette <- c("A"="goldenrod","C"="forestgreen","B"="purple" )
# my_plot_forest(forest = forest,highlight_sample = T,color_map = color_palette,horizontal = F,color_sample = F)+theme(legend.position = "none")

# setwd("/orfeo/cephfs/scratch/cdslab/ggandolfi/process_on_the_fly_muts/reference")
m_engine <- MutationEngine(setup_code = "GRCh38",tumour_type = "COADREAD", context_sampling = 20)


mu_SNV = 1e-8
# mu_CNA = 5e-10
mu_INDELs = 1e-9

CNA_Clone2 = ProCESS::CNA(type = "D", "5",
                          from = 107707518, len = 2e7)
CNA_Clone1 = ProCESS::CNA(type = "A", "12",
                          from = 10000000, len = 4e7)
CNA_Clone3 = ProCESS::CNA(type = "A", "1",
                          from = 10000000, len = 2e7)
## Drivers for the tumors
m_engine$add_mutant(mutant_name = "A",
                    passenger_rates = c(SNV = mu_SNV, CNA = 1e-12,indel=mu_INDELs),drivers = list("APC R1450*",CNA_Clone1))
m_engine$add_mutant(mutant_name = "B",passenger_rates = c(SNV = mu_SNV, CNA = 1e-11,indel=mu_INDELs),
                    drivers = list(CNA_Clone2))
m_engine$add_mutant(mutant_name = "C",passenger_rates = c(SNV = mu_SNV, CNA = 5e-11,indel=mu_INDELs),drivers = list("KRAS G12D",CNA_Clone3))

m_engine$add_exposure(time = 0,coefficients = c(SBS1 = 0.15,SBS5 = 0.40,
                                                SBS18 = 0.15,SBS17b = 0.20,ID1 = 0.40,ID2 = 0.40,ID18=0.2,SBS88 = 0.10))
phylo_forest <- m_engine$place_mutations(forest, num_of_preneoplatic_SNVs=800, num_of_preneoplatic_indels=200)
phylo_forest$save(paste0(outdir,"phylo_forest.sff"))