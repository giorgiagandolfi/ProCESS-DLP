rm(list=ls())
library(ProCESS)
library(dplyr)
library(ggplot2)
source("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-DLP/scripts/utils/DLP.R")

set.seed(12345)
sim <- TissueSimulation(name = "DLP",save_snapshots = F)
sim$add_mutant(name = "A", list(duplication = 0.1, death = 0))
sim$place_cell("A", 500, 500)
sim$run_up_to_size(species = 'A', 3000)

# sampling tissue
n_w <- n_h <- 20
ncells <- 40

# adding second mutant
sim$add_mutant(name = "B", list(duplication =  0.5, death = 0.0))
sim$mutate_progeny(sim$choose_cell_in("A"), "B")
sim$run_up_to_size(species = 'B', 2000)

# adding third mutant
sim$add_mutant(name = "C", list(duplication = 0.8, death= 0.0))
sim$set_rates(list("A" = list(duplication = 0.01, death = .1)))
sim$set_rates(list("B" = list(duplication = 0.3, death = .1)))

sim$mutate_progeny(sim$choose_cell_in("A"), "C")
sim$run_up_to_size(species = 'C', 2000)
# find a tissue rectangle containing 5 cells of type A and B at least
bbox1 <- sim$search_sample(c("B" = ncells/2,"A"=ncells/2), n_w, n_h)
bbox2 <- sim$search_sample(c("C" = ncells/2), n_w, n_h)

#### run this in case you want to have a single cell sample
DLP.sample(sim, bbox1$lower_corner, bbox1$upper_corner, sample_prefix="DLP2_A")
DLP.sample(sim, bbox2$lower_corner, bbox2$upper_corner, sample_prefix="DLP2_B")
sim$get_samples_info()$tumour_cells %>% table()

### run this in case you want a bulk sample
sim$sample_cells("DLP2_A", bbox1$lower_corner, bbox1$upper_corner)
sim$sample_cells("DLP2_B", bbox2$lower_corner, bbox2$upper_corner)

plot_tissue(sim,color_map = c("A"="darkseagreen4","B"="lightpink2","C"="cadetblue3"))+
  geom_rect(xmin = bbox1$lower_corner[1], xmax = bbox1$upper_corner[1],
            ymin = bbox1$lower_corner[2], ymax = bbox1$upper_corner[2],
            fill = NA, color = "goldenrod2")+
  geom_rect(xmin = bbox2$lower_corner[1], xmax = bbox2$upper_corner[1],
            ymin = bbox2$lower_corner[2], ymax = bbox2$upper_corner[2],
            fill = NA, color = "forestgreen")

forest <- sim$get_sample_forest()
forest$save("/path/to/your/sample_forest.sff")
### run this to visualize forest
plot_forest(forest,highlight_sample = T,color_map = c("A"="darkseagreen4","B"="lightpink2","C"="cadetblue3"))+theme(legend.position = "none")

### set the working directory where you generated the reference
### if it is the first time generating the reference you do not need to set
### any wd but run directly the MutationEngine set up
### For doind so you need to be subscribed to COSMIC

setwd("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/scATAC_project/ProCESS-scATAC/scripts/simulate_fragments/peak_based/0_process_simulations/process_references_v1.3.5/")
m_engine <- MutationEngine(setup_code = "GRCh38",tumour_type = "COADREAD", context_sampling = 20,
                           COSMIC_account = list("email"="email address","password"="passwd"))


mu_SNV = 1e-8
mu_CNA = 1e-9
mu_INDELs = 1e-9

CNA_Clone2 = ProCESS::CNA(type = "D", "5",
                          from = 107707518, len = 2e7,allele = 0)

## Drivers for the tumors
m_engine$add_mutant(mutant_name = "A",
                    passenger_rates = c(SNV = mu_SNV, CNA = mu_CNA,indel=mu_INDELs),drivers = list(list("APC R1450*", allele = 1)))
m_engine$add_mutant(mutant_name = "B",passenger_rates = c(SNV = mu_SNV, CNA = 0,indel=mu_INDELs),drivers = list(CNA_Clone2))
m_engine$add_mutant(mutant_name = "C",passenger_rates = c(SNV = mu_SNV, CNA = 0,indel=mu_INDELs),drivers = list("KRAS G12D"))

m_engine$add_exposure(time = 0,coefficients = c(SBS1 = 0.15,SBS5 = 0.40,
                                                SBS18 = 0.15,SBS17b = 0.20,ID1 = 0.40,ID2 = 0.40,ID18=0.2,SBS88 = 0.10))
phylo_forest <- m_engine$place_mutations(forest, num_of_preneoplatic_SNVs=800, num_of_preneoplatic_indels=200)
phylo_forest$save("/path/to/your/phylo_forest.sff")



####### Where you can retrieve information about CNA?
### get allele specific CN segments only for sampled cells
phylo_forest$get_cell_allelic_fragmentation() %>% head()

##### this function gives you all the nodes in the forest (both sampled and internal)
phylo_forest$get_nodes() %>% head()

#### get CN history of an internal node
cell_id=6 ## cell id from the previous function
genome <- phylo_forest$get_node(cell_id)$get_genome()
head(genome$get_CNAs())


