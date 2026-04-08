rm(list=ls())
library(ProCESS)
library(dplyr)
library(ggplot2)
source("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-DLP/scripts/utils/DLP.R")
source("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-DLP/scripts/utils/utils_plot.R")
set.seed(12345)

setwd("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-DLP/dlp_simulation_4/")
sim <- TissueSimulation(name = "DLP4",save_snapshots = T)
sim$add_mutant(name = "A", growth_rates = 0.1, death_rates = 0)
sim$place_cell("A", 500, 500)
sim$run_up_to_size(species = 'A', 3000)

# sampling tissue
n_w <- n_h <- 20
ncells <- 50

# adding second mutant
sim$add_mutant(name = "B", growth_rates = 0.5, death_rates = 0.0)
#sim$update_rates(species = "A", rates = c(growth = 0.01, death = 0.1))
sim$mutate_progeny(sim$choose_cell_in("A"), "B")
sim$run_up_to_size(species = 'B', 500)

# adding third mutant
sim$add_mutant(name = "C", growth_rates = 0.8, death_rates = 0.0)
sim$mutate_progeny(sim$choose_cell_in("A"), "C")
sim$run_up_to_size(species = 'C', 1800)
# find a tissue rectangle containing 5 cells of type A and B at least
bbox1 <- sim$search_sample(c("A"=ncells), n_w, n_h)
bbox2 <- sim$search_sample(c("C" = ncells), n_w, n_h)
#'
DLP.sample(sim, bbox1$lower_corner, bbox1$upper_corner, sample_prefix="DLP4_A")
DLP.sample(sim, bbox2$lower_corner, bbox2$upper_corner, sample_prefix="DLP4_B")

sim$add_mutant(name = "D", growth_rates = 1.5, death_rates = 0.0)


sim$mutate_progeny(sim$choose_cell_in("B"), "D")
sim$update_rates(species = "A", rates = c(growth = 0.06, death = 0.01))

sim$run_up_to_size(species = 'D', 3000)
bbox3 <- sim$search_sample(c("D" = ncells), n_w, n_h)


DLP.sample(sim, bbox3$lower_corner, bbox3$upper_corner, sample_prefix="DLP4_C")
sim$get_samples_info()$tumour_cells %>% table()

forest <- sim$get_sample_forest()
# all_nodes <- forest$get_nodes() %>% 
#   filter(!is.na(sample))
# sampled_mutants <- all_nodes %>% 
#   separate(col = sample,into = c("p1","s1","c1"),sep = "_") %>% 
#   mutate(sample_name=paste(p1,s1,sep="_")) %>% 
#   group_by(sample_name,mutant) %>% 
#   summarise(prop=n())

forest$save("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-DLP/dlp_simulation_4/sample_forest_4.sff")
## 

color_palette = c("A"="darkseagreen4","B"="lightpink2","C"="cadetblue3")
plot_forest(forest,highlight_sample = T)+theme(legend.position = "none")
my_plot_forest(forest = forest,highlight_sample = T,color_map = color_palette,horizontal = F,color_sample = F)+
  theme(legend.position = "none")
## 
#
#
setwd("/orfeo/cephfs/scratch/cdslab/ggandolfi/process_on_the_fly_muts/reference")
m_engine <- MutationEngine(setup_code = "GRCh38",tumour_type = "COADREAD", context_sampling = 20)


mu_SNV = 1e-8
mu_CNA = 1e-12
mu_INDELs = 1e-9

CNA_Clone2 = ProCESS::CNA(type = "D", "5",
                          from = 107707518, len = 2e7,allele = 0)

CNA_Clone4 = ProCESS::CNA(type='A', chr='8', from=46000001, len=99138636)



## Drivers for the tumors
m_engine$add_mutant(mutant_name = "A",
                    passenger_rates = c(SNV = mu_SNV, CNA = mu_CNA,indel=mu_INDELs),drivers = list(list("APC R1450*", allele = 1)))
m_engine$add_mutant(mutant_name = "B",passenger_rates = c(SNV = mu_SNV, CNA = 0,indel=mu_INDELs),drivers = list(CNA_Clone2))
m_engine$add_mutant(mutant_name = "C",passenger_rates = c(SNV = mu_SNV, CNA = mu_CNA,indel=mu_INDELs),drivers = list("KRAS G12D"))
m_engine$add_mutant(mutant_name = "D",passenger_rates = c(SNV = mu_SNV, CNA = 0,indel=mu_INDELs),drivers = list(CNA_Clone4))

m_engine$add_exposure(time = 0,coefficients = c(SBS1 = 0.15,SBS5 = 0.40,
                                                SBS18 = 0.15,SBS17b = 0.20,ID1 = 0.40,ID2 = 0.40,ID18=0.2,SBS88 = 0.10))
phylo_forest <- m_engine$place_mutations(forest, num_of_preneoplatic_SNVs=800, num_of_preneoplatic_indels=200)
phylo_forest$save("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-DLP/dlp_simulation_4/phylo_forest_4.sff")

cov <- DLP.coverage(phylo_forest, 100, sample_prefix = "DLP4_")


no_error_seq <- ErrorlessIlluminaSequencer()

chromosomes <- phylo_forest$get_absolute_chromosome_positions()$chr

results <- parallel::mclapply(
  chromosomes,
  function(chr) {
    simulate_seq(
      phylo_forest,
      coverage = cov,
      write_SAM = TRUE,
      with_normal_sample = FALSE,chromosomes = c(chr),
      read_size = 150,
      sequencer = no_error_seq,
      insert_size_mean = 350,
      insert_size_stddev = 10,
      missed_SID_statistics=F, germline_statistics=F,
      wide_format=FALSE,
      output_dir = "/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-DLP/sequencing_dlp4",
      update_SAM = TRUE
    )
  },
  mc.cores = 4
)

saveRDS(object = results,file = "/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-DLP/dlp_simulation_4/seq_res_4.rds")

