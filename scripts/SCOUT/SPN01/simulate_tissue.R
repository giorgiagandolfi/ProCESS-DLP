rm(list = ls())
library(ProCESS)
library(dplyr)
library(ggplot2)
library(dplyr)
library(tidyr)
source("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-DLP/scripts/utils/DLP.R")
source("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-DLP/scripts/utils/utils_plot.R")


dir <- getwd()
outdir <- "/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-DLP/SCOUT/SPN01/process/"
setwd(outdir)
set.seed(06117)
sim <- TissueSimulation(name = 'SPN01', seed = 12345, save_snapshot=T)
sim$history_delta <- 1
sim$death_activation_level <- 50



# Clone 1
sim$add_mutant(name = "Clone 1", growth_rates = 0.08, death_rates = 0.01)
sim$place_cell("Clone 1", 500, 500)
sim$run_up_to_size("Clone 1", 500)

# Clone 2
sim$add_mutant(name = "Clone 2", growth_rates = 0.3, death_rates = 0.01)
sim$update_rates(species = "Clone 1", rates = c(growth = 0.06, death = 0.01))
sim$mutate_progeny(sim$choose_cell_in("Clone 1"), "Clone 2")
sim$run_up_to_size(species = 'Clone 2', 700)
sim$update_rates(species = "Clone 1", rates = c(growth = 0.02, death = 0.03))
sim$update_rates(species = "Clone 2", rates = c(growth = 0.4, death = 0.005))
sim$run_up_to_size(species = 'Clone 2', 1000)


# Grow 2
sim$update_rates(species = "Clone 1", rates = c(growth = 0.001, death = 0.3))
sim$run_up_to_size(species = 'Clone 2', 2000)


# Clone 3
sim$add_mutant(name = "Clone 3", growth_rates = 1, death_rates = 0.01)
sim$update_rates(species = "Clone 2", rates = c(growth = 0.1, death = 0.01))
sim$mutate_progeny(sim$choose_cell_in("Clone 2"), "Clone 3")
sim$run_up_to_size("Clone 3", 4000)
sim$update_rates(species = "Clone 2", rates = c(growth = 0.05, death = 0.5))
sim$run_up_to_size("Clone 3", 6000)
sim$update_rates(species = "Clone 2", rates = c(growth = 0.01, death = 0.8))
sim$run_up_to_size("Clone 3", 8000) ## orignal 15000

sim$update_rates(species = "Clone 2", rates = c(growth = 0.01, death = 1))
sim$run_up_to_size("Clone 3", 10000)

# Clone 4
sim$add_mutant(name = "Clone 4", growth_rates = 2.5 , death_rates = 0.01)
sim$update_rates(species = "Clone 2", rates = c(growth = 0.005, death = 1))
sim$update_rates(species = "Clone 3", rates = c(growth = 0.5, death = 0.02))
sim$mutate_progeny(sim$choose_cell_in("Clone 3"), "Clone 4")
sim$run_up_to_size("Clone 4", 4000) ## original 6000
sim$update_rates(species = "Clone 3", rates = c(growth = 0.02, death = 0.02))
sim$run_up_to_size("Clone 4", 6000) ##60000 ## orignal 8000
#### try this
sim$update_rates(species = "Clone 3", rates = c(growth = 0.01, death = 0.01))
sim$run_up_to_size("Clone 4", 8000)
sim$run_up_to_size("Clone 4", 15000)

print("End simulation")
## Sample 3
bboxC_lower_corner <- c(300, 400)
bboxC_upper_corner <- c(340, 440)
current <- plot_tissue(sim)
current +
  geom_rect(xmin = bboxC_lower_corner[1], xmax = bboxC_upper_corner[1],
            ymin = bboxC_lower_corner[2], ymax = bboxC_upper_corner[2],
            fill = NA, color = "black")


DLP.sample(sim, bboxC_lower_corner, bboxC_upper_corner, sample_prefix="SPN01_1.3")

sim$get_samples_info()$tumour_cells %>% table()

## Sample 1
bboxA <- sim$search_sample(c("Clone 4" = 200, 'Clone 3' = 200), 40,40)
DLP.sample(sim, bboxA$lower_corner, bboxA$upper_corner, sample_prefix="SPN01_1.1")


bboxB_lower_corner <- c(420,500)
bboxB_upper_corner <- c(460,540)

DLP.sample(sim, bboxB_lower_corner, bboxB_upper_corner, sample_prefix="SPN01_1.2")



# # Forest
forest <- sim$get_sample_forest()
#pie_chart <- plot_DLP_state(sample_forest = forest)
forest$save(paste0(outdir,"sample_forest.sff"))
