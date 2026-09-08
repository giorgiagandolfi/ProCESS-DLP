rm(list = ls())
library(ProCESS)
library(dplyr)
set.seed(06117)


outdir <- "/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-DLP/SCOUT/SPN01/process/"
forest <- load_sample_forest(paste0(outdir,"sample_forest.sff"))


setwd("/orfeo/cephfs/scratch/cdslab/ggandolfi/process_on_the_fly_muts/reference")
m_engine <- MutationEngine(setup_code = "GRCh38",tumour_type = "COADREAD", context_sampling = 20)
                           # COSMIC_account = list("email"="giorgia.gandolfi@phd.units.it","password"="2*db!XQ4sgQ!dbg"))


mu_SNV = 1e-8
mu_CNA = 2e-10
mu_INDELs = 1e-9
##112707518-112846239 
CNA_Clone2 = ProCESS::CNA(type = "D", "5",
                         from = 107707518, len = 2e7,allele = 0)

## Drivers for the tumors
m_engine$add_mutant(mutant_name = "Clone 1",
                    passenger_rates = c(SNV = mu_SNV, CNA = 0,indel=mu_INDELs),drivers = list(list("APC R1450*", allele = 1)))
m_engine$add_mutant(mutant_name = "Clone 2",passenger_rates = c(SNV = mu_SNV, CNA = mu_CNA,indel=mu_INDELs),drivers = list(CNA_Clone2))
mu_SNV = 1e-8
mu_CNA = 1e-12
m_engine$add_mutant(mutant_name = "Clone 3",passenger_rates = c(SNV = mu_SNV, CNA = mu_CNA,indel=mu_INDELs),drivers = list("TP53 R175H"))

m_engine$add_mutant(mutant_name = "Clone 4",passenger_rates = c(SNV = mu_SNV, CNA = mu_CNA,indel=mu_INDELs),drivers = list(WGD,"PIK3CA R88Q"))


# Mutational signatures
m_engine$add_exposure(time = 0,coefficients = c(SBS1 = 0.15,SBS5 = 0.40,
                                                SBS18 = 0.15,SBS17b = 0.20,ID1 = 0.40,ID2 = 0.40,ID18=0.2,SBS88 = 0.10))
print("Mutation engine created")
phylo_forest <- m_engine$place_mutations(forest, num_of_preneoplatic_SNVs=800, num_of_preneoplatic_indels=200)
phylo_forest$save(paste0(outdir,"phylo_forest.sff"))

