validate_sample_prefix <- function(sample_prefix) {
  
  if (is.null(sample_prefix)) {
    return("")
  }
  
  if (!is.character(sample_prefix)) {
    stop("The parameter \"sample_prefix\" must be either NULL or a string.")
  }

  if (length(sample_prefix) > 0) {
    last_char <- substr(sample_prefix, nchar(sample_prefix),
                        nchar(sample_prefix))
    if (last_char != "_" && last_char != " ") {
      sample_prefix <- paste0(sample_prefix, "_")
    }
  }
  
  return(sample_prefix)
}


#' Simulate DLP+ sampling 
#'
#' @description
#' This function simulates DLP+ sampling by collecting a sample for
#' each of the cells in a tissue rectangle. ProCESS supports multiple
#' DLP+ samples, however, the overall number of collected cells must
#' be less than 2^16. 
#'
#' @param tissue_simulation The ProCESS tissue simulation from which the 
#'   DLP+ sample must be collected.
#' @param bottom_left The coordinate of the bottom left corner of the
#'   DLP+ sample.
#' @param top_right The coordinate of the top right corner of the
#'   DLP+ sample.
#' @param sample_prefix A prefix for the DLP+ sample (default: `NULL`).
#' @export
#'
#' @examples
#' set.seed(0)
#' 
#' sim <- TissueSimulation()
#' sim$add_mutant(name = "A", growth_rates = 0.1, death_rates = 0.0)
#' sim$place_cell("A", 500, 500)
#' sim$run_up_to_time(100)
#'
#' # sampling tissue
#' n_w <- n_h <- 10
#' ncells <- 5
#'
#' # adding second mutant
#' sim$add_mutant(name = "B", growth_rates = 0.3, death_rates = 0.0)
#' sim$mutate_progeny(sim$choose_cell_in("A"), "B")
#' sim$run_up_to_time(300)
#'
#' # find a tissue rectangle containing 5 cells of type A and B at least
#' bbox <- sim$search_sample(c("A" = ncells, "B" = ncells), n_w, n_h)
#'
#' # collect a DLP+ sample from the tissue simulation `sim`
#' DLP.sample(sim, bbox$lower_corner, bbox$upper_corner, sample_prefix="test")
#'
#' @seealso DLP.coverage()
DLP.sample <- function(tissue_simulation, bottom_left, top_right,
                       sample_prefix = NULL) {
  if (class(sim)[1] != "Rcpp_TissueSimulation") {
    stop(paste("The parameter \"tissue_simulation\" must be a",
               "ProCESS's TissueSimulation object."))
  }
  
  sample_prefix <- validate_sample_prefix(sample_prefix)
  
  if (length(bottom_left) != 2 || bottom_left[1] < 0 || bottom_left[2] < 0) {
    stop(paste("The parameter \"bottom_left\" must be a 2-dimensional vector",
               "containing the coordinates of the bottom left cell of",
               "the DLP sample."))
  }

  if (length(top_right) != 2 || top_right[1] < 0 || top_right[2] < 0) {
    stop(paste("The parameter \"top_right\" must be a 2-dimensional vector",
               "containing the coordinates of the top right cell of",
               "the DLP sample."))
  }
  
  num_of_samples <- nrow(tissue_simulation$get_samples_info())
  
  num_of_future_samples <- ((top_right[1]-bottom_left[1]+1)
                            * (top_right[2]-bottom_left[2]+1))
  
  if (num_of_samples+num_of_future_samples >= 2^(16)) {
    stop(paste("Too many cells in the sample. Currently, ProCESS supports",
               "up to 2^16 cells in DLP samples."))
  }
  
  for (y in seq(bottom_left[2], top_right[2])) {
    for (x in seq(bottom_left[1], top_right[1])) {
      cell_pos <- c(x, y)
      cells <- tissue_simulation$get_cells(c(x, y), c(x,y))
      
      if (nrow(cells)>0) { # if the cell is a tumour cell
        cell_sample_name <- paste0(sample_prefix, x, "-", y)
        tissue_simulation$sample_cells(cell_sample_name, cell_pos, cell_pos)
      }
    }
  }
}

#' Compute DLP+ adjusted coverage 
#'
#' @description
#' This function computes the coverage to use in `simulate_seq()` to
#' properly simulate a DLP+ sequencing in ProCESS.DLP+ sequencing
#' is simulated in ProCESS by associating each DPL+ sample cell to a
#' sample. However, the coverage parameter in `simulate_seq()` refers to
#' the coverage per sample, i.e., the coverage per cell. Thus, we need
#' to compute the coverage of the ProCESS samples before sequencing a
#' full DPL+ sample. Luckly, the coverage is uniform over cells in DLP+.
#' Hence, we can compute the cell coverage as the DPL+ sample coverage
#' divided by the number of cells in the DLP+ sample.
#'
#' @param phylo_forest The ProCESS phylogenetic forest from which the 
#'   DLP+ sample must be collected.
#' @param overall_coverage The aimed coverage for the entire DPL+ sample.
#' @param sample_prefix The prefix of a DLP+ sample (default: `NULL`).
#' @export
#'
#' @examples
#' set.seed(0)
#'
#' sim <- TissueSimulation()
#' sim$add_mutant(name = "A", growth_rates = 0.1, death_rates = 0.0)
#' sim$place_cell("A", 500, 500)
#' sim$run_up_to_time(100)
#'
#' # sampling tissue
#' n_w <- n_h <- 10
#' ncells <- 5
#'
#' # adding second mutant
#' sim$add_mutant(name = "B", growth_rates = 0.3, death_rates = 0.0)
#' sim$mutate_progeny(sim$choose_cell_in("A"), "B")
#' sim$run_up_to_time(300)
#'
#' # find a tissue rectangle containing 5 cells of type A and B at least
#' bbox <- sim$search_sample(c("A" = ncells, "B" = ncells), n_w, n_h)
#'
#' # collect a DLP+ sample from the tissue simulation `sim`
#' DLP.sample(sim, bbox$lower_corner, bbox$upper_corner, sample_prefix="test")
#'
#' forest <- sim$get_sample_forest()
#'
#' rm(sim)
#'
#' # placing mutations
#' m_engine <- MutationEngine(setup_code = "demo")
#'
#' m_engine$add_mutant(mutant_name="A", passenger_rates=c(SNV=5e-8),
#'                     drivers = list(SNV("22", 46510210, "C", "A", allele = 1),
#'                                    "DGCR8 P26L"))
#' m_engine$add_mutant(mutant_name="B", passenger_rates=c(SNV=5e-9),
#'                     drivers = list(list("DGCR8 A18V", allele = 1)))
#' m_engine$add_exposure(c(SBS1 = 0.2, SBS5 = 0.8))
#'
#' phylo_forest <- m_engine$place_mutations(forest, 10, 10)
#'
#' rm(forest)
#'
#' # compute adjusted coverage for simulating DLP+ sequencing
#' cov <- DLP.coverage(phylo_forest, 10, sample_prefix = "test")
#' cov
#'
#' # simulating DLP+ sequencing
#' seq_results <- simulate_seq(phylo_forest, coverage = cov, write_SAM = F,
#'                             with_normal_sample = FALSE)
#'
#' head(seq_results$mutations)
#' @seealso DLP.sample()
DLP.coverage <- function(phylo_forest, overall_coverage,
                         sample_prefix = NULL) {
  
  if (class(phylo_forest)[1] != "Rcpp_PhylogeneticForest") {
    stop(paste("The parameter \"phylo_forest\" must be a",
               "ProCESS's PhylogeneticForest object."))
  }
  
  if (class(overall_coverage) != "numeric") {
    stop("The parameter \"overall_coverage\" must be a numeric value.")
  }

  sample_prefix <- validate_sample_prefix(sample_prefix)

  sample_info <- phylo_forest$get_samples_info()
  
  num_of_DLP_samples <- sum(grepl(paste0("^", sample_prefix),
                                  sample_info$name))
  rm(sample_info)
  
  return(overall_coverage/num_of_DLP_samples)
}