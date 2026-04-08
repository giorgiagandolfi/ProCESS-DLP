library(tidyr)
my_plot_forest = function (forest, highlight_sample = NULL, color_map = NULL, horizontal = FALSE,color_sample=NULL,color_real_sample=NULL) {
  stopifnot(inherits(forest, "Rcpp_SampleForest"))
  nodes <- forest$get_nodes()
  if (nrow(nodes) == 0) {
    warning("The forest does not contain any node")
    return(ggplot2::ggplot())
  } else {
    forest_data <- forest$get_nodes()
    sampled_cells <- forest$get_nodes() %>%
      filter(!is.na(sample)) %>% 
      separate(col = sample,into = c("p1","s1","c1"),sep = "_",remove = F) %>%
      mutate(real_sample=paste(p1,s1,sep="_"))
    # sampled_cells$color <- color_real_sample[sampled_cells$real_sample]
    forest_data[nrow(forest_data) + 1, ] <- c(NA, NA, NA, NA, NA, 0)
    forest_data <- forest_data %>%
      dplyr::as_tibble() %>%
      dplyr::rename(from = .data$ancestor, to = .data$cell_id) %>%
      dplyr::select(.data$from, .data$to, .data$mutant, .data$epistate, .data$sample, .data$birth_time) %>%
      dplyr::mutate(from = ifelse(is.na(.data$from), "WT", .data$from),
                    to   = ifelse(is.na(.data$to), "WT", .data$to),
                    species = paste0(.data$mutant, .data$epistate),
                    sample  = ifelse(is.na(.data$sample), "N/A", .data$sample),
                    highlight = FALSE)
    
    if (!is.null(highlight_sample)) {
      highlight <- ProCESS:::paths_to_sample(forest_data, highlight_sample)
      forest_data$highlight <- forest_data$to %in% highlight$to
    }
    
    edges <- forest_data %>% dplyr::select("from", "to", "highlight")
    graph <- tidygraph::as_tbl_graph(edges, directed = TRUE)
    graph <- graph %>%
      tidygraph::activate("nodes") %>%
      dplyr::left_join(forest_data %>%
                         dplyr::rename(name = .data$to) %>%
                         dplyr::mutate(name = as.character(.data$name)),
                       by = "name")
    
    layout <- ggraph::create_layout(graph, layout = "tree", root = "WT")
    max_Y <- max(layout$birth_time, na.rm = TRUE)
    
    layout$y <- layout$birth_time
    
    if (is.null(color_map)) {
      color_map <- ProCESS:::get_species_colors(forest$get_species_info())
    }
    
    nsamples <- forest$get_samples_info() %>% nrow()
    labels_every <- max_Y / 10
    point_size <- c(0.5, rep(1, nsamples))
    names(point_size) <- c("N/A", forest$get_samples_info() %>% dplyr::pull(.data$name))
    
    group_name <- ProCESS:::get_group_cell_name(forest)
    
    # --- Base tree edges
    graph_plot <- ggraph::ggraph(layout, "tree") +
      ggraph::geom_edge_link(edge_width = 0.1,
                             ggplot2::aes(edge_color = ifelse(highlight, "indianred3", "black")))
    
    # --- Base nodes: colored by species
    p <- graph_plot +
      ggraph::geom_node_point(
        ggplot2::aes(color = .data$species,
                     # shape = ifelse(is.na(.data$sample), "N/A", .data$sample),
                     size  = .data$sample)
      ) +
      # ggplot2::scale_color_manual(values = c(0:nsamples + 1)) +
      ggplot2::scale_color_manual(values = color_map) +
      # ggplot2::scale_color_manual(values = color_real_sample) +
      ggplot2::theme_minimal() +
      ggplot2::theme(legend.position = "none") +
      ggplot2::labs(color = group_name, shape = "Sample", x = NULL, y = "Time") +
      ggplot2::guides(size = "none",
                      shape = ggplot2::guide_legend("Sample"),
                      color = ggplot2::guide_legend(group_name)) +
      ggplot2::scale_size_manual(values = point_size) +
      ggplot2::scale_y_continuous(labels = seq(0, max_Y, labels_every) %>% round,
                                  breaks = seq(0, max_Y, labels_every) %>% round) +
      ggplot2::theme(axis.line.x = ggplot2::element_blank(),
                     axis.text.x = ggplot2::element_blank(),
                     axis.ticks.x = ggplot2::element_blank())
    
    # --- NEW LAYER: outline circles for Sample != "N/A"
    # p <- p +
    #   ggraph::geom_node_point(
    #     data = dplyr::filter(layout, sample != "N/A"),
    #     inherit.aes = FALSE,
    #     ggplot2::aes(x = x, y = y, fill = sample),  # <-- map sample to fill
    #     pch = 21,          # circle with fill
    #     stroke = 0,      # border thickness
    #     color = "black",   # border color stays black
    #     size = 2          # fixed circle size
    #   ) #+
      # ggplot2::scale_fill_manual(
      #   values = setNames(sampled_cells$color,
      #                     sampled_cells$sample)
      #) 
    
    
    p = p +
      ggplot2::labs(color = "Clone", fill = "Sample")
    
    if (horizontal==F) {
      p = p + scale_y_reverse()
    } else {
      p = p +
        coord_flip() +
        theme_minimal() +
        theme(
          axis.ticks.y = element_blank(),
          axis.text.y  = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank(),
          legend.position = "none"
        )
    }
    
    return(p)
  }
}

get_clone_map <- function(sample_forest){
  n_clones <- nrow(sample_forest$get_species_info())
  clones <- sample_forest$get_species_info() %>% 
    pull(mutant)
  clone_colors <- RColorBrewer::brewer.pal(n_clones, "Dark2")
  names(clone_colors) <- clones
  return(clone_colors)
}

plot_DLP_state <- function(sample_forest){
  all_nodes <- sample_forest$get_nodes() %>%
    filter(!is.na(sample))
  sampled_mutants <- all_nodes %>%
    separate(col = sample,into = c("p1","s1","c1"),sep = "_") %>%
    mutate(sample_name=paste(p1,s1,sep="_")) %>%
    group_by(sample_name,mutant) %>%
    summarise(prop=n())
  
  p <- ggplot(sampled_mutants, aes(x = "", y = prop, fill = mutant)) +
    geom_bar(stat = "identity", width = 1) +
    scale_fill_manual(values=get_clone_map(sample_forest = sample_forest))+
    coord_polar("y") +
    facet_wrap(~ sample_name,scales = "free",nrow = 1) +
    theme_void() +
    labs(title = "Proportion of Mutants per Sample") +
    theme(
      plot.title = element_text(hjust = 0.5),
      legend.title = element_blank()
    )
  return(p)
}


my_annotate_forest <- function (tree_plot, forest, samples = TRUE, MRCAs = TRUE, exposures = FALSE, 
                                facet_signatures = TRUE, drivers = TRUE, add_driver_label = TRUE) 
{
  samples_info <- forest$get_samples_info()
  if (samples) {
    max_Y <- max(tree_plot$data$y, na.rm = TRUE)
    tree_plot <- tree_plot + ggplot2::geom_hline(yintercept = max_Y - 
                                                   samples_info$time, color = "indianred3", linetype = "dashed", 
                                                 linewidth = 0.3)
  }
  if (MRCAs) {
    sample_names <- samples_info %>% dplyr::pull(.data$name)
    MRCAs_cells <- lapply(sample_names, function(s) {
      forest$get_coalescent_cells(forest$get_nodes() %>% 
                                    dplyr::filter(sample %in% s) %>% dplyr::pull(.data$cell_id)) %>% 
        dplyr::mutate(sample = s)
    }) %>% Reduce(f = dplyr::bind_rows) %>% dplyr::group_by(.data$cell_id) %>% 
      dplyr::mutate(cell_id = paste(.data$cell_id)) %>% 
      dplyr::summarise(label = paste0("    ", .data$sample, 
                                      collapse = "\n"))
    layout <- tree_plot$data %>% dplyr::select(.data$x, .data$y, 
                                               .data$name) %>% dplyr::mutate(cell_id = paste(.data$name)) %>% 
      dplyr::filter(.data$name %in% MRCAs_cells$cell_id) %>% 
      dplyr::left_join(MRCAs_cells, by = "cell_id")
    tree_plot <- tree_plot + ggplot2::geom_point(data = layout, 
                                                 ggplot2::aes(x = .data$x, y = .data$y), color = "purple3", 
                                                 size = 3, pch = 21) + ggplot2::geom_text(data = layout, 
                                                                                          ggplot2::aes(x = .data$x, y = .data$y, label = .data$label), 
                                                                                          color = "purple3", size = 3, hjust = 0, vjust = 1)
  }
  if (exposures) {
    if (inherits(forest, "Rcpp_PhylogeneticForest")) {
      max_Y <- max(tree_plot$data$y, na.rm = TRUE)
      exposures <- forest$get_exposures()
      exposure_colors <- get_colors_for(exposures %>% dplyr::pull(signature) %>% 
                                          unique)
      times <- exposures$time %>% unique() %>% sort()
      exposures <- exposures %>% dplyr::rowwise() %>% dplyr::mutate(t_end = dplyr::case_when(time == 
                                                                                               max(times) ~ Inf, .default = min(times[times >= 
                                                                                                                                        time]))) %>% dplyr::mutate(signature = factor(signature, 
                                                                                                                                                                                      levels = exposures %>% dplyr::arrange(time) %>% 
                                                                                                                                                                                        dplyr::pull(signature) %>% unique()))
      breaks <- sort(unique(exposures$exposure))
      tree_plot <- tree_plot + ggplot2::geom_rect(data = exposures, 
                                                  ggplot2::aes(xmin = -Inf, xmax = Inf, ymin = ifelse(is.infinite(t_end), 
                                                                                                      0, max_Y - t_end), ymax = max_Y - time, fill = signature, 
                                                               alpha = exposure)) + ggplot2::scale_fill_manual(values = exposure_colors) + 
        ggplot2::scale_alpha_continuous(range = c(0.25, 
                                                  0.75), breaks = breaks) + ggplot2::guides(fill = ggplot2::guide_legend(title = "Signature"), 
                                                                                            alpha = ggplot2::guide_legend(title = "Exposure"))
      if (facet_signatures) {
        tree_plot <- tree_plot + ggplot2::facet_wrap(~signature)
      }
      layers_new <- list(tree_plot$layers[[length(tree_plot$layers)]])
      layers_new <- c(layers_new, tree_plot$layers[1:(length(tree_plot$layers) - 
                                                        1)])
      tree_plot$layers <- layers_new
    }
  }
  if (drivers) {
    if (inherits(forest, "Rcpp_PhylogeneticForest")) {
      drivers_mutations = drivers_CNAs = data.frame()
      try(expr = {
        drivers_mutations <- forest$get_sampled_cell_mutations() %>% 
          dplyr::filter(class == "driver") %>% dplyr::mutate(driver_id = paste0(chr, 
                                                                                ":", from, ":", ref, ">", alt), driver_type = type) %>% 
          dplyr::select(cell_id, driver_id, driver_type)
      })
      try(expr = {
        drivers_CNAs <- forest$get_sampled_cell_CNAs() %>% 
          dplyr::filter(class == "driver") %>% dplyr::mutate(driver_id = paste0(chr, 
                                                                                ":", begin, "-", end, ":", allele), driver_type = "CNA") %>% 
          dplyr::select(cell_id, driver_id, driver_type)
      })
      drivers <- dplyr::bind_rows(drivers_mutations, drivers_CNAs)
      drivers_start_nodes <- lapply(unique(drivers$driver_id), 
                                    function(d) {
                                      nodes_with_driver <- drivers %>% dplyr::filter(driver_id == 
                                                                                       d) %>% dplyr::pull(cell_id)
                                      d_type <- drivers %>% dplyr::filter(driver_id == 
                                                                            d) %>% dplyr::pull(driver_type) %>% unique()
                                      forest$get_coalescent_cells(nodes_with_driver) %>% 
                                        dplyr::mutate(driver_id = d, driver_type = d_type)
                                    }) %>% dplyr::bind_rows() %>% dplyr::mutate(cell_id = as.character(cell_id)) %>% 
        dplyr::group_by(cell_id) %>% dplyr::summarise(driver_id = paste0(driver_id, 
                                                                         collapse = "\n"))
      layout <- tree_plot$data %>% dplyr::select(x, y, 
                                                 mutant,
                                                 name) %>% dplyr::mutate(cell_id = paste(name), 
                                                                         has_driver = TRUE) %>% dplyr::filter(name %in% 
                                                                                                                drivers_start_nodes$cell_id) %>% dplyr::left_join(drivers_start_nodes, 
                                                                                                                                                                  by = "cell_id")
      layout_sample <- tree_plot$data %>% dplyr::select(x, y, 
                                                 mutant,
                                                 name)
      tree_plot <- tree_plot + ggplot2::geom_point(data = layout, 
                                                   ggplot2::aes(x = .data$x, y = .data$y,fill=.data$mutant),
                                                   color = "#FF000000", size = 2, pch = 21)+
        scale_fill_manual(values = color_map)
      if (add_driver_label) {
        nudge_x <- (max(tree_plot$data$x) - min(tree_plot$data$x)) * 
          0.15
        tree_plot <- tree_plot + ggrepel::geom_label_repel(data = layout,
                                                           ggplot2::aes(x = .data$x, y = .data$y, label = .data$driver_id,fill=.data$mutant), 
                                                           size = 2.5, hjust = 0, nudge_x = -nudge_x, 
                                                           direction = "x")+
          scale_fill_manual(values = color_map)
      }
    }
  }
  if (passenger){
    
  }
  tree_plot
}


# passenger_mutations = passenger_CNAs = data.frame()
# 
# try(expr = {
#   passenger_CNAs <- forest$get_sampled_cell_CNAs() %>% 
#     dplyr::filter(class == "passenger") %>% dplyr::mutate(passenger_id = paste0(chr, 
#                                                                           ":", begin, "-", end, ":", allele), passenger_type = "CNA") %>% 
#     dplyr::select(cell_id, passenger_id, passenger_type)
# })
# # drivers <- dplyr::bind_rows(drivers_mutations, drivers_CNAs)
# passenger_start_nodes <- lapply(unique(passenger_CNAs$passenger_id), 
#                               function(d) {
#                                 nodes_with_driver <- passenger_CNAs %>% dplyr::filter(passenger_id == 
#                                                                                  d) %>% dplyr::pull(cell_id)
#                                 d_type <- passenger_CNAs %>% dplyr::filter(passenger_id == 
#                                                                       d) %>% dplyr::pull(passenger_type) %>% unique()
#                                 forest$get_coalescent_cells(nodes_with_driver) %>% 
#                                   dplyr::mutate(passenger_id = d, passenger_type = d_type)
#                               }) %>% dplyr::bind_rows() %>% dplyr::mutate(cell_id = as.character(cell_id)) %>% 
#   dplyr::group_by(cell_id) %>% dplyr::summarise(passenger_id = paste0(passenger_id, 
#                                                                    collapse = "\n"))
# layout <- tree_plot$data %>% dplyr::select(x, y, 
#                                            mutant,
#                                            name) %>% dplyr::mutate(cell_id = paste(name), 
#                                                                    has_passenger = TRUE) %>% dplyr::filter(name %in% 
#                                                                                                           passenger_start_nodes$cell_id) %>% dplyr::left_join(passenger_start_nodes, 
#                                                                                                                                                             by = "cell_id")
# layout_sample <- tree_plot$data %>% dplyr::select(x, y, 
#                                                   mutant,
#                                                   name)
# tree_plot <- tree_plot + ggplot2::geom_point(data = layout, 
#                                              ggplot2::aes(x = .data$x, y = .data$y,fill=.data$mutant),
#                                              color = "#FF000000", size = 2, pch = 21)+
#   scale_fill_manual(values = color_map)
# nudge_x <- (max(tree_plot$data$x) - min(tree_plot$data$x)) * (-10)
#   #0.15
# 
#   tree_plot <- tree_plot + ggrepel::geom_label_repel(data = layout,
#                                                    ggplot2::aes(x = .data$x, y = .data$y, label = .data$passenger_id,color=.data$mutant), 
#                                                    size = 2.5, hjust = 0, nudge_x = -nudge_x, 
#                                                    fill="transparent",
#                                                    direction = "x")+
#   scale_color_manual(values = color_map)

