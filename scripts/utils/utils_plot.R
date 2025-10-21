my_plot_forest = function (forest, highlight_sample = NULL, color_map = NULL, horizontal = FALSE,color_sample=NULL) {
  stopifnot(inherits(forest, "Rcpp_SampleForest"))
  nodes <- forest$get_nodes()
  if (nrow(nodes) == 0) {
    warning("The forest does not contain any node")
    return(ggplot2::ggplot())
  } else {
    forest_data <- forest$get_nodes()
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
                     shape = ifelse(is.na(.data$sample), "N/A", .data$sample),
                     size  = .data$sample)
      ) +
      ggplot2::scale_shape_manual(values = c(0:nsamples + 1)) +
      ggplot2::scale_color_manual(values = color_map) +
      ggplot2::theme_minimal() +
      ggplot2::theme(legend.position = "bottom") +
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
    p <- p +
      ggraph::geom_node_point(
        data = dplyr::filter(layout, sample != "N/A"),
        inherit.aes = FALSE,
        ggplot2::aes(x = x, y = y, fill = sample),  # <-- map sample to fill
        pch = 21,          # circle with fill
        stroke = 0,      # border thickness
        color = "black",   # border color stays black
        size = 3           # fixed circle size
      ) +
      ggplot2::scale_fill_manual(
        values = setNames(color_sample,
                          forest$get_samples_info() %>% dplyr::pull(.data$name))
      )
    
    p = p +
      labs(color = "Clone", fill = "Sample")
    
    if (horizontal) {
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
          legend.position = "bottom"
        )
    }
    
    return(p)
  }
}