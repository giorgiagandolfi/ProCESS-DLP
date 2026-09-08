library(dplyr)
library(tidyverse)
library(GenomicRanges)
library(circlize)
library(ComplexHeatmap)
library(ProCESS)

source("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-examples/report/plotting/utils.R")

plot_ProCESS_heatmap<-function(sample_forest,phylo_forest,output_file){
  samples <- phylo_forest$get_samples_info()$name
  bulk_cell_cna <- phylo_forest$get_cell_allelic_fragmentation() 
  sampled_cell_info <- sample_forest$get_nodes() %>% 
    dplyr::filter(!is.na(sample))
  
  bulk_cell_cna_annot <- left_join(bulk_cell_cna,sampled_cell_info,by="cell_id")
  
  x <- bulk_cell_cna_annot %>% 
    dplyr::mutate(chr=paste0("chr",chr)) %>% 
    dplyr::rename("from"="begin") %>% 
    dplyr::rename("to"="end") %>% 
    dplyr::rename(Major=major) %>% 
    dplyr::rename(sample_id=sample) %>% 
    dplyr::mutate(segment_id=paste(chr,from,to,sep=":"))
  
 
  
  out = lapply(x$chr %>% unique() , function(chr) {
    
    cli::cli_alert_info("Iterating on {.field {chr}}")
    cli::cli_h2("Iterating on {.field {chr}}")
    
    old_segments = sapply(x$sample_id %>% unique(), function(s) {
      x %>%
        dplyr::filter(sample_id == s) %>%
        dplyr::filter(chr == !!chr) %>%
        dplyr::pull(segment_id) %>%
        unique() %>%
        length()
    })
    
    cli::cli_alert_info("Number of original segments in individual {.cls CNAqc} objects in {.field {chr}}:")
    cli::cli_ul(paste(names(old_segments), old_segments, sep = " = "))
    cat("\n")
    
    # Chromosome-specific new breakpoints
    new_breakpoints = c(
      x %>%
        dplyr::filter(chr == !!chr) %>%
        dplyr::pull(from),
      x %>%
        dplyr::filter(chr == !!chr) %>%
        dplyr::pull(to)) %>% 
      unique() %>%
      sort()
    
    cli::cli_alert_info("Found {.field {length(new_breakpoints)}} breakpoints")
    #cat("\n")
    
    #  Separate new breakpoints into segment from and to values
    new_from = new_breakpoints[ !new_breakpoints == dplyr::last(new_breakpoints)] # last element will not be included in the from column 
    new_to = new_breakpoints[!new_breakpoints == dplyr::first(new_breakpoints)] # first element will not be included in the to column 
    
    # iterate on the new breakpoints to subset the cna piled up 
    lapply(new_breakpoints[-1] %>% seq_along(), function(i) {
      
      # iterate for each sample 
      lapply(x$sample_id %>% unique(), function(s) {
        
        # define which columns must be kept in the new table
        not_wanted = c("from", "to", "length", "size", "segment_id", "chr", "sample_id", "n")
        wanted = setdiff(colnames(x), not_wanted)
        
        # get the information for the sample in the new segment 
        tmp = get_segment_info(x, 
                               chr = chr, 
                               sample = s, 
                               new_from = new_from[i], 
                               new_to = new_to[i], 
                               keep_columns = wanted)
        
        # do some checking on the result
        if (nrow(tmp) == 0) { # there is no information on copy number on the new segment: insert NA as value of all the columns, except from, to, segment_id and sample_id
          
          tmp_v2 = rep(NA, ncol(tmp))  
          names(tmp_v2) = colnames(tmp)
          
          tmp = tmp_v2 %>% 
            tibble::as_tibble_row()
        }
        
        # create a tibble with the information on the new breakpoints and include the previously retrieved information 
        tidyr::tibble(chr = chr, 
                      from = new_from[i], 
                      to = new_to[i], 
                      sample_id = s) %>% 
          dplyr::bind_cols(tmp) %>% 
          dplyr::mutate(segment_id = paste(chr, from, to, sep = ":")) 
        
      }) %>% do.call(bind_rows, .)
    }) %>% do.call(bind_rows, .)
  })  %>% do.call(bind_rows, .) # create a unique big tibble with the new segmentation

  
  wide_df <- out %>%
    dplyr::mutate(total_cn=Major+minor) %>% 
    na.omit() %>%
    dplyr::select(chr, from, to, cell_id,total_cn) %>% # Select only relevant columns
    pivot_wider(names_from = cell_id, values_from = total_cn,values_fill = list(total_cn = 0),values_fn = list(total_cn = sum)) %>% as.data.frame() # Reshape to wide format
  
  wide_df_gr1 <- GRanges(seqnames = wide_df[, 1], ranges = IRanges(wide_df[, 2], wide_df[, 3]))
  
  
  
  chr_df = read.chromInfo()$df
  chr_df = chr_df[chr_df$chr %in% paste0("chr", 1:22), ]
 
  chr_gr = GRanges(seqnames = chr_df[, 1], ranges = IRanges(chr_df[, 2] + 1, chr_df[, 3]))
  
  library(EnrichedHeatmap)
  chr_window = makeWindows(chr_gr, w = 1e6)
  
  chr = as.vector(seqnames(chr_window))
  chr_level = paste0("chr", 1:22)
  
  chr = factor(chr, levels = chr_level)
  
  
  num_mat_mine = average_in_window(chr_window, wide_df_gr1, wide_df[, -(1:3)])
  num_mat_mine[is.na(num_mat_mine)] <- 0
  
  subgroup = out %>% select(cell_id,mutant,sample_id) %>% na.omit() %>% unique()
  
  
  subgroup <- subgroup %>%
    group_by(sample_id) %>%
    mutate(
      clonal_status = ifelse(n_distinct(mutant) > 1, "polyclonal", "monoclonal")
    ) %>%
    ungroup() %>% 
    separate(col = sample_id,into = c("p1","p2","p3"),sep = "_") %>% 
    mutate(sample_name=paste(p1,p2,sep="_"))
  ######## for karotype
  wide_df_kar <- out %>%
    # mutate(karyotype=paste(Major,minor,sep=":")) %>% 
    mutate(CN=Major+minor) %>% 
    na.omit() %>% 
    mutate(cell_id=as.character(cell_id)) %>% 
    select(chr, from, to, cell_id, CN) %>% # Select only relevant columns
    pivot_wider(names_from = cell_id, values_from = CN) %>% as.data.frame() # Reshape to wide format
  char_mat = NULL
  n_cells <- ncol(wide_df_kar)-3
  
  for(i in 1:n_cells) {
    reg = wide_df_kar[,1:3]
    j=i+3
    reg = cbind(reg,wide_df_kar[,j])
    # bed = bed[sample(nrow(bed), 20), , drop = FALSE]
    gr_cnv = GRanges(seqnames = reg[, 1], ranges = IRanges(reg[, 2], reg[, 3]))
    
    char_mat = cbind(char_mat, average_in_window(chr_window, gr_cnv, reg[, 4]))
    print(paste0(i,"/",n_cells))
  }
  
  
  
  # col_clone = get_clone_map(sample_forest = sample_forest)
  col_clone =c("A"="darkseagreen4","B"="lightpink2","C"="firebrick4")
  col_sample = c("DPL7_A"="goldenrod2")
  n_samples <- length(unique(subgroup$sample_name))
  #col_sample <- RColorBrewer::brewer.pal(n_samples, "Set2")
  #names(col_sample) = unique(subgroup$sample_id)
  #col_clonality = c("monoclonal" = "#ffff99","polyclonal"="#b15928")
  
  row_ha = rowAnnotation(clone = subgroup$mutant,
                         # sample = subgroup$sample_name,
                         # clonality=subgroup$clonal_status,
                         col=list(clone=col_clone))
  
  
  
  colors <- rep(c("gray", "black"), length.out = length(levels(chr)))
  
  # Create a named vector to assign colors to each chromosome level
  color_map <- setNames(colors, levels(chr))
  
  # Apply the color map to the factor to get the color for each chromosome in `chr`
  chr_colors <- color_map[chr]
  column_ha = HeatmapAnnotation(chromosome = chr, col=list(chromosome=chr_colors),
                                # label = anno_mark(at = at, labels = labels),
                                show_legend = F)
  my_palette <- c(
    "0"="#3590d2",
    "1"="#7fbbe3",
    "2" = "#d4cfcd",
    "3" = "#fdc38dff",
    "4" = "#fca16cff",
    "5" = "#f67b51ff",
    "6" = "#e7533aff",
    "7" = "#cf2518ff",
    "8" = "#ad0000ff",
    "9" = "#c4198aff"
  )
  
  #### COLOR FOR KARYOTYPE
  # col = c("1:1" = "#228B22CC", "1:0" = "steelblue",
  #         "0:0"="darkblue","2:0"="turquoise4",
  #         "2:1"="#FFA500CC","2:2"="firebrick3","Other"="grey")
  # karyotypes <- x %>% mutate(karyotype=paste(Major,minor,sep=":"))%>%
  #   pull(karyotype) %>% unique()
  # col <- get_karyotypes_colors(karyotypes = karyotypes)
  # names(col)%in%karyotypes

  ## eventually we can annotated drivers on the heatmap
  ## if the getter from the phylo_forest is present
  # drivers <- seq$mutations %>% filter(classes=="driver") %>% 
  #   rename(chr_pos="from") %>% 
  #   mutate(chr=paste0("chr",chr))
  # 
  # drivers_table <- read.table("/Users/gandolfi.giorgia/demo/drivers.txt",header = T,sep = "\t") %>% 
  #   mutate(chr=paste0("chr",chr))
  # gene_id <- inner_join(x = drivers,y = drivers_table,by=c("chr","from","ref","alt")) %>% 
  #   select(chr,from,driver_gene) %>% 
  #   unique()
  # gr3 = GRanges(seqnames = gene_id[, 1], ranges = IRanges(gene_id[, 2], gene_id[, 2]))
  # gr3$gene = gene_id$driver_gene
  # 
  # mtch = as.matrix(findOverlaps(chr_window, gr3))
  # at = mtch[, 1]
  # labels = mcols(gr3)[mtch[, 2], 1]
  
  ht <-Heatmap(t(char_mat), name = "CNV", 
               col=my_palette,
               # col=col[names(col)%in%karyotypes],
               top_annotation = column_ha,
               right_annotation = row_ha,
               column_split = chr, cluster_columns = FALSE, show_row_dend = TRUE,
               column_title_gp = gpar(fontsize = 10), 
               border = TRUE,
               column_gap = unit(0, "points"),
               column_title = ifelse(1:22 %% 2 == 0, paste0("\n", chr_level), paste0(chr_level, "\n")),
               heatmap_legend_param = list(direction = "horizontal", title_position = "lefttop"))
  pdf(output_file,width = 20,height = 20)
  draw(ht)
  dev.off()
  return(ht)
}
sample_forest <- load_sample_forest("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-DLP/dlp_simulation_7/sample_forest_7.sff")
phylo_forest <- load_phylogenetic_forest("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-DLP/dlp_simulation_7/phylo_forest_7.sff")
ht_process <- plot_ProCESS_heatmap(sample_forest = forest,phylo_forest = phylo_forest, output_file="/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-DLP/dlp_simulation_7/cna_heatmap.pdf")
# all_nodes <- forest$get_nodes() %>%
#   filter(!is.na(sample))
# sampled_mutants <- all_nodes %>%
#   separate(col = sample,into = c("p1","s1","c1"),sep = "_") %>%
#   mutate(sample_name=paste(p1,s1,sep="_")) %>%
#   group_by(sample_name,mutant) %>%
#   summarise(prop=n())
# 
# ggplot(sampled_mutants, aes(x = "", y = prop, fill = mutant)) +
#   geom_bar(stat = "identity", width = 1) +
#   scale_fill_manual(values=c("A"="darkseagreen4","B"="lightpink2","C"="cadetblue3"))+
#   coord_polar("y") +
#   facet_wrap(~ sample_name,scales = "free",nrow = 1) +
#   theme_void() +
#   labs(title = "Proportion of Mutants per Sample") +
#   theme(
#     plot.title = element_text(hjust = 0.5),
#     legend.title = element_blank()
#   )
