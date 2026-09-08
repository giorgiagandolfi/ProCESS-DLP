## ---------------------------------------------------------------
## Cell x genomic-position copy-number heatmap from segment data
## Input columns expected: chr,start,end,state,median,multiplier,cell_id
## ---------------------------------------------------------------

library(data.table)
library(GenomicRanges)
library(ComplexHeatmap)
library(circlize)
setwd("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-DLP/results_mondrian_qc_spn01")
# ================= PARAMETERS (edit these) =================
input_file <- "SPN01_hmmcopy_segments.csv.gz"     # path to your segments CSV
output_pdf <- "cnv_heatmap_hmmcopy.pdf"
output_png <- "cnv_heatmap.png"

bin_size   <- 5e5                # genomic bin width in bp (500kb default)
chr_order  <- paste0("chr", c(1:22, "X", "Y"))   # adjust if you have other contigs
value_col  <- "state"            # column to plot: "state" (integer CN) or "median"
# cells_to_use <- readRDS("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-DLP/SCOUT/SPN01/process/cells_subclonal_cn.rds")
# # ==============================================================
# # cells_to_use_subsampled <- sample(x = cells_to_use,size = 100)
# cells_to_use_notsub <- readRDS("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-DLP/SCOUT/SPN01/process/cells_not_subclonal.rds")
cells_metadata <- readRDS("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-DLP/SCOUT/SPN01/process/sampled_cells_spn01_info.rds")
# cells_to_use_all <- c(cells_to_use_notsub,cells_to_use)
# cells_metadata_filtered <- cells_metadata %>% 
#   filter(sample%in%cells_to_use_all)
# ---- Read segments ----
# cells_to_use <- readRDS("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-DLP/SCOUT/SPN01/process/sampled_cells_balanced.rds")
# cells_to_use <- cells_to_use %>% 
#   mutate(sample=paste(patient,sample_prefix,cell_pos,sep="_")) %>% 
#   filter(sample_prefix=="1.2")

segs <- fread(input_file)
segs <- segs %>% 
  filter(startsWith(cell_id, "SPN01_1.2"))
correct_called_cells <- segs %>% 
  filter(chr=="chr5") %>% filter(start==107500001) %>% 
  filter(state==1) %>% 
  pull(cell_id)
segs <- segs %>% 
  filter(cell_id%in%correct_called_cells)
segs[, chr := factor(chr, levels = chr_order)]
n_dropped <- sum(is.na(segs$chr))
if (n_dropped > 0) {
  message(sprintf("Dropping %d rows with chromosomes not in chr_order", n_dropped))
}
segs <- segs[!is.na(chr)]

# ---- Build genome-wide bins (per-chromosome, sized from data) ----
chr_max <- segs[, .(max_end = max(end)), by = chr]
setorder(chr_max, chr)

bins_list <- lapply(seq_len(nrow(chr_max)), function(i) {
  ch  <- chr_max$chr[i]
  len <- chr_max$max_end[i]
  starts <- seq(1, len, by = bin_size)
  data.table(chr = ch, start = starts, end = pmin(starts + bin_size - 1, len))
})
bins <- rbindlist(bins_list)
setorder(bins, chr, start)
bins[, bin_id := paste0(chr, ":", formatC(start, width = 12, flag = "0"))]
bins[, bin_id := factor(bin_id, levels = bin_id)]   # preserves genome order

# ---- Overlap segments onto bins with GenomicRanges ----
seg_gr <- GRanges(
  seqnames = segs$chr,
  ranges   = IRanges(segs$start, segs$end),
  cell_id  = segs$cell_id,
  value    = segs[[value_col]]
)
bin_gr <- GRanges(
  seqnames = bins$chr,
  ranges   = IRanges(bins$start, bins$end),
  bin_id   = bins$bin_id
)

hits <- findOverlaps(bin_gr, seg_gr)
overlap_dt <- data.table(
  bin_id  = bin_gr$bin_id[queryHits(hits)],
  cell_id = seg_gr$cell_id[subjectHits(hits)],
  value   = seg_gr$value[subjectHits(hits)]
)

# A bin can straddle a segment breakpoint and overlap two segments; keep the
# first hit (by widest overlap) so each (bin, cell) pair has a single value.
overlap_dt <- unique(overlap_dt, by = c("bin_id", "cell_id"))

# ---- Cast to matrix: rows = cell_id, columns = bin_id (genome order) ----
mat_dt <- dcast(overlap_dt, cell_id ~ bin_id, value.var = "value")
cell_ids <- mat_dt$cell_id
mat_dt[, cell_id := NULL]
mat <- as.matrix(mat_dt)
rownames(mat) <- cell_ids

# reorder columns to genome order, keep only bins actually present
ordered_bins <- levels(bins$bin_id)[levels(bins$bin_id) %in% colnames(mat)]
mat <- mat[, ordered_bins, drop = FALSE]

# ---- Chromosome split for column annotation / gaps ----
col_chr <- bins$chr[match(colnames(mat), bins$bin_id)]
col_chr <- droplevels(col_chr)

# ---- Copy-number color scale ----
# Adjust breakpoints to match your ploidy expectations; this centers on CN=2.
max_state <- suppressWarnings(max(mat, na.rm = TRUE))
if (!is.finite(max_state)) max_state <- 6
cn_colors <-  c("0"="#3590d2",
"1"="#7fbbe3",
"2" = "#d4cfcd",
"3" = "#fdc38dff",
"4" = "#fca16cff",
"5" = "#f67b51ff",
"6" = "#e7533aff",
"7" = "#cf2518ff",
"8" = "#ad0000ff",
"9" = "#c4198aff",
"10" = '#990071ff',
"11"='#78056cff')
cells_metadata_filtered = cells_metadata %>% 
  filter(startsWith(sample, "SPN01_1.2")) %>% 
  filter(sample%in%correct_called_cells)
row_ha = rowAnnotation(clone = cells_metadata_filtered$mutant,
                       # bulk_sample = cells_metadata_filtered$,
                       col = list(clone = c("Clone 3" = "#c7c6c4ff", "Clone 4" = "#9b9b9bff")))
                                  # bulk_sample = c("SPN01_1.1"="goldenrod","SPN01_1.2"="#5c3c88ff","SPN01_1.3"="magenta4")))

# ---- Draw heatmap ----
pdf(output_pdf, width = 14, height = 10)
# png(output_png, width = 14, height = 10)
ht <- Heatmap(
  mat,
  name = if (value_col == "state") "Copy\nnumber" else value_col,
  col = cn_colors,
  cluster_rows = TRUE,               # cluster cells by CN profile similarity
  cluster_columns = FALSE,           # keep genomic order within each chr
  cluster_column_slices = FALSE,
  column_split = col_chr,            # one block per chromosome
  column_gap = unit(0.5, "mm"),
  column_title_gp = gpar(fontsize = 8),
  column_title_rot = 90,
  show_column_names = FALSE,
  show_row_names = nrow(mat) <= 100,
  row_names_gp = gpar(fontsize = 6),
  right_annotation=row_ha,
  na_col = "white",
  border = TRUE,
  use_raster = TRUE
)
draw(ht)
dev.off()