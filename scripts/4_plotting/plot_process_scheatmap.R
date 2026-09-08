## ---------------------------------------------------------------
## Cell x genomic-position TOTAL copy-number heatmap from
## allele-specific segment data (total CN = major + minor).
##
## Input columns expected: cell_id,chr,begin,end,major,minor,sample
## ---------------------------------------------------------------

library(data.table)
library(GenomicRanges)
library(ComplexHeatmap)
library(circlize)

# ================= PARAMETERS (edit these) =================
input_file <- "/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-DLP/SCOUT/SPN01/process/single_cell_cna.rds"
output_pdf <- "cnv_heatmap.pdf"
bin_size   <- 5e5                # genomic bin width in bp (500kb default)
cell_col   <- "cell_id"          # column identifying each cell (or use "sample")
# ==============================================================

# ---- Read segments ----
segs <- readRDS(input_file)
setDT(segs)   # ensure data.table class even if segs came from elsewhere (e.g. read.csv)
setnames(segs, cell_col, "cell_id")

cells_metadata <- readRDS("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-DLP/SCOUT/SPN01/process/sampled_cells_spn01_info.rds")

# cells_to_use <- readRDS("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-DLP/SCOUT/SPN01/process/sampled_cells_balanced.rds")
# cells_to_use <- cells_to_use %>% 
#   mutate(sample=paste(patient,sample_prefix,cell_pos,sep="_")) %>% 
#   filter(sample_prefix=="1.2")
# ==============================================================
# cells_to_use_subsampled <- sample(x = cells_to_use,size = 100)
# cells_to_use_notsub <- readRDS("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-DLP/SCOUT/SPN01/process/cells_not_subclonal.rds")
# cells_metadata <- readRDS("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-DLP/SCOUT/SPN01/process/sampled_cells_spn01_info.rds")
# cells_to_use_all <- c(cells_to_use_notsub,cells_to_use)
# cells_metadata_filtered <- cells_metadata %>% 
#   filter(sample%in%cells_to_use_all)

segs <- segs %>% 
  filter(startsWith(sample, "SPN01_1.2"))
  # filter(sample%in%cells_to_use$sample)
# chr may come as bare integers/strings ("1","2",...,"X","Y") -- normalize
segs[, chr := as.character(chr)]
chr_levels <- unique(segs$chr)
num_chrs <- chr_levels[!is.na(suppressWarnings(as.integer(chr_levels)))]
num_chrs <- as.character(sort(as.integer(num_chrs)))
other_chrs <- sort(setdiff(chr_levels, num_chrs))  # e.g. "X","Y"
chr_order <- c(num_chrs, other_chrs)
segs[, chr := factor(chr, levels = chr_order)]

n_dropped <- sum(is.na(segs$chr))
if (n_dropped > 0) message(sprintf("Dropping %d rows with unrecognized chr", n_dropped))
segs <- segs[!is.na(chr)]

# ---- Total copy number = major + minor ----
segs[, value := major + minor]

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
  ranges   = IRanges(segs$begin, segs$end),
  cell_id  = segs$cell_id,
  value    = segs$value
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
overlap_dt <- unique(overlap_dt, by = c("bin_id", "cell_id"))  # resolve breakpoint dupes

# ---- Cast to matrix: rows = cell_id, columns = bin_id (genome order) ----
mat_dt <- dcast(overlap_dt, cell_id ~ bin_id, value.var = "value")
cell_ids <- mat_dt$cell_id
mat_dt[, cell_id := NULL]
mat <- as.matrix(mat_dt)
rownames(mat) <- cell_ids

ordered_bins <- levels(bins$bin_id)[levels(bins$bin_id) %in% colnames(mat)]
mat <- mat[, ordered_bins, drop = FALSE]

# ---- Chromosome split for column blocks ----
col_chr <- bins$chr[match(colnames(mat), bins$bin_id)]
col_chr <- droplevels(col_chr)

# ---- Copy-number color scale (centered on diploid = 2) ----
max_val <- suppressWarnings(max(mat, na.rm = TRUE))
if (!is.finite(max_val)) max_val <- 6
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
# ---- Draw heatmap ----
# cells_metadata_filtered <- cells_to_use %>% 
#   filter(sample%in%cells_to_use_all)
cells_metadata_filtered = cells_metadata %>% 
  filter(startsWith(sample, "SPN01_1.2"))
row_ha = rowAnnotation(clone = cells_metadata_filtered$mutant,
                       # bulk_sample = cells_metadata_filtered$,
                       col = list(clone = c("Clone 3" = "#c7c6c4ff", "Clone 4" = "#9b9b9bff")))
# bulk_sample = c("SPN01_1.1"="goldenrod","SPN01_1.2"="#5c3c88ff","SPN01_1.3"="magenta4")))pdf(output_pdf, width = 14, height = 10)
pdf(output_pdf, width = 14, height = 10)
ht <- Heatmap(
  mat,
  name = "Total\nCN",
  col = cn_colors,
  cluster_rows = TRUE,               # cluster cells by CN profile similarity
  cluster_columns = FALSE,
  cluster_column_slices = FALSE,
  right_annotation = row_ha,
  column_split = col_chr,            # one block per chromosome
  column_gap = unit(0.5, "mm"),
  column_title_gp = gpar(fontsize = 8),
  column_title_rot = 90,
  show_column_names = FALSE,
  show_row_names = nrow(mat) <= 100,
  row_names_gp = gpar(fontsize = 6),
  na_col = "white",
  border = TRUE,
  use_raster = TRUE
)
draw(ht)
dev.off()

message(sprintf("Wrote heatmap for %d cells x %d bins to %s", nrow(mat), ncol(mat), output_pdf))