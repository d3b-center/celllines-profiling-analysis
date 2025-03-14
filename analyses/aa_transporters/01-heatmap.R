suppressPackageStartupMessages({
  library(tidyverse)
  library(DESeq2)
  library(pheatmap)
  library(NOISeq)
})

# set directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
data_dir <- file.path(root_dir, "data")
analysis_dir <- file.path(root_dir, "analyses", "aa_transporters")
plots_dir <- file.path(analysis_dir, "plots")
dir.create(plots_dir, showWarnings = FALSE, recursive = TRUE)

# source functions
source(file.path(analysis_dir, "utils", "plot_umap.R"))
source(file.path(analysis_dir, "utils", "plot_pca.R"))

# read count data
count_dat <- readRDS(file.path(data_dir, "20250115", "20250115_release.gene-counts-rsem-expected_count-collapsed.all.rds"))

# read histologies
hist_df <- read_tsv(file.path(data_dir, "20250115", "20250115_release.annotated_histologies.tsv"))

# subset to HGAT cell line data
hist_df <- hist_df %>%
  filter(short_histology == "HGAT",
         sample_type == "Tumor",
         composition == "Derived Cell Line",
         Kids_First_Biospecimen_ID %in% colnames(count_dat)) %>%
  dplyr::select(Kids_First_Biospecimen_ID, cell_line_composition, tumor_descriptor, molecular_subtype, CNS_region, RNA_library)
count_dat <- count_dat %>%
  dplyr::select(any_of(hist_df$Kids_First_Biospecimen_ID))

# subset to protein coding genes only
gencode_gtf <- rtracklayer::import(con = file.path(data_dir, "gencode.v39.primary_assembly.annotation.gtf.gz"))
gencode_gtf <- as.data.frame(gencode_gtf)
gencode_gtf <- gencode_gtf %>%
  dplyr::select(gene_id, gene_name, gene_type) %>%
  dplyr::filter(gene_type == "protein_coding") %>%
  unique()
count_dat <- count_dat %>%
  filter(rownames(count_dat) %in% gencode_gtf$gene_name)

# remove low expression genes
count_dat <- DGCA::filterGenes(
  inputMat = count_dat,
  filterTypes = c("central", "dispersion"),
  filterDispersionType = "cv",
  filterDispersionPercentile = 0.2,
  sequential = TRUE)

# DESeq2 dds object
dds <- DESeq2::DESeqDataSetFromMatrix(
  countData = round(count_dat),
  colData = hist_df,
  design = ~ cell_line_composition
)
dds <- DESeq2::DESeq(dds)
counts_norm <- DESeq2::counts(dds, normalized = T) # normalized counts

# UMAP: log transformed counts
plot_umap(expr_mat = counts_norm,
          log_normalized = F,
          hist_df,
          prefix = "dds_norm_data",
          plots_dir = plots_dir)

# PCA: log transformed counts
plot_pca(expr_mat = counts_norm,
         log_normalized = F,
         hist_df,
         prefix = "dds_norm_data",
         plots_dir = plots_dir)

# If batch correction is needed, run NOISeq's Arsynseq,
# specifying factor = RNA_library and batch = TRUE and
hist_df <- hist_df %>%
  column_to_rownames("Kids_First_Biospecimen_ID")
noiseq_data <- readData(
  data = as.data.frame(counts_norm),
  factors = hist_df %>%
    mutate("batch" = "RNA_library") # define batch
)
# phenoData rownames should match column names of matrix
stopifnot(identical(colnames(noiseq_data), rownames(hist_df)))

# add a small value to expression data before running ARSyNseq
exprs(noiseq_data) <- exprs(noiseq_data) + 0.01
noiseq_batch_true <- NOISeq::ARSyNseq(
  data = noiseq_data,
  factor = "RNA_library",
  batch = TRUE,
  norm = "uqua",
  logtransf = T
)
hist_df <- hist_df %>%
  rownames_to_column("Kids_First_Biospecimen_ID")

# UMAP: after NOISeq + Batch
plot_umap(
  expr_mat = exprs(noiseq_batch_true),
  log_normalized = T,
  hist_df,
  prefix = "noiseq_batch_true",
  plots_dir = plots_dir
)

# PCA: after NOISeq + Batch
plot_pca(
  expr_mat = exprs(noiseq_batch_true),
  log_normalized = T,
  hist_df,
  prefix = "noiseq_batch_true",
  plots_dir = plots_dir
)

# # Run Deseq2's variance stabilizing transformation (vst) on the expression data
# count_dat <- DESeq2::varianceStabilizingTransformation(round(as.matrix(count_dat)), fitType = "parametric", blind = T)
# count_dat <- as.data.frame(count_dat)

# use noiseq batch corrected data for heatmap
noiseq_dat <- exprs(noiseq_batch_true)

# filter to genes of interest
reactome_dat <- msigdbr::msigdbr(category = "C2", subcategory = "REACTOME")
genes_of_interest <- reactome_dat %>% 
  filter(gs_exact_source == "R-HSA-352230") %>%
  pull(gene_symbol) %>%
  unique()
print(length(genes_of_interest)) # 33 genes in REACTOME_AMINO_ACID_TRANSPORT_ACROSS_THE_PLASMA_MEMBRANE
genes_of_interest <- unique(c(genes_of_interest, "ACLY", "ACSS1", "ACSS2", "ACACA", "FASN")) # 38 genes total

# 28/38 genes are expressed
noiseq_dat <- noiseq_dat %>%
  as.data.frame() %>%
  filter(rownames(noiseq_dat) %in% genes_of_interest)

# annotation colors
colors <- rlist::list.load(file.path(data_dir, 'colors.yaml'))
mycolors <- lapply(colors, function(x) unlist(x))

# set NA to N/A
hist_df[is.na(hist_df)] <- "N/A"

# generate heatmap
matrix_colors = colorRampPalette(c("cadetblue3", "white", "coral3"))(25)
pdf(file = file.path(plots_dir, "hgat_cellline_heatmap.pdf"), width = 14)
pheatmap(noiseq_dat, 
         scale = "row", 
         color = matrix_colors,
         annotation_col = hist_df %>%
           column_to_rownames("Kids_First_Biospecimen_ID") %>%
           dplyr::select(molecular_subtype, CNS_region, tumor_descriptor, cell_line_composition),
         angle_col = 45, 
         annotation_colors = mycolors, 
         main = paste0("HGAT Derived Cell Lines (n = ", ncol(noiseq_dat), ")", "\nVST normalized Counts"))
dev.off()
