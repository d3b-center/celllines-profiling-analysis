suppressPackageStartupMessages({
  library(tidyverse)
  library(DESeq2)
  library(pheatmap)
})

# set directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
data_dir <- file.path(root_dir, "data")
analysis_dir <- file.path(root_dir, "analyses", "aa_transporters")
output_dir <- file.path(analysis_dir, "results")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

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
  dplyr::select(Kids_First_Biospecimen_ID, tumor_descriptor, molecular_subtype, CNS_region)
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

# # remove low expression genes
# count_dat <- DGCA::filterGenes(
#   inputMat = count_dat,
#   filterTypes = c("central", "dispersion"),
#   filterDispersionType = "cv",
#   filterDispersionPercentile = 0.2,
#   sequential = TRUE)

# Run Deseq2's variance stabilizing transformation (vst) on the expression data
count_dat <- DESeq2::varianceStabilizingTransformation(round(as.matrix(count_dat)), fitType = "parametric", blind = T)
count_dat <- as.data.frame(count_dat)

# filter to genes of interest
reactome_dat <- msigdbr::msigdbr(category = "C2", subcategory = "REACTOME")
genes_of_interest <- reactome_dat %>% 
  filter(gs_exact_source == "R-HSA-352230") %>%
  pull(gene_symbol) %>%
  unique()
print(length(genes_of_interest)) # 33 genes in REACTOME_AMINO_ACID_TRANSPORT_ACROSS_THE_PLASMA_MEMBRANE
genes_of_interest <- unique(c(genes_of_interest, "ACLY", "ACSS1", "ACSS2", "ACACA", "FASN")) # 38 genes total

# all genes are used
count_dat <- count_dat %>%
  filter(rownames(count_dat) %in% genes_of_interest)

# annotation colors
colors <- rlist::list.load(file.path(data_dir, 'colors.yaml'))
mycolors <- lapply(colors, function(x) unlist(x))

# set NA to N/A
hist_df[is.na(hist_df)] <- "N/A"

# without filtering genes perform manual scaling as scale = "row" does not work in the pheatmap function due to zero standard deviation issue
count_dat_scaled <- t(count_dat) %>%
  as.data.frame() %>%
  mutate_at(c(rownames(count_dat)), ~ (scale(.) %>% as.vector))
count_dat_scaled <- t(count_dat)
count_dat_scaled <- count_dat_scaled %>%
  as.data.frame() %>%
  mutate_at(c(colnames(count_dat_scaled)), ~ (scale(.) %>% as.vector))

# generate heatmap
matrix_colors = colorRampPalette(c("cadetblue3", "white", "coral3"))(25)
pdf(file = file.path(output_dir, "hgat_cellline_heatmap_all.pdf"), width = 14)
pheatmap(t(count_dat_scaled),
         color = matrix_colors,
         annotation_col = hist_df %>%
           column_to_rownames("Kids_First_Biospecimen_ID") %>%
           dplyr::select(molecular_subtype, CNS_region, tumor_descriptor),
         angle_col = 45,
         annotation_colors = mycolors,
         main = paste0("HGAT Derived Cell Lines (n = ", ncol(count_dat), ")", "\nVST normalized Counts"))
dev.off()
