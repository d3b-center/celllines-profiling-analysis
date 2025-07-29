# function for assessing batch corrected data 
# PCA plot colored by batch and a separate PCA colored by molecular_subtype
suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
  library(ggpubr)
  library(gridExtra)
})

plot_pca <- function(expr_mat, log_normalized = T, hist_df, prefix, plots_dir){
  set.seed(100)
  
  if(log_normalized){
    # do nothing
  } else {
    expr_mat <- log2(expr_mat + 1)
  }
  prData <- prcomp(expr_mat)
  var_explained = round(prData$sdev^2/sum(prData$sdev^2)*100, 1)
  pca_data <- prData$rotation
  pca_data <- data.frame(pca_data)[1:4]
  
  # merge with annotations
  pca_data <- pca_data %>%
    rownames_to_column("Kids_First_Biospecimen_ID") %>%
    inner_join(hist_df, by = "Kids_First_Biospecimen_ID")
  
  # color by cell line composition and shape by RNA library
  p1 <- ggplot(pca_data, aes(PC1, PC2)) +
    geom_point(size = 5, alpha = 0.5, aes(color = as.character(cell_line_composition), shape = as.character(RNA_library))) +
    theme_bw() +
    ggtitle("PCA Clustering") +
    ggpubr::theme_pubr(legend = "right") + 
    xlab(paste0('PC1 (', var_explained[1], '%)')) + 
    ylab(paste0('PC2 (', var_explained[2], '%)')) 
  ggsave(filename = file.path(plots_dir, paste0(prefix, "_pca.pdf")), plot = p1, width = 8, height = 8)
}
