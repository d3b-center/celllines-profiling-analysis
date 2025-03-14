# function for assessing batch corrected data 
# UMAP plot colored by batch and a separate UMAP colored by molecular_subtype
suppressPackageStartupMessages({
  library(tidyverse)
  library(uwot)
  library(ggplot2)
  library(ggpubr)
  library(gridExtra)
})

plot_umap <- function(expr_mat, log_normalized = T, hist_df, prefix, plots_dir){
  set.seed(100)
  
  if(log_normalized){
    expr_mat <- t(expr_mat)
  } else {
    expr_mat <- log2(t(expr_mat)+1)
  }
  umap_out <- uwot::umap(X = expr_mat, n_neighbors = 15, n_components = 2, metric = "correlation", ret_nn = TRUE, n_sgd_threads = 1)
  umap_out <- umap_out$embedding
  
  # add colnames/rownames to embeddings
  colnames(umap_out) <- c("UMAP1", "UMAP2")
  rownames(umap_out) <- hist_df$Kids_First_Biospecimen_ID
  
  # merge with annotations
  umap_out <- umap_out %>%
    as.data.frame() %>%
    rownames_to_column("Kids_First_Biospecimen_ID") %>%
    inner_join(hist_df, by = "Kids_First_Biospecimen_ID")
  
  # color by cell line composition and shape by RNA library
  p1 <- ggplot(umap_out, aes(UMAP1, UMAP2)) +
    geom_point(size = 5, alpha = 0.5, aes(color = as.character(cell_line_composition), shape = as.character(RNA_library))) +
    theme_bw() +
    ggtitle("UMAP Clustering") +
    ggpubr::theme_pubr(legend = "right") + 
    xlab("X1") + 
    ylab("X2") 
  ggsave(filename = file.path(plots_dir, paste0(prefix, "_umap.pdf")), plot = p1, width = 8, height = 8)
}

