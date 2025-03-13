# Author: Komal S. Rathi
# Date: 11/26/2019
# Function: Brain cell type proportion analysis using BRETIGEA 

# load libraries
library(BRETIGEA, quietly = TRUE)
library(reshape2)
library(ggplot2)
library(ggpubr)

# load data
source('R/pubTheme.R')
load('data/fpkm-matrix.RData') # FPKM

# create function
plot.bretigea <- function(expr, title){
  ct_res <- brainCells(expr, nMarker = 50)
  ct_res <- cbind(ct_res, colsplit(rownames(ct_res), pattern = "_", names = c("sample", "type")))
  ct_res <- melt(ct_res)

  # change factor for type
  ct_res$type <- factor(ct_res$type, levels = c("s", "a", "solid_tissue"))
  levels(ct_res$type) <- c("Suspension\n(7)", "Adhesion\n(9)", "Solid Tissue\n(12)")
  
  # change factor for variable
  levels(ct_res$variable) <- c("Astrocytes", "Endothelial Cells", "Microglia",
                               "Neurons", "Oligodendrocytes", "Oligodendrocyte PC") 
  
  #ct_res$variable <- : astrocytes (ast), endothelial cells
  #(end), microglia (mic), neurons (neu), oligodendrocytes (oli), and oligodendrocyte precursor cells (opc)
  my_comparisons <- list(c(levels(ct_res$type)[1], levels(ct_res$type)[2]),
                         c(levels(ct_res$type)[1], levels(ct_res$type)[3]),
                         c(levels(ct_res$type)[2], levels(ct_res$type)[3]))
  
  p <- ggplot(ct_res, aes(x = type, y = value, fill = type)) + 
    stat_boxplot(geom ='errorbar', width = 0.2) +
    geom_boxplot(lwd = 0.5, fatten = 0.7, outlier.shape = 1, width = 0.5, outlier.size = 1) +
    facet_wrap(~variable) +
    theme_Publication() + xlab("") + ylab("BRETIGEA score") +
    guides(fill = FALSE) + ggtitle(title) +
    stat_compare_means(comparisons = my_comparisons, color = 'red', size = 3)
  return(p)
}

# use FPKM
p <- plot.bretigea(expr = expr.collapsed, title = "Brain cell type proportions") # fpkm
ggsave(filename = 'results/bretigea_plot.png', plot = p,width = 10, height = 7)
