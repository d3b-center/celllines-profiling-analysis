# Author: Komal S. Rathi
# Date: 11/26/2019
# Function: Clustering of RNA-seq data 

library(Rtsne)
library(reshape2)
library(ggplot2)

# load libraries and source plotting theme
load('data/input.RData')
source('R/pubTheme.R')

# format data for t-SNE
df <- data.frame(snames = colnames(expr.collapsed), colsplit(colnames(expr.collapsed), pattern = "_", names = c("sample","type")))
rownames(df) <- df$snames
rownames(expr.collapsed) <- NULL
expr.collapsed <- unique(expr.collapsed)

# t-SNE and plot
set.seed(42)
tsneOut <- Rtsne(t(expr.collapsed), initial_dims = 50, perplexity = 8, max_iter = 1000)
tsneOut <- data.frame(tsneOut$Y, df)
p <- ggplot(tsneOut, aes(X1, X2, color = type)) +
  geom_point(size = 3, alpha = 0.5) +
  theme_bw() +
  ggtitle("T-SNE of HGG cell lines") +
  theme_Publication2() + xlab("PC1") + ylab("PC2") +
  ggrepel::geom_text_repel(aes(label = sample)) +
  scale_color_discrete(name = "Type", labels = c("a" = "Adhesion", 
                                  "s" = "Suspension", 
                                  "solid_tissue" = "Solid Tissue"))
p
ggsave(filename = 'results/tSNE_plot.png', plot = p, width = 8, height = 6)

# PCA
pcaOut <- prcomp(t(expr.collapsed), center = TRUE,scale. = TRUE)
pcaOut <- cbind(pcaOut$x[,c("PC1","PC2")], df)
p <- ggplot(pcaOut, aes(PC1, PC2, color = type)) +
  geom_point(size = 3, alpha = 0.5) +
  theme_bw() +
  ggtitle("PCA of HGG cell lines") +
  theme_Publication2() + xlab("PC1") + ylab("PC2") +
  ggrepel::geom_text_repel(aes(label = sample)) +
  scale_color_discrete(name = "Type", labels = c("a" = "Adhesion", 
                                                 "s" = "Suspension", 
                                                 "solid_tissue" = "Solid Tissue"))
p
ggsave(filename = 'results/PCA_plot.png', plot = p, width = 8, height = 6)
