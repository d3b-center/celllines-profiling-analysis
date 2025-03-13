# Author: Komal S. Rathi
# Date: 02/24/2020
# GSVA using proliferative genes

library(GSVA)
library(ggpubr)
source('R/pubTheme.R')

# z-score
zscore <- function(x){
  z <- (x - mean(x)) / sd(x)
  return(z)
}

# z-score FPKM data
expr.fpkm <- get(load('data/fpkm-matrix.RData'))
expr.fpkm <- log2(expr.fpkm+1)
expr.fpkm <- apply(expr.fpkm, 1, zscore)
expr.fpkm <- t(expr.fpkm)

# GSVA
sig <- read.delim('data/cell-cycle-genes.txt', header = F, stringsAsFactors = F)
sig <- list(cell_cycle = sig$V1)
GeneSetExprsMat <- gsva(expr.fpkm, gset.idx.list = sig,
                        abs.ranking = F, min.sz = 1,
                        max.sz = 1500, parallel.sz = 4, mx.diff = F,
                        parallel.type = 'SOCK')
GeneSetExprsMat <- melt(GeneSetExprsMat)
GeneSetExprsMat <- cbind(GeneSetExprsMat, type = gsub(".*_","",GeneSetExprsMat$Var2))
GeneSetExprsMat$type  <- ifelse(GeneSetExprsMat$type  == "tissue", "Solid_Tissue", 
       ifelse(GeneSetExprsMat$type  == "s", "Suspension", "Adherent"))
GeneSetExprsMat$type <- factor(GeneSetExprsMat$type, levels = c("Adherent", "Solid_Tissue", "Suspension"))

# plot
my_comparisons <- list(c("Suspension", "Adherent"), c("Suspension", "Solid_Tissue"), c("Adherent", "Solid_Tissue"))
p <- ggplot(GeneSetExprsMat, aes(x = type, y = value, fill = type)) + 
  stat_boxplot(geom ='errorbar', width = 0.2) +
  geom_boxplot(outlier.shape = 21, outlier.fill = "white", outlier.color = "white",
               lwd = 0.5, fatten = 0.7, width = 0.5) +
  geom_jitter(position=position_jitter(width=.25), shape = 21) +
  theme_Publication() + guides(fill = F) +
  ggtitle("GSVA: Cell cycle genes") + 
  ylab("Proliferative Index") + xlab("Type") +
  stat_compare_means(comparisons = my_comparisons, color = "darkred", size = 3) +
  stat_compare_means(color = "darkred", label.y = 1)
p
ggsave(filename = "results/proliferative-index.pdf", plot = p, device = "pdf", width = 8, height = 6)
