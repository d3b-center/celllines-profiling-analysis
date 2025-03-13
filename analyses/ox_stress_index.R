# Author: Komal S. Rathi
# Date: 02/24/2020
# GSVA using oxidative stress genes
library(msigdbr)
library(tidyverse)
library(GSVA)
library(ggpubr)
library(reshape2)
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

# get emt gene set
hallmark_geneset <- msigdbr::msigdbr(species = "Homo sapiens", category = "C5", subcategory = "BP")
oxidative_stress_geneset <- hallmark_geneset %>% 
  filter(gs_name %in% c("GO_RESPONSE_TO_OXIDATIVE_STRESS",
                        "GO_RESPONSE_TO_HYDROGEN_PEROXIDE",
                        "GO_CELLULAR_RESPONSE_TO_HYDROGEN_PEROXIDE"))

run_ssgsea <- function(expr, geneset){
  geneset <- list(geneset = geneset$human_gene_symbol)
  GeneSetExprsMat <- gsva(expr = expr, 
                          gset.idx.list = geneset,
                          method = "ssgsea",
                          min.sz = 1, max.sz = 1500, 
                          mx.diff = F)
}

GeneSetExprsMat <- plyr::ddply(.data = oxidative_stress_geneset, 
                               .variables = 'gs_name', 
                               .fun = function(x) run_ssgsea(expr = expr.fpkm, geneset = x))
GeneSetExprsMat <- melt(GeneSetExprsMat)
GeneSetExprsMat <- cbind(GeneSetExprsMat, type = gsub(".*_","",GeneSetExprsMat$variable))
GeneSetExprsMat$type  <- ifelse(GeneSetExprsMat$type  == "tissue", "Solid_Tissue", 
                                ifelse(GeneSetExprsMat$type  == "s", "Suspension", "Adherent"))
GeneSetExprsMat$type <- factor(GeneSetExprsMat$type, levels = c("Adherent", "Solid_Tissue", "Suspension"))

# plot
my_comparisons <- list(c("Suspension", "Adherent"), 
                       c("Suspension", "Solid_Tissue"), 
                       c("Adherent", "Solid_Tissue"))

p <- ggplot(GeneSetExprsMat, aes(x = type, y = value, fill = factor(gs_name, unique(as.character(gs_name))))) + 
  stat_boxplot(geom ='errorbar', width = 0.7, lwd = 0.3) +
  geom_boxplot(outlier.shape = 21, outlier.fill = "white", outlier.color = "white",
               lwd = 0.3, fatten = 0.7, width = 0.7) +
  theme_Publication2() + 
  ggtitle("GSVA: Oxidative stress index") + 
  ylab("Oxidative stress index") + xlab("") +
  stat_compare_means(comparisons = my_comparisons, color = "darkred", size = 3) +
  stat_compare_means(color = "darkred", label.y = 1.6) +
  labs(fill = "GO genesets")
p
ggsave(filename = "results/oxidative-stress-index.pdf", plot = p, device = "pdf", width = 12, height = 6)
