# Author: Komal S. Rathi
# Date: 02/17/2019
# Function: Volcano plots of TMM proteins

library(dplyr)
library(ggplot2)
library(scales)

source('R/pubTheme.R')

# read DEGs
s.vs.tissue <- read.delim('results/limma/s_vs_tissue_limma.txt', stringsAsFactors = F)
a.vs.tissue <- read.delim('results/limma/a_vs_tissue_limma.txt', stringsAsFactors = F)

# common genes
common <- merge(s.vs.tissue, a.vs.tissue, by = "Gene")
common <- common %>%
  filter(adj.P.Val.x < 0.05 & adj.P.Val.y < 0.05) %>%
  filter((logFC.x > 0 & logFC.y > 0) | 
           logFC.x < 0 & logFC.y < 0) %>%
  .$Gene %>%
  unique()

# surface proteins
load('~/Projects/marislab-webportal/data/TMlist.RData')

# create volcano plot
plotVolcano <- function(x, TM = TRUE, title = "Volcano Plot", otherCol = c("blue","gray","red"), pval.cutoff = 0.05, fc.cutoff =  1, goi) {
  
  # filter by pvalue and logFC
  x <- x %>%
    mutate(sig = ifelse(adj.P.Val < pval.cutoff, T, F),
           DEGAnnot = ifelse(adj.P.Val < pval.cutoff & logFC > 1*fc.cutoff, 
                             yes = "Up", 
                             no = ifelse(test = adj.P.Val < pval.cutoff & logFC < -1*fc.cutoff, 
                                         yes = "Down", 
                                         no = "Unchanged")),
           tmm = ifelse(Gene %in% TMlist$hgnc_symbol, T, F))
  
  # filter to TMM
  if(TM == TRUE){
    x <- x %>%
      filter(tmm == TRUE)
  }
  
  # count number of each type and add as labels
  nums <- plyr::count(x$DEGAnnot)
  x <- merge(x, nums, by.x = 'DEGAnnot', by.y = 'x')
  x$DEGAnnot <- paste0(x$DEGAnnot," (n = ", x$freq, ")")
  
  # plot volcano
  p <- ggplot(x, aes(x = logFC, y = -log10(adj.P.Val), color = DEGAnnot)) +
    geom_point(alpha = 0.5) +  
    ggtitle("") + 
    scale_colour_manual(values = c("navy","gray","darkred")) +
    theme_Publication() + 
    labs(color = "DEG")  +  ylab("-log10Pvalue") + 
    geom_vline(xintercept = c(-1*(fc.cutoff), 1*(fc.cutoff)), linetype="dashed") +
    geom_hline(yintercept = -log10(pval.cutoff), linetype="dashed") +
    ggtitle(title)
  
  return(p)
}

pdf(file = "results/limma/summary/volcano-plots.pdf", width = 8, height = 6)
plotVolcano(x = s.vs.tissue, title = "Suspension vs Solid Tissue\nTMM Genes", TM = TRUE, goi = common)
plotVolcano(x = a.vs.tissue, title = "Adherent vs Solid Tissue\nTMM Genes", TM = TRUE, goi = common)
dev.off()

# version 2 (combined plots and define shape by common genes)
plotVolcano <- function(x, TM = TRUE, title = "Volcano Plot", otherCol = c("blue","gray","red"), pval.cutoff = 0.05, fc.cutoff =  1, goi) {
  
  # filter by pvalue and logFC
  x <- x %>%
    mutate(sig = ifelse(adj.P.Val < pval.cutoff, T, F),
           DEGAnnot = ifelse(adj.P.Val < pval.cutoff & logFC > 1*fc.cutoff, 
                             yes = "Up", 
                             no = ifelse(test = adj.P.Val < pval.cutoff & logFC < -1*fc.cutoff, 
                                         yes = "Down", 
                                         no = "Unchanged")),
           tmm = ifelse(Gene %in% TMlist$hgnc_symbol, T, F))
  
  # filter to TMM
  if(TM == TRUE){
    x <- x %>%
      filter(tmm == TRUE)
  }
  
  # add info about common genes
  x <- x %>%
    group_by(Gene, DEGAnnot, sig) %>%
    mutate(n = n(),
           Genes = ifelse(n == 1, "Unique", "Common"))
  
  # count number of each type and add as labels
  x <- x %>%
    group_by(label, DEGAnnot) %>%
    mutate(freq = n()) %>%
    ungroup() %>%
    mutate(DEGAnnot.lab = paste0(DEGAnnot," (n = ", freq, ")"))
  
  # plot volcano
  p <- ggplot(x, aes(x = logFC, y = -log10(adj.P.Val), color = DEGAnnot)) +
    geom_point(alpha = 0.5, aes(shape = Genes), size = 3) +  
    facet_wrap(~label) +
    ggtitle("") + 
    scale_colour_manual(values = c("navy","gray","darkred")) +
    theme_Publication(base_size = 10) + theme(legend.box = "vertical") +
    labs(color = "DEG")  +  ylab("-log10Pvalue") + 
    geom_vline(xintercept = c(-1*(fc.cutoff), 1*(fc.cutoff)), linetype="dashed") +
    geom_hline(yintercept = -log10(pval.cutoff), linetype="dashed") +
    ggtitle(title)
  
  return(p)
}

pdf(file = "results/limma/summary/volcano-plots-combined.pdf", width = 12, height = 6)
s.vs.tissue$label <- 'S vs T'
a.vs.tissue$label <- 'A vs T'
x <- rbind(s.vs.tissue, a.vs.tissue)
plotVolcano(x = x, title = "Suspension and Adhesion vs Solid Tissue\nTMM Genes", TM = TRUE)
dev.off()

