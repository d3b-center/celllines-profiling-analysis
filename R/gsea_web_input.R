# Author: Komal S. Rathi
# Date: 11/26/2019
# Function: create GSEA input files (web version)
# The input files should be first normalized using the DESeq2 module of GenePattern

# libraries
library(tidyverse)
library(reshape2)

# expected counts
load('data/count-matrix.RData')
# counts.collapsed <- log2(counts.collapsed + 1)  

# counts.collapsed
# groups <- c("s", "a")
# gct.file <- 'data/gsea_input/input_mat_s_vs_a.gct'
# cls.file <- 'data/gsea_input/phenotype_s_vs_a.cls'
# control <- "#Adhesion"
# treat <- "Suspension"

gsea.input <- function(counts.collapsed, groups, gct.file, cls.file, control, treat) {
  # make file GSEA web version
  # expression file  -  counts file 
  df <- data.frame(snames = colnames(counts.collapsed), colsplit(colnames(counts.collapsed), pattern = "_", names = c("sample","type")))
  df <- df[order(df$type),]
  df <- df %>%  
    filter(type %in% groups)
  rownames(df) <- df$snames
  gct <- counts.collapsed[,rownames(df)]
  add <- data.frame(NAME = c("#1.2",nrow(gct),"NAME"), Description = c('', ncol(gct), "Description"))
  total.cols <- ncol(gct) + 2
  add[,3:total.cols] <- ''
  colnames(add)[3:total.cols] <- colnames(gct)
  add[3,3:total.cols] <- colnames(gct)
  annot <- data.frame(NAME = rownames(gct), Description = 'na')
  annot <- merge(annot, gct, by.x = 'NAME', by.y = 'row.names')
  add <- rbind(add, annot)
  write.table(add, file = gct.file, quote = F, sep = "\t", col.names = F, row.names = F)
  
  # phenotype file
  ph <- data.frame(V1 = c(ncol(gct), control,''), V2 = c('2', treat, ''), V3 = c('1','',''), stringsAsFactors = F)
  ph[,4:ncol(gct)] <- ''
  ph[3,] <- df$type
  write.table(ph, file = cls.file, quote = F, sep = " ", col.names = F, row.names = F)
  
  # chp file
  # chp <- data.frame('Probe Set ID' = rownames(gct), 'Gene Symbol' = rownames(gct), 'Gene Title' = 'na', check.names = F)
  # write.table(chp, file = 'data/gsea_input/mapping.chip', quote = F, sep = "\t", row.names = F)
}

gsea.input(counts.collapsed, groups = c("s","a"), 
           gct.file = 'data/gsea_input/raw/input_mat_s_vs_a.gct',
           cls.file = 'data/gsea_input/phenotype_s_vs_a.cls',
           control = '# Adhesion',
           treat = 'Suspension')

gsea.input(counts.collapsed, groups = c("s","solid_tissue"), 
           gct.file = 'data/gsea_input/raw/input_mat_s_vs_st.gct',
           cls.file = 'data/gsea_input/phenotype_s_vs_st.cls',
           control = '# Suspension',
           treat = 'Solid_Tissue')

gsea.input(counts.collapsed, groups = c("solid_tissue","a"), 
           gct.file = 'data/gsea_input/raw/input_mat_st_vs_a.gct',
           cls.file = 'data/gsea_input/phenotype_st_vs_a.cls',
           control = '# Adhesion',
           treat = 'Solid_Tissue')
