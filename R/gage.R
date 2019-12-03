# Author: Komal S. Rathi
# Date: 11/26/2019
# Function: GAGE using MSigDB 

# Use GAGE to perform enrichment using the following signatures
# MSigDB TF Targets
# MSigDB BP
# Oligodendrocyte Progenitor Cells

# load libraries
library(gage)
library(msigdbr)

# load data
load('data/input.RData')
df <- data.frame(snames = colnames(expr.collapsed), colsplit(colnames(expr.collapsed), pattern = "_", names = c("sample","type")))
rownames(df) <- df$snames
ref.idx = grep('^s$', df$type, invert = T)
samp.idx = grep('^s$', df$type)

# TFT
tft <- msigdbr(category = "C3", subcategory = "TFT")
tft <- split(tft$human_gene_symbol, tft$gs_name)
tft.gsea <- gage(exprs = expr.collapsed, 
                 gsets = tft, 
                 ref = ref.idx, samp = samp.idx, 
                 compare = "unpaired")
tft.gsea <- as.data.frame(tft.gsea$greater[,1:5])
tft.gsea <- tft.gsea[which(tft.gsea$p.val < 0.05),]

# BP
bp <- msigdbr(category = "C5", subcategory = "BP")
bp <- split(bp$human_gene_symbol, bp$gs_name)
bp.gsea <- gage(exprs = expr.collapsed, 
                 gsets = bp, 
                 ref = ref.idx, samp = samp.idx, 
                 compare = "unpaired")
bp.gsea <- as.data.frame(bp.gsea$greater[,1:5])
bp.gsea <- bp.gsea[which(bp.gsea$p.val < 0.05),]

# OPC: oligodendrocyte progenitor cells signature
# names are not visible in Figure 4
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4152602/figure/F4/
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4152602/table/T2/?report=objectonly
opc.genes <- read.delim('data/OPC_genes.txt', header = F, stringsAsFactors = F)
opc.genes <- list(OPC = toupper(opc.genes$V1))
opc.genes <- gage(exprs = expr.collapsed, 
                gsets = opc.genes, 
                ref = ref.idx, samp = samp.idx, 
                compare = "unpaired")
opc.genes <- as.data.frame(opc.genes$greater[,1:5])

