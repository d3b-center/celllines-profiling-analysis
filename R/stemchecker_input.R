# Author: Komal S. Rathi
# Date: 10/17/2019
# Function: Differential expression 
# unpaired suspension vs adhesion
# unpaired suspension vs solid tissue
# Code to create input for stemchecker and cibersort
# References: 
# Paper: https://www.cell.com/cell/pdf/S0092-8674(18)30358-1.pdf
# Tutorial: http://tcgabiolinks.fmrp.usp.br/PanCanStem/mRNAsi.html

library(xlsx)
library(dplyr)
library(tidyverse)
library(limma)

# create input expected count matrix
if(file.exists('data/input-counts.RData')){
  print("Counts file exists")
  load('data/input-counts.RData')
} else {
  # read meta file
  meta <- read.xlsx('data/201910_Cell_line_byDerivedType_PBTA_KF_ID.xlsx', 3)
  meta <- meta %>%
    filter(experimental_strategy == 'RNA-Seq') %>%
    dplyr::select(c(Kids_First_Biospecimen_ID,composition.type,tumor_descriptor,primary_site,reported_gender,sample_id))
  meta$label <- paste0(meta$sample_id,'_', tolower(meta$composition.type))
  rownames(meta) <- meta$Kids_First_Biospecimen_ID
  
  # read expected counts
  dat <- readRDS('data/pbta-gene-counts-rsem-expected_count.stranded.rds')
  common.cols <- intersect(colnames(dat), meta$Kids_First_Biospecimen_ID)
  meta <- meta <- meta %>%
    filter(Kids_First_Biospecimen_ID %in% common.cols)
  rownames(meta) <- meta$Kids_First_Biospecimen_ID
  # meta <- meta[which(meta$Kids_First_Biospecimen_ID %in% common.cols),] # "BS_BYCX6VK1" "BS_T3DGY9J9" missing tumors
  dat <- dat[,colnames(dat) %in% c('gene_id', common.cols)] # subset to meta file
  meta <- meta[colnames(dat)[2:ncol(dat)],]
  if(identical(colnames(dat)[2:ncol(dat)], as.character(meta$Kids_First_Biospecimen_ID))){
    colnames(dat)[2:ncol(dat)] <- meta$label
  } else {
    print("Please check mapping")
  }
  
  # reduce data
  dat <- dat[which(rowSums(dat[,2:ncol(dat)])>0),] # remove all rows with only 0
  dat <- dat[grep('_PAR_', dat$gene_id, invert = T),] # remove PAR_Y chr genes
  
  # convert into a matrix of gene symbol x sample name
  counts.collapsed <- dat %>% 
    separate(gene_id, c("gene_id", "gene_symbol"), sep = "\\_", extra = "merge") %>%
    dplyr::mutate(means = rowMeans(dplyr::select(.,-gene_id,-gene_symbol))) %>%
    arrange(desc(means)) %>%
    distinct(gene_symbol, .keep_all = TRUE) %>% 
    ungroup() %>% 
    dplyr::select(-means,-gene_id) %>%
    column_to_rownames(var="gene_symbol") %>%
    as.data.frame()
  
  save(counts.collapsed, file = 'data/input-counts.RData')
}

# perform differential expression
# counts.collapsed <- counts.collapsed[,grep('solid_tissue', colnames(counts.collapsed), invert = T)]
meta <- data.frame(id = colnames(counts.collapsed), type = gsub(".*_","",colnames(counts.collapsed)))
rownames(meta) <- meta$id

# voom transformation
mydesign <- model.matrix(~0+type, meta)
colnames(mydesign) <- sub('type','',colnames(mydesign))
tmp.voom <- voom(counts = as.matrix(counts.collapsed), design = mydesign)
tmp.voom <- tmp.voom$E

# 1. unpaired analysis
# differential expression
fit <- lmFit(tmp.voom, design = mydesign)
myConts <- c(paste0('s','-','a'), paste0('s','-','tissue'))
print(paste0("Contrast: ", myConts))
contrast.matrix = makeContrasts(contrasts = myConts, levels = mydesign)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

# Suspension vs Adhesion
tmpOut <- topTable(fit2, coef = 1, number = Inf, p.value = 0.05)[,c("logFC", "P.Value", "adj.P.Val")]
tmpOut$Gene <- rownames(tmpOut)
tmpOut <- tmpOut[,c("Gene","logFC","adj.P.Val")]
up <- tmpOut[which(tmpOut$logFC > 0),] # 2966
down <- tmpOut[which(tmpOut$logFC < 0),] # 832
# for Metacore/IPA
write.table(tmpOut[which(abs(tmpOut$logFC) > 1),], file = "results/limma/s_vs_a_limma_logfc1.txt", quote = F, row.names = F, sep = "\t")
write.table(tmpOut[which(tmpOut$logFC > 1),], file = "results/limma/s_vs_a_up_logfc1.txt", quote = F, row.names = F, sep = "\t")
write.table(tmpOut[which(tmpOut$logFC < -1),], file = "results/limma/s_vs_a_down_logfc1.txt", quote = F, row.names = F, sep = "\t")
# for stemchecker
write.table(up$Gene, file = 'data/stemchecker-input-svsa-up.txt', quote = F, row.names = F, col.names = F)
write.table(down$Gene, file = 'data/stemchecker-input-svsa-down.txt', quote = F, row.names = F, col.names = F)

# Suspension vs Solid Tissue
tmpOut <- topTable(fit2, coef = 2, number = Inf, p.value = 0.05)[,c("logFC", "P.Value", "adj.P.Val")]
tmpOut$Gene <- rownames(tmpOut)
tmpOut <- tmpOut[,c("Gene","logFC","adj.P.Val")]
up <- tmpOut[which(tmpOut$logFC > 0),] # 655
down <- tmpOut[which(tmpOut$logFC < 0),] # 2271
# for Metacore/IPA
write.table(tmpOut[which(abs(tmpOut$logFC) > 1),], file = "results/limma/s_vs_tissue_limma_logfc1.txt", quote = F, row.names = F, sep = "\t")
write.table(tmpOut[which(tmpOut$logFC > 1),], file = "results/limma/s_vs_tissue_up_logfc1.txt", quote = F, row.names = F, sep = "\t")
write.table(tmpOut[which(tmpOut$logFC < -1),], file = "results/limma/s_vs_tissue_down_logfc1.txt", quote = F, row.names = F, sep = "\t")
# for stemchecker
write.table(rownames(up), file = 'data/stemchecker-input-svstissue-up.txt', quote = F, row.names = F, col.names = F)
write.table(rownames(down), file = 'data/stemchecker-input-svstissue-down.txt', quote = F, row.names = F, col.names = F)

# # 2. paired analysis - not feasible
# meta$pairs <- gsub('_.*','',meta$id)
# meta <- meta[-which(meta$pairs %in% c('7316-1769','7316-2176','7316-2151','7316-2189')),]
# counts.paired <- counts.collapsed[,rownames(meta)]
# if(identical(colnames(counts.paired), rownames(meta))){
#   print("Matched samples")
# } else {
#   return(NULL)
# }
# 
# # voom transformation
# mydesign <- model.matrix(~0+type, meta)
# colnames(mydesign) <- sub('type','',colnames(mydesign))
# tmp.voom <- voom(counts = as.matrix(counts.paired), design = mydesign)
# tmp.voom <- tmp.voom$E
# 
# # differential expression
# pairs <- factor(meta$pairs)
# type <- factor(meta$type, levels = c("a","s"))
# mydesign <- model.matrix(~pairs+type)
# fit <- lmFit(tmp.voom, mydesign)
# fit2 <- eBayes(fit)
# 
# # format output and split into up and down genes
# tmpOut <- topTable(fit2, coef = "types", number = Inf, p.value = 0.05) # 0 rows
