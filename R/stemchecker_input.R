# Author: Komal S. Rathi
# Date: 10/17/2019
# Function: Differential expression (unpaired suspension vs adhesion)
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
  meta <- read.xlsx('data/201910_Cell_line_byDerivedType_PBTA_KF_ID.xlsx', 2)
  meta <- meta[meta$experimental_strategy == 'RNA-Seq',]
  meta <- meta[,c('Kids_First_Biospecimen_ID','composition.type','tumor_descriptor','primary_site','reported_gender','sample_id')]
  plyr::count(meta$composition.type) # 9 A, 7 S, 12 solid_tissue
  meta$label <- paste0(meta$sample_id,'_', tolower(meta$composition.type))
  rownames(meta) <- meta$Kids_First_Biospecimen_ID
  
  # read expected counts
  dat <- readRDS('data/pbta-gene-counts-rsem-expected_count.stranded.rds')
  common.cols <- intersect(colnames(dat), meta$Kids_First_Biospecimen_ID)
  meta <- meta[which(meta$Kids_First_Biospecimen_ID %in% common.cols),] # "BS_BYCX6VK1" "BS_T3DGY9J9" missing tumors
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
    dplyr::mutate(means = rowMeans(select(.,-gene_id,-gene_symbol))) %>%
    arrange(desc(means)) %>%
    distinct(gene_symbol, .keep_all = TRUE) %>% 
    ungroup() %>% 
    dplyr::select(-means,-gene_id) %>%
    column_to_rownames(var="gene_symbol") %>%
    as.data.frame()
  
  save(counts.collapsed, file = 'data/input-counts.RData')
}

# perform differential expression
counts.collapsed <- counts.collapsed[,grep('solid_tissue', colnames(counts.collapsed), invert = T)]
meta <- data.frame(id = colnames(counts.collapsed), type = gsub(".*_","",colnames(counts.collapsed)))

# voom transformation
mydesign <- model.matrix(~0+type, meta)
rownames(mydesign) <- meta$id
colnames(mydesign) <- sub('type','',colnames(mydesign))
tmp.voom <- voom(counts = as.matrix(counts.collapsed), design = mydesign)
tmp.voom <- tmp.voom$E

# differential expression
fit <- lmFit(tmp.voom, design = mydesign)
myConts <- paste0('s','-','a')
print(paste0("Contrast: ", myConts))
contrast.matrix = makeContrasts(contrasts = myConts, levels = mydesign)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

# format output and split into up and down genes
tmpOut <- topTable(fit2, number = Inf, p.value = 0.05)[,c("logFC", "P.Value", "adj.P.Val")]
up <- tmpOut[which(tmpOut$logFC > 0),]
down <- tmpOut[which(tmpOut$logFC < 0),]
write.table(rownames(up), file = 'data/stemchecker-input-up.txt', quote = F, row.names = F, col.names = F)
write.table(rownames(down), file = 'data/stemchecker-input-down.txt', quote = F, row.names = F, col.names = F)

