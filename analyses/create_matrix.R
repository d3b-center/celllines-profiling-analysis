# Author: Komal S. Rathi
# Date: 02/03/2020
# Function: Create Count and FPKM matrices

library(dplyr)
library(xlsx)

# create input expected count matrix
create.mat <- function(type = c("fpkm", "counts")){
  # read meta file
  meta <- read.xlsx('data/201910_Cell_line_byDerivedType_PBTA_KF_ID.xlsx', 3)
  meta <- meta %>%
    filter(experimental_strategy == 'RNA-Seq') %>%
    dplyr::select(c(Kids_First_Biospecimen_ID,composition.type,tumor_descriptor,primary_site,reported_gender,sample_id))
  meta$label <- paste0(meta$sample_id,'_', tolower(meta$composition.type))
  rownames(meta) <- meta$Kids_First_Biospecimen_ID
  
  # read expected counts
  if(type == "counts") {
    dat <- readRDS('data/pbta-gene-counts-rsem-expected_count.stranded.rds')
    fname <- 'data/count-matrix.RData'
  } else {
    dat <- readRDS('data/pbta-gene-expression-rsem-fpkm.stranded.rds')
    fname <- 'data/fpkm-matrix.RData'
  }
  
  common.cols <- intersect(colnames(dat), meta$Kids_First_Biospecimen_ID)
  meta <- meta <- meta %>%
    filter(Kids_First_Biospecimen_ID %in% common.cols)
  rownames(meta) <- meta$Kids_First_Biospecimen_ID
  write.table(meta, file = 'data/metadata_hgg.txt', quote = F, sep = '\t', row.names = F)
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
  
  save(counts.collapsed, file = fname)
}

create.mat(type = "fpkm")
create.mat(type = "counts")
