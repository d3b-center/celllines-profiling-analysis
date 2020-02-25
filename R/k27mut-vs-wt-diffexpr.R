# Author: Komal S. Rathi
# Date: 02/24/2020
# Differential expression
# k27 mutated vs wild type comparison

library(limma)
library(xlsx)

# add k27 mutation information to meta file
k27 <- read.delim('data/K27_mutated_samples.txt')

meta <- read.xlsx('data/201910_Cell_line_byDerivedType_PBTA_KF_ID.xlsx', 3)
meta <- meta %>%
  filter(experimental_strategy == "WGS") %>%
  dplyr::select(c(Kids_First_Biospecimen_ID, sample_id, composition.type)) %>%
  mutate(sample_id = paste0(sample_id,'_',tolower(composition.type))) %>%
  unique
meta$composition.type <- ifelse(meta$composition.type == "S", "Suspension", ifelse(meta$composition.type == "A", "Adherent", "Solid_Tissue"))
meta <- meta %>%
  mutate(k27.mut = ifelse(Kids_First_Biospecimen_ID %in% k27$Tumor_Sample_Barcode,"Y","N"))

# perform differential expression on k27 mutated vs wild type samples
load('data/count-matrix.RData')

# group <- 'Suspension'
# countData <- counts.collapsed
# mData <- meta
diffexpr <- function(countData, mData, group) {
    mData <- mData %>%
      filter(composition.type == group) %>%
      column_to_rownames("sample_id")
    countData <- countData[,rownames(mData)]
    
    if(identical(rownames(mData), colnames(countData))){
      print(dim(countData))
    } else {
      break
    }
    
    # voom transformation
    mydesign <- model.matrix(~0+k27.mut, mData)
    colnames(mydesign) <- gsub('k27.mut','',colnames(mydesign))
    tmp.voom <- voom(counts = as.matrix(countData), design = mydesign)
    tmp.voom <- tmp.voom$E
    print(dim(tmp.voom))
    
    # differential expression
    fit <- lmFit(tmp.voom, design = mydesign)
    myConts <- c(paste0('Y','-','N'))
    print(paste0("Contrast: ", myConts))
    contrast.matrix = makeContrasts(contrasts = myConts, levels = mydesign)
    fit2 <- contrasts.fit(fit, contrast.matrix)
    fit2 <- eBayes(fit2)
    
    # topTable
    tmpOut <- topTable(fit2, coef = 1, number = Inf)[,c("logFC", "adj.P.Val")]
    tmpOut <- tmpOut %>% 
      rownames_to_column("Gene")
    tmpOut <- tmpOut %>%
      filter(adj.P.Val < 0.05)
    if(nrow(tmpOut) > 1) {
      return(tmpOut)
    } else {
      return(NULL)
    }
}

diffexpr(countData = counts.collapsed, mData = meta, group = "Suspension")
diffexpr(countData = counts.collapsed, mData = meta, group = "Adherent")
diffexpr(countData = counts.collapsed, mData = meta, group = "Solid_Tissue")
