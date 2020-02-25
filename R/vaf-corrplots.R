# Author: Komal S. Rathi
# Date: 02/17/2020
# Function: VAF Correlation plots between Patient and Cell-lines

setwd('~/Projects/celllines-profiling-analysis/')

library(data.table)
library(tidyverse)
library(reshape2)
library(xlsx)
library(ggrepel)
source('R/pubTheme.R')

meta <- read.xlsx('data/201910_Cell_line_byDerivedType_PBTA_KF_ID.xlsx', 3)
meta <- meta %>%
  filter(experimental_strategy == "WGS") %>%
  dplyr::select(c(Kids_First_Biospecimen_ID, sample_id, composition.type))  %>%
  unique
meta$composition.type <- ifelse(meta$composition.type == "S", "Suspension", ifelse(meta$composition.type == "A", "Adhesion", "Solid_Tissue"))

# consensus maf
consensus.maf <- fread('data/pbta-snv-consensus-mutation.maf.tsv.gz')

# brain gene list
brain.goi <- read.delim('data/brain-goi-list-new.txt', stringsAsFactors = F, header = F)

# subset to brain goi and output select columns
consensus.maf <- consensus.maf %>%
  filter(Tumor_Sample_Barcode %in% meta$Kids_First_Biospecimen_ID) %>% # subset to samples
  select(Hugo_Symbol, HGVSp, HGVSp_Short, VAF, Variant_Classification, Tumor_Sample_Barcode) %>%
  inner_join(meta, by = c("Tumor_Sample_Barcode" = "Kids_First_Biospecimen_ID"))

samples <- unique(as.character(consensus.maf$sample_id))
create.scatter <- function(maf, type, sid) {
  maf <- maf %>%
    filter(composition.type %in% c("Solid_Tissue", type),
           HGVSp_Short != ".",
           Hugo_Symbol != "Unknown") %>%
    mutate(label = paste0(Hugo_Symbol, "_", HGVSp_Short),
           sym = ifelse(Hugo_Symbol %in% brain.goi$V1, T, F)) %>%
    pivot_wider(id_cols = c(label, sym), names_from = composition.type, values_from = VAF) %>%
    mutate_if(is.numeric, funs(replace_na(., 0))) %>%
    as.data.frame()
  
  if(nrow(maf) == 0) {
    return(NULL)
  }
  # patient or cell-line specific variants
  maf$group <- ifelse(maf[,type] != maf[,'Solid_Tissue'], 
                      ifelse(maf[,type] == 0, 'Solid_Tissue', 
                             ifelse(maf[,'Solid_Tissue'] == 0, type, 'Common')),'Common')
  
  # convert labels
  maf$label <- ifelse(maf$sym == TRUE, maf$label, '')
  
  # plot
  maf$group <- factor(maf$group, levels = c("Adhesion","Suspension","Common","Solid_Tissue"))
  p <- ggplot(maf, aes_string(x = "Solid_Tissue", y = type, color = "group")) +
    geom_point(size = 10,fill = 4, alpha = 1/6) + 
    scale_colour_manual(values = c("Common" = "gray34", 
                                   "Solid_Tissue" = "dodgerblue3", 
                                   "Suspension" = "firebrick3",
                                   "Adhesion" = "firebrick3")) +
    ggtitle(sid) +
    geom_vline(xintercept = 0.1, linetype = "dashed") + # Adding vertical intercept
    geom_hline(yintercept = 0.1, linetype = "dashed") + # Adding horizontal intercept line
    geom_text_repel(aes(label = label), size = 3.5, hjust = 0,vjust = 0, nudge_x = 0.005, point.padding = NA, segment.color = NA, show.legend = FALSE, xlim = c(0.02,NA),ylim = c(0.025,0.96)) + # if label column is one, then label the point with gene name and it's corresponding Protein Change value
    theme_Publication() + 
    theme(plot.title = element_text(hjust = 0.5))+ 
    xlim(0,1) + 
    ylim(0,1)
  return(p)
}


system('mkdir -p results/vaf-corrplots/')
for(i in 1:length(samples)){
  print(i)
  maf.sub <- consensus.maf %>%
    filter(sample_id == samples[i])
  types <- unique(maf.sub$composition.type)
  if(length(grep("Solid_Tissue", types)) == 0){
    next
  }
  types <- grep("Solid_Tissue", types, invert = TRUE, value = TRUE)
  if(length(types) == 0){
    next
  }
  
  # create scatter with each cell line type
  for(j in 1:length(types)){
    fname <- paste0('results/vaf-corrplots/',samples[i],'-',types[j],'-vaf.pdf')
    print(fname)
    p <- create.scatter(maf = maf.sub, type = types[j], sid = samples[i])
    if(is.null(p)) {
      print("No plot available")
      next
    }
    pdf(file = fname, width = 10, height = 8)
    print(p)
    dev.off()
  }
}
