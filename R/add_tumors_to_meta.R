# Author: Komal S. Rathi
# Date: 12/23/2019
# Function: Add metadata to excel file

library(xlsx)
library(tidyverse)

# Add patient tumors info to meta file
dat <- read.delim('~/Projects/OpenPBTA-analysis/data/pbta-histologies.tsv', stringsAsFactors = F)
meta <- read.xlsx('data/201910_Cell_line_byDerivedType_PBTA_KF_ID.xlsx', 1)
dat <- dat[which(dat$Kids_First_Participant_ID %in% meta$Kids_First_Participant_ID),]
dat <- dat %>% 
  filter(sample_type == "Tumor" & composition == "Solid Tissue")
dat$composition.type <- "solid_tissue"
meta <- meta[,intersect(colnames(dat), colnames(meta))]
dat <- dat[,intersect(colnames(dat), colnames(meta))]
dat <- rbind(meta, dat)
dat <- unique(dat)
dat[is.na(dat)] <- ''
write.xlsx(dat, file = 'data/201910_Cell_line_byDerivedType_PBTA_KF_ID.xlsx', sheetName = "v12", append = TRUE, row.names = F)

