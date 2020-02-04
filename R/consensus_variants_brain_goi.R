# Author: Komal S. Rathi
# Date: 02/03/2020
# Function: Simple Variants Summary

library(data.table)
library(dplyr)

# meta data
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
maf <- consensus.maf %>%
  filter(Hugo_Symbol %in% brain.goi$V1,
         Tumor_Sample_Barcode %in% meta$Kids_First_Biospecimen_ID) %>% # subset to brain goi
  select(Hugo_Symbol, HGVSp, HGVSp_Short, VAF, Variant_Classification, Tumor_Sample_Barcode) %>%
  inner_join(meta, 
             by = c("Tumor_Sample_Barcode" = "Kids_First_Biospecimen_ID"))
write.table(maf, file = 'results/variant_summary/consensus_variant_brain_goi.txt', quote = F,  sep = "\t", row.names = F)
