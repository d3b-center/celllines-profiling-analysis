# Author: Komal S. Rathi
# Date: 12/02/2019
# Function: Simple and Structural Variant Summary

library(data.table)
library(tidyverse)

meta <- read.xlsx('data/201910_Cell_line_byDerivedType_PBTA_KF_ID.xlsx', 3)
meta <- meta %>%
  filter(experimental_strategy == "WGS") %>%
  dplyr::select(c(Kids_First_Biospecimen_ID,sample_id,composition.type))  %>%
  unique

# summarize structural variants
manta <- fread('data/pbta-sv-manta.tsv.gz')
manta <- manta  %>%
  filter(Kids.First.Biospecimen.ID.Tumor %in% meta$Kids_First_Biospecimen_ID)
manta <- merge(manta, meta, by.x = 'Kids.First.Biospecimen.ID.Tumor', by.y = 'Kids_First_Biospecimen_ID')
manta.summary <- manta %>% 
  group_by(sample_id, composition.type, SV.type) %>%
  dplyr::summarise(sum = n())  %>%
  pivot_wider(id_cols = c(sample_id, composition.type), names_from = SV.type, values_from = sum)
write.table(manta.summary, file = 'results/variant_summary/manta_count_per_class.txt', quote = F, sep = "\t", row.names = F)

# summarize simple variants
snv <- fread('data/pbta-snv-consensus-mutation.maf.tsv.gz')
snv <- snv  %>%
  filter(Tumor_Sample_Barcode %in% meta$Kids_First_Biospecimen_ID)
snv <- merge(snv, meta, by.x = 'Tumor_Sample_Barcode', by.y = 'Kids_First_Biospecimen_ID')
snv.summary <- snv %>% 
  group_by(sample_id, composition.type, Variant_Classification) %>%
  dplyr::summarise(sum = n())  %>%
  pivot_wider(id_cols = c(sample_id, composition.type), names_from = Variant_Classification, values_from = sum)
write.table(snv.summary, file = 'results/variant_summary/consensus_snv_count_per_class.txt', quote = F, sep = "\t", row.names = F)
