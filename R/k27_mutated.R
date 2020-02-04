# Author: Komal S. Rathi
# Date: 02/03/2020
# Function: K27 Mutated samples for
# Proteasomal Pathway enrichment

library(data.table)
library(xlsx)
library(dplyr)

# consensus maf
dat <- fread('data/pbta-snv-consensus-mutation.maf.tsv.gz')

# meta data
meta <- read.xlsx('data/201910_Cell_line_byDerivedType_PBTA_KF_ID.xlsx', 3)
meta <- meta %>%
  filter(experimental_strategy == "WGS") %>%
  dplyr::select(c(Kids_First_Biospecimen_ID, sample_id, composition.type))  %>%
  unique
meta$composition.type <- ifelse(meta$composition.type == "S", "s", ifelse(meta$composition.type == "A", "a", "solid_tissue"))

# get samples which have K27 mutation
dat.sub <- dat %>%
  filter(Tumor_Sample_Barcode %in% meta$Kids_First_Biospecimen_ID,
         HGVSp != ".",
         Hugo_Symbol %in% c("H3F3A", "HISH1H3B")) %>%
  inner_join(meta, by = c("Tumor_Sample_Barcode" = "Kids_First_Biospecimen_ID")) %>%
  select(Tumor_Sample_Barcode, sample_id, composition.type) %>%
  unique
write.table(dat.sub, file = "data/K27_mutated_samples.txt", quote = F, row.names = F, sep = "\t")
