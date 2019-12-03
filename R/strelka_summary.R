# Author: Komal S. Rathi
# Date: 12/02/2019
# Function: Summary of somatic variants
# Median count of Variant Classes per Composition Type

library(tidyverse)

# meta file for cell line RNA-seq
meta <- read.xlsx('data/201910_Cell_line_byDerivedType_PBTA_KF_ID.xlsx', 2)
meta <- meta %>% filter(meta$experimental_strategy == 'WGS') %>%
  dplyr::select(c(Kids_First_Participant_ID, Kids_First_Biospecimen_ID, 
           composition.type, tumor_descriptor,
           primary_site, reported_gender, sample_id))

# meta file for strelka
strelka.meta <- read.csv('data/1575324497393-manifest-strelka.csv')
strelka.meta <- merge(strelka.meta[,c("name", "Kids.First.Participant.ID", "Kids.First.Biospecimen.ID", "sample_id")], 
                      meta[,c("Kids_First_Participant_ID","Kids_First_Biospecimen_ID","sample_id","composition.type")], by.x = c("Kids.First.Participant.ID","Kids.First.Biospecimen.ID","sample_id"), by.y = c("Kids_First_Participant_ID", "Kids_First_Biospecimen_ID", "sample_id"))

# subset mutation files
if(file.exists('data/strelka/mutData_strelka.RData')){
  load('data/strelka/mutData_strelka.RData')
} else {
  lf <- list.files('data/strelka', pattern = '*.maf',  full.names = T)
  lf <- lf[gsub('data/strelka/','',lf) %in% strelka.meta$name]
  mutFiles <- lapply(lf, data.table::fread, skip = 1, stringsAsFactors = F)
  mutData <- data.table::rbindlist(mutFiles)
  mutData <- mutData %>% 
    as.data.frame() %>%
    filter(FILTER == "PASS")  %>%
    unique()
  mutData <- mutData[,colSums(is.na(mutData)) < nrow(mutData)]
  save(mutData, file = 'data/strelka/mutData.RData')
  system('rm data/strelka/*.maf')
}

# merge with meta data and summarise by composition type
mutData <- merge(mutData, strelka.meta, by = 'Tumor_Sample_Barcode', by.y = 'Kids.First.Biospecimen.ID')
mutData.med.count <- mutData %>% 
  group_by(Tumor_Sample_Barcode, composition.type, Variant_Classification) %>%
  summarise(sum = n())  %>%
  group_by(composition.type, Variant_Classification) %>%
  summarise(median = median(sum)) %>%
  pivot_wider(id_cols = Variant_Classification, names_from = composition.type, values_from = median)
colnames(mutData.med.count)[2:4] <- c("Adhesion", "Suspension", "Solid_Tissue")
write.table(mutData.med.count, file = 'results/strelka_median_count_per_class.txt', quote = F, sep = "\t", row.names = F)
