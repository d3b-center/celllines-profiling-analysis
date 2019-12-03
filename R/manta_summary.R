# Author: Komal S. Rathi
# Date: 12/02/2019
# Function: Summary of structural variants
# Median count of Str Variant Class per Composition Type

library(StructuralVariantAnnotation)

# meta file for cell line RNA-seq
meta <- read.xlsx('data/201910_Cell_line_byDerivedType_PBTA_KF_ID.xlsx', 2)
meta <- meta %>% filter(meta$experimental_strategy == 'WGS') %>%
  dplyr::select(c(Kids_First_Participant_ID, Kids_First_Biospecimen_ID, 
                  composition.type, tumor_descriptor,
                  primary_site, reported_gender, sample_id))

# meta file for manta
manta.meta <- read.csv('data/1575397670816-manifest-manta.csv')
manta.meta <- merge(manta.meta[,c("name", "Kids.First.Participant.ID", "Kids.First.Biospecimen.ID", "sample_id")], 
                      meta[,c("Kids_First_Participant_ID","Kids_First_Biospecimen_ID","sample_id","composition.type")], 
                    by.x = c("Kids.First.Participant.ID","Kids.First.Biospecimen.ID","sample_id"), 
                    by.y = c("Kids_First_Participant_ID", "Kids_First_Biospecimen_ID", "sample_id"))

breakpointgr2bedpe <- function(gr){
  bedpe <- data.frame(chrom1 = GenomeInfoDb::seqnames(gr), 
                      start1 = start(gr) - 1, end1 = end(gr), chrom2 = GenomeInfoDb::seqnames(partner(gr)), 
                      start2 = start(partner(gr)) - 1, end2 = end(partner(gr)), 
                      name = names(gr), partner.name = names(partner(gr)), 
                      score = gr$QUAL, strand1 = strand(gr), strand2 = strand(partner(gr)), 
                      svtype = gr$svtype)
  bedpe <- bedpe[(as.character(bedpe$chrom1) < as.character(bedpe$chrom2)) | 
                   (bedpe$chrom1 == bedpe$chrom2 & bedpe$start1 < bedpe$start2), 
                 -c(8)]
  return(bedpe)
}

# subset str. variant files
if(file.exists('data/manta/strData_manta.RData')){
  load('data/manta/strData_manta.RData')
} else {
  lf <- list.files('data/manta', pattern = '*.vcf',  full.names = T)
  lf <- lf[gsub('data/manta/','',lf) %in% manta.meta$name]
  
  # now apply on each file
  vcf.file <- lf[1]
  read.vcf <- function(x){
    vcf <- VariantAnnotation::readVcf(x, genome = "hg38")
    gr <- breakpointRanges(vcf)
    gr.bedpe <- breakpointgr2bedpe(gr)
    gr.bedpe$filename <- gsub('data/manta/', '', x)
    return(gr.bedpe)
  }
  
  strFiles <- lapply(lf, read.vcf)
  strData <- data.table::rbindlist(strFiles)
  strData <- merge(strData, manta.meta, by.x = 'filename', by.y = 'name')
  save(strData, file = 'data/manta/strData.RData')
  system('rm data/manta/*.vcf.gz')
}

strData.med.count <- strData %>% 
  group_by(Kids.First.Biospecimen.ID, composition.type, svtype) %>%
  summarise(sum = n())  %>%
  group_by(composition.type, svtype) %>%
  summarise(median = median(sum)) %>%
  pivot_wider(id_cols = svtype, names_from = composition.type, values_from = median)
colnames(strData.med.count)[2:4] <- c("Adhesion", "Suspension", "Solid_Tissue")
write.table(strData.med.count, file = 'results/manta_median_count_per_class.txt', quote = F, sep = "\t", row.names = F)
