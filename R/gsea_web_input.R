# Author: Komal S. Rathi
# Date: 11/26/2019
# Function: create GSEA input files (web version)

load('data/input-counts.RData')
counts.collapsed <- log2(counts.collapsed +1)

# make file GSEA web version
# expression file
df <- data.frame(snames = colnames(counts.collapsed), colsplit(colnames(counts.collapsed), pattern = "_", names = c("sample","type")))
df <- df[order(df$type),]
rownames(df) <- df$snames
gct <- counts.collapsed[,rownames(df)]
add <- data.frame(NAME = c("#1.2",nrow(gct),"NAME"), Description = c('', ncol(gct), "Description"))
add[,3:28] <- ''
colnames(add)[3:28] <- colnames(gct)
add[3,3:28] <- colnames(gct)
annot <- data.frame(NAME = rownames(gct), Description = 'na')
annot <- merge(annot, gct, by.x = 'NAME', by.y = 'row.names')
add <- rbind(add, annot)
write.table(add, file = 'data/gsea_input/input_mat.gct', quote = F, sep = "\t", col.names = F, row.names = F)

# phenotype file
ph <- data.frame(V1 = c(ncol(gct),"#Adhesion",''), V2 = c('3', 'Suspension',''), V3 = c(1,'Solid_Tissue',''), stringsAsFactors = F)
ph[,4:26] <- ''
ph[3,] <- df$type
write.table(ph, file = 'data/gsea_input/phenotype.cls', quote = F, sep = "\t", col.names = F, row.names = F)

# chp file
# chp <- data.frame('Probe Set ID' = rownames(gct), 'Gene Symbol' = rownames(gct), 'Gene Title' = 'na', check.names = F)
# write.table(chp, file = 'data/gsea_input/mapping.chip', quote = F, sep = "\t", row.names = F)

