# Author: Komal S. Rathi
# Date: 10/18/2019
# Function: Plots for CIBERSORT results

library(pheatmap)

# sort samples cibersort input
dat <- get(load('data/fpkm-matrix.RData'))
# dat <- read.delim('data/input-counts-cibersort.txt',  check.names = F)
annot <- data.frame(type = sub('.*_','', colnames(dat)), id = colnames(dat))
rownames(annot) <- colnames(dat)
annot <- annot[which(annot$type != "title"),]
annot$type <- ifelse(annot$type == "s", "Suspension", 
                     ifelse(annot$type == "a", "Adherence", "Solid Tissue"))
annot <- as.data.frame(annot[order(annot$type),])
dat <- cbind("!Sample_title" = rownames(dat), dat)
write.table(dat, file = 'data/input-cibersort.txt',  quote = F, sep = "\t", row.names = F)

# cibersort output matrix
cibersort <- read.delim('results/cibersort/CIBERSORT.Output_Job8.txt')
rownames(cibersort) <- cibersort$Input.Sample
mat <- cibersort[,2:23]

# annotation
annot <- data.frame(type = sub('.*_','', rownames(mat)), id = rownames(mat))
rownames(annot) <- rownames(mat)
annot$type <- ifelse(annot$type == "s", "Suspension", 
                     ifelse(annot$type == "a", "Adherence", "Solid Tissue"))
annot <- as.data.frame(annot[order(annot$type),])
mat <- mat[rownames(annot),]

# heatmap
pdf(file = 'results/cibersort/heatmap.pdf', width = 10, height = 8)
pheatmap(mat = t(mat), annotation_col = annot[,1, drop=F])
dev.off()

