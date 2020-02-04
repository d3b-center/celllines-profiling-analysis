# Author: Komal S. Rathi
# Date: 02/03/2019
# Function: Venn diagram of DEGs with Surface proteins

library(dplyr)
library(VennDiagram)
library(gridExtra)

# read DEGs
s.vs.a <- read.delim('results/limma/sus_vs_adh_limma.txt')
s.vs.tissue <- read.delim('results/limma/sus_vs_tissue_limma.txt')
a.vs.tissue <- read.delim('results/limma/adh_vs_tissue_limma.txt')

# surface proteins
load('~/Projects/marislab-webportal/data/TMlist.RData')
s.vs.a$TM <-  ifelse(s.vs.a$Gene %in% TMlist$hgnc_symbol, "Yes", "No")
s.vs.a <- s.vs.a %>% filter(TM == "Yes")
s.vs.tissue$TM <-  ifelse(s.vs.tissue$Gene %in% TMlist$hgnc_symbol, "Yes", "No")
s.vs.tissue <- s.vs.tissue %>% filter(TM == "Yes")
a.vs.tissue$TM <-  ifelse(a.vs.tissue$Gene %in% TMlist$hgnc_symbol, "Yes", "No")
a.vs.tissue <- a.vs.tissue %>% filter(TM == "Yes")

# venn diagram
one <- s.vs.a$Gene
two <- s.vs.tissue$Gene
three <- a.vs.tissue$Gene
pdf(file = 'results/DEG_TM_Venn.pdf', width = 6, height = 5, onefile = FALSE)
p <- draw.triple.venn(area1 = length(one),
                      area2 = length(two),
                      area3 = length(three),
                      n12 = length(intersect(one, two)),
                      n23 = length(intersect(two, three)),
                      n13 = length(intersect(one, three)),
                      n123 = length(intersect(intersect(one, two), three)),
                      scaled = F, euler.d = F,
                      category = c(paste0("S vs A\n(",length(one),")"),
                                   paste0("S vs ST\n(",length(two),")"),
                                   paste0("A vs ST\n(",length(three),")")),
                      fill = c("skyblue", "pink1", "mediumorchid"), lty = "blank")
grid.arrange(gTree(children=p), top="Differentially Expressed Surface Genes")
dev.off()
