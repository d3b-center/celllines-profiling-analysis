# Author: Komal S. Rathi
# Date: 02/03/2019
# Function: Venn diagram of DEGs with Surface proteins

library(dplyr)
library(VennDiagram)
library(gridExtra)

# read DEGs
s.vs.a <- read.delim('results/limma/s_vs_a_limma.txt', stringsAsFactors = F)
s.vs.tissue <- read.delim('results/limma/s_vs_tissue_limma.txt', stringsAsFactors = F)
tissue.vs.a <- read.delim('results/limma/tissue_vs_a_limma.txt', stringsAsFactors = F)

# surface proteins
load('~/Projects/marislab-webportal/data/TMlist.RData')
s.vs.a$TM <-  ifelse(s.vs.a$Gene %in% TMlist$hgnc_symbol, "Yes", "No")
s.vs.a <- s.vs.a %>% filter(TM == "Yes")
s.vs.a$Label <- "Suspension_vs_Adhesion"

s.vs.tissue$TM <-  ifelse(s.vs.tissue$Gene %in% TMlist$hgnc_symbol, "Yes", "No")
s.vs.tissue <- s.vs.tissue %>% filter(TM == "Yes")
s.vs.tissue$Label <- "Suspension_vs_SolidTissue"

tissue.vs.a$TM <-  ifelse(tissue.vs.a$Gene %in% TMlist$hgnc_symbol, "Yes", "No")
tissue.vs.a <- tissue.vs.a %>% filter(TM == "Yes")
tissue.vs.a$Label <- "SolidTissue_vs_Adhesion"

create.venn <- function(direction, fname){
  
  # table
  total <- rbind(s.vs.a, s.vs.tissue, tissue.vs.a)
  
  if(direction == "down") {
    one <- s.vs.a %>% filter(logFC < 0) %>% .$Gene
    two <- s.vs.tissue %>% filter(logFC < 0) %>% .$Gene
    three <- tissue.vs.a %>% filter(logFC < 0) %>% .$Gene
    title <- "Downregulated Surface Genes"
    total <- total %>% filter(logFC < 0) %>% 
      group_by(Gene) %>% 
      summarise(Label = toString(Label), N = n()) %>%
      arrange(desc(N)) %>%
      as.data.frame()
  } else if(direction == "up") {
    one <- s.vs.a %>% filter(logFC > 0) %>% .$Gene
    two <- s.vs.tissue %>% filter(logFC > 0) %>% .$Gene
    three <- tissue.vs.a %>% filter(logFC > 0) %>% .$Gene
    title <- "Upregulated Surface Genes"
    total <- total %>% filter(logFC > 0) %>% 
      group_by(Gene) %>% 
      summarise(Label = toString(Label), N = n()) %>%
      arrange(desc(N)) %>%
      as.data.frame()
  } else {
    one <- s.vs.a %>% .$Gene
    two <- s.vs.tissue %>% .$Gene
    three <- tissue.vs.a %>% .$Gene
    title  <- "Differentially Expressed Surface Genes"
    total <- total %>% group_by(Gene) %>% 
      mutate(N = n()) %>%
      select(-c(adj.P.Val, TM)) %>%
      arrange(desc(N), Gene) %>%
      as.data.frame()
  }
  
  # table
  write.xlsx(x = total, file = file.path("results/limma/summary", paste0(fname, "_Table.xlsx")), append = FALSE, row.names = FALSE)
  
  # venn diagram
  pdf(file = file.path("results/limma/summary", paste0(fname, "_Venn.pdf")), width = 6, height = 5, onefile = FALSE)
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
                                     paste0("ST vs A\n(",length(three),")")),
                        fill = c("skyblue", "pink1", "mediumorchid"), lty = "blank")
  grid.arrange(gTree(children=p), top = title)
  dev.off()
  
}

create.venn(direction = "up", fname = "DEG_UP_TM")
create.venn(direction = "down", fname = "DEG_DOWN_TM")
create.venn(direction = "", fname = "DEG_TM")


