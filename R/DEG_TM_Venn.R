# Author: Komal S. Rathi
# Date: 02/17/2019
# Function: Pairwise venn diagram of DEGs with Surface proteins

library(dplyr)
library(VennDiagram)
library(gridExtra)

# read DEGs
s.vs.tissue <- read.delim('results/limma/s_vs_tissue_limma_sig.txt', stringsAsFactors = F)
a.vs.tissue <- read.delim('results/limma/a_vs_tissue_limma_sig.txt', stringsAsFactors = F)

# surface proteins
load('~/Projects/marislab-webportal/data/TMlist.RData')
s.vs.tissue$TM <-  ifelse(s.vs.tissue$Gene %in% TMlist$hgnc_symbol, "Yes", "No")
s.vs.tissue <- s.vs.tissue %>% filter(TM == "Yes")
s.vs.tissue$Label <- "Suspension_vs_SolidTissue"

a.vs.tissue$TM <-  ifelse(a.vs.tissue$Gene %in% TMlist$hgnc_symbol, "Yes", "No")
a.vs.tissue <- a.vs.tissue %>% filter(TM == "Yes")
a.vs.tissue$Label <- "SolidTissue_vs_Adhesion"

create.venn <- function(direction, fname){
  
  # table
  total <- rbind(s.vs.tissue, a.vs.tissue)
  
  if(direction == "down") {
    one <- s.vs.tissue %>% filter(logFC < 0) %>% .$Gene
    two <- a.vs.tissue %>% filter(logFC < 0) %>% .$Gene
    title <- "Downregulated Surface Genes"
    total <- total %>% filter(logFC < 0) %>% 
      group_by(Gene) %>% 
      summarise(Label = toString(Label), N = n()) %>%
      arrange(desc(N)) %>%
      as.data.frame()
  } else if(direction == "up") {
    one <- s.vs.tissue %>% filter(logFC > 0) %>% .$Gene
    two <- a.vs.tissue %>% filter(logFC > 0) %>% .$Gene
    title <- "Upregulated Surface Genes"
    total <- total %>% filter(logFC > 0) %>% 
      group_by(Gene) %>% 
      summarise(Label = toString(Label), N = n()) %>%
      arrange(desc(N)) %>%
      as.data.frame()
  } else {
    one <- s.vs.tissue %>% .$Gene
    two <- a.vs.tissue %>% .$Gene
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
  p <- draw.pairwise.venn(area1 = length(one),
                          area2 = length(two), 
                          cross.area = length(intersect(one, two)),
                          scaled = F, euler.d = F, cat.pos = 0,
                          category = c(paste0("S vs ST\n(",length(one),")"),
                                       paste0("A vs ST\n(",length(two),")")),
                          fill = c("skyblue", "pink1"), lty = "blank")
  # p <- draw.triple.venn(area1 = length(one),
  #                       area2 = length(two),
  #                       area3 = length(three),
  #                       n12 = length(intersect(one, two)),
  #                       n23 = length(intersect(two, three)),
  #                       n13 = length(intersect(one, three)),
  #                       n123 = length(intersect(intersect(one, two), three)),
  #                       scaled = F, euler.d = F,
  #                       category = c(paste0("S vs A\n(",length(one),")"),
  #                                    paste0("S vs ST\n(",length(two),")"),
  #                                    paste0("ST vs A\n(",length(three),")")),
  #                       fill = c("skyblue", "pink1", "mediumorchid"), lty = "blank")
  grid.arrange(gTree(children=p), top = title)
  dev.off()
  
}

create.venn(direction = "up", fname = "DEG_UP_TM")
create.venn(direction = "down", fname = "DEG_DOWN_TM")
create.venn(direction = "", fname = "DEG_TM")


