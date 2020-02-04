# Author: Komal S. Rathi
# Date: 02/03/2019
# Function: Script for boxplot
# Stratified by group

# Run
# Rscript create_boxplot.R -e fpkm-matrix.RData -g MYC -l TRUE -o MYC_expr.pdf

# install packages if not present
list.of.packages <- c("ggplot2", "reshape2", "ggpubr")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) BiocManager::install(new.packages)

# load packages
library(ggplot2)
library(reshape2)
library(ggpubr)
library(optparse)

# arguments
option_list <- list(
  make_option(c("-e", "--expr"), type = "character",
              help = "Expression data: HUGO symbol x Sample identifiers (.RData)"),
  make_option(c("-g", "--gene"), type = "character", 
              help = "Gene Name"),
  make_option(c("-l", "--log"), type = "character",
              help = "Log (TRUE or FALSE)"),
  make_option(c("-o","--outputfile"), type = "character",
              help = "Output PDF")
)

# load data
opt <- parse_args(OptionParser(option_list = option_list))
dat <- opt$expr
dat <- get(load(dat))
gene <- opt$gene
log <- opt$log
output <- opt$outputfile

create.box <- function(dat, gene, log, fname) {
  my_comparisons <- list(c("s", "a"), c("s", "solid_tissue"), c("a", "solid_tissue"))
  dat.sub <- melt(dat[rownames(dat) == gene,])
  dat.sub <-  cbind(dat.sub, colsplit(dat.sub$variable, pattern = "_", names = c("sample_name","type")))
  
  if(log) {
    dat.sub$value <- log2(dat.sub$value + 1)
    y.lab <- "log2FPKM"
  } else {
    y.lab <- "FPKM"
  }
  pdf(file = fname, width = 8, height = 8)
  p <- ggplot(dat.sub, aes(x = type, y = value, fill = type)) +
    stat_boxplot(geom ='errorbar', width = 0.2) +
    geom_boxplot(outlier.shape = 21, outlier.fill = "white", outlier.color = "white",
                 lwd = 0.5, fatten = 0.7, width = 0.5) + 
    ggtitle(paste0('Gene Expression: ', gene)) +
    theme_bw() + ylab(y.lab) + xlab("") + 
    guides(fill = F) + geom_jitter(position=position_jitter(width=.25), shape = 21) +
    scale_x_discrete(labels = c("s" = "Suspension", 
                                "a" = "Adherent",
                                "solid_tissue" = "Patient Tumor")) +
    scale_fill_manual(values = c("s" = "#F8766D",
                                 "a" = "#00BFC4",
                                 "solid_tissue" = "#00BA38")) +
    stat_compare_means(comparisons = my_comparisons, color = "darkred", size = 3) +
    stat_compare_means(color = "darkred") +
    labs(fill = "Type")
  print(p)
  dev.off()
}

print("Create boxplot...")
create.box(dat = dat, gene = gene, log = log, fname = output)
