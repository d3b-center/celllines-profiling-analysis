# Author: Komal S. Rathi
# Date: 10/11/2019
# Function: Run OCLR based Stemness profiling
# References: 
  # Paper: https://www.cell.com/cell/pdf/S0092-8674(18)30358-1.pdf
  # Tutorial: http://tcgabiolinks.fmrp.usp.br/PanCanStem/mRNAsi.html

setwd('~/Projects/celllines-profiling-analysis/')
library(xlsx)
library(plyr)
library(dplyr)
library(reshape2)
library(biomaRt)
library(gelnet)
library(tidyverse)
library(ggpubr)
library(synapser) # use synapser instead of synapseClient

# Synapse login required to get PCBC data
SYNAPSE_ID <- Sys.getenv("SYNAPSE_ID")
SYNAPSE_PWD <- Sys.getenv("SYNAPSE_PWD")
synLogin(email = SYNAPSE_ID, password = SYNAPSE_PWD) # your id and password from https://www.synapse.org/

# source R scripts
source('R/pubTheme.R') # themes
source('R/genes2hugo.R') # convert entrez ids to hugo symbols
source('R/main.train.R') # train using PCBC dataset
source('R/main.predict.R') # predict using vector generated in main.train.R

# run training set (the output of this should be vector of length 78 containing only 1s)
# stem/progenitor cells from the Progenitor Cell Biology Consortium
if(file.exists('data/pcbc-stemsig.tsv')){
  print("Training set ready")
} else {
  print("Run training set")
  auc <- main.train(fnOut = 'data/pcbc-stemsig.tsv', fnGenes = NULL)
}

# load FPKM matrix
load('data/fpkm-matrix.RData')

# use predict function
if(file.exists('data/mRNA_StemScore.tsv')){
  print("Predicted scores ready")
} else {
  main.predict(fnSig = 'data/pcbc-stemsig.tsv', expr = 'data/fpkm-matrix.RData', fnOut = "data/mRNA_StemScore.tsv")
}

# visualization
predicted.stemness.scores <- read.delim('data/mRNA_StemScore.tsv', header = F)
colnames(predicted.stemness.scores) <- c("sample_name", "stemness_score")
predicted.stemness.scores <- cbind(colsplit(predicted.stemness.scores$sample_name, pattern = "_", names = c("sample_name","type")), 
                                   stemness_score = predicted.stemness.scores$stemness_score)

# boxplot
shapiro.test(x = predicted.stemness.scores[predicted.stemness.scores$type == "a",3]) # p-value 0.001936
shapiro.test(x = predicted.stemness.scores[predicted.stemness.scores$type == "s",3]) # p-value 0.3465
shapiro.test(x = predicted.stemness.scores[predicted.stemness.scores$type == "solid_tissue",3]) # p-value 0.05523

# let's use Kruskal Wallis here as type a is not normally distributed
# Adherent outliers: 7316-3058, 7316-913
# Solid tissue outliers: 7316-913
predicted.stemness.scores$type <- factor(predicted.stemness.scores$type, levels = c("a", "solid_tissue", "s"))
my_comparisons <- list(c("s", "a"), c("s", "solid_tissue"), c("a", "solid_tissue"))
predicted.stemness.scores$label <- ifelse(predicted.stemness.scores$type %in% c("a","solid_tissue") & predicted.stemness.scores$stemness_score > 0.5, predicted.stemness.scores$sample_name, "")
p <- ggplot(predicted.stemness.scores, aes(x = type, y = stemness_score, fill = type)) +
  stat_boxplot(geom ='errorbar', width = 0.2) +
  geom_boxplot(outlier.shape = 21, outlier.fill = "white", outlier.color = "white",
               lwd = 0.5, fatten = 0.7, width = 0.5) + ggtitle('OCLR predicted Stemness index\nAdherent vs Suspension Cell types') +
  theme_bw() + theme_Publication(base_size = 10) + ylab('Stemness Index') + xlab("") + 
  geom_text(aes(label = label), size = 3, hjust = 0.7, vjust = -0.5) +
  guides(fill = F) + geom_jitter(position=position_jitter(width=.25), shape = 21) +
  scale_x_discrete(labels = c("s" = "Suspension", 
                              "a" = "Adherent",
                              "solid_tissue" = "Patient Tumor")) +
  scale_fill_manual(values = c("s" = "#F8766D",
                               "a" = "#00BFC4",
                               "solid_tissue" = "#00BA38")) +
  stat_compare_means(comparisons = my_comparisons, color = "darkred", size = 3) +
  stat_compare_means(label.y = 1.4, color = "darkred") +
  labs(fill = "Type")
p
ggsave(plot = p, filename = 'results/oclr_stemness_index.png', device = 'png', width = 5, height = 5)
