# Author: Komal S. Rathi
# Date: 10/11/2019
# Function: Run OCLR based Stemness profiling

setwd('~/Projects/celllines-profiling-analysis/')
library(xlsx)
library(plyr)
library(dplyr)
library(reshape2)
library(biomaRt)
library(gelnet)
library(tidyverse)
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

# prepare test set
if(file.exists('data/input.RData')){
  print("Test set ready")
} else {
  print("Prepate test set")
  # read meta file
  meta <- read.xlsx('data/201910_Cell_line_byDerivedType_PBTA_KF_ID.xlsx', 1)
  meta <- meta[meta$experimental_strategy == 'RNA-Seq',]
  meta <- meta[,c('Kids_First_Biospecimen_ID','composition.type','tumor_descriptor','primary_site','reported_gender','sample_id')]
  plyr::count(meta$composition.type) # 9 A, 7 S
  meta$label <- paste0(meta$sample_id,'_', tolower(meta$composition.type))
  
  # read expression data (stranded)
  dat <- readRDS('data/pbta-gene-expression-rsem-fpkm.stranded.rds')
  dat <- dat[,colnames(dat) %in% c('gene_id', as.character(meta$Kids_First_Biospecimen_ID))] # subset to meta file
  colnames(dat)[2:ncol(dat)] <- meta$label
  
  # reduce data
  dat <- dat[apply(dat[,-1], 1, function(x) !all(x==0)),] # remove all rows with only 0
  dat <- dat[grep('_PAR_', dat$gene_id, invert = T),] # remove PAR_Y chr genes
  
  # convert into a matrix of gene symbol x sample name
  expr.collapsed <- dat %>% 
    separate(gene_id, c("gene_id", "gene_symbol"), sep = "\\_", extra = "merge") %>%
    pivot_longer(-c(gene_id, gene_symbol), names_to = "sample_name", values_to = "fpkm") %>% 
    group_by(gene_id) %>%
    mutate(means = mean(fpkm)) %>% 
    group_by(gene_symbol) %>%
    filter(means == max(means)) %>% 
    dplyr::select(-gene_id) %>% unique() %>%
    pivot_wider(id_cols = "gene_symbol", names_from = "sample_name", values_from = "fpkm") %>%
    column_to_rownames(var="gene_symbol") %>%
    as.data.frame()
  
  # save input matrix (test data)
  save(expr.collapsed, file = 'data/input.RData')
}

# use predict function
if(file.exists('data/mRNA_StemScore.tsv')){
  print("Predicted scores ready")
} else {
  main.predict(fnSig = 'data/pcbc-stemsig.tsv', expr = 'data/input.RData', fnOut = "data/mRNA_StemScore.tsv")
}

# visualization
predicted.stemness.scores <- read.delim('data/mRNA_StemScore.tsv', header = F)
colnames(predicted.stemness.scores) <- c("sample_name", "stemness_score")
predicted.stemness.scores <- cbind(colsplit(predicted.stemness.scores$sample_name, pattern = "_", names = c("sample_name","type")), 
                                   stemness_score = predicted.stemness.scores$stemness_score)

# boxplot
p <- ggplot(predicted.stemness.scores, aes(x = type, y = stemness_score, alpha = 0.5, fill = type)) +
  geom_boxplot(outlier.shape = 21, outlier.fill = "white", outlier.color = "white") + ggtitle('OCLR predicted Stemness index\nAdherent vs Suspension Cell types') +
  theme_bw() + theme_Publication(base_size = 10) + ylab('Stemness Index') + xlab("") + 
  guides(alpha = F, fill = F) + geom_jitter(position=position_jitter(width=.25), shape = 21) +
  scale_x_discrete(labels = c("s" = "Suspension", 
                              "a" = "Adherent")) +
  scale_fill_manual(values = c("s" = "#F8766D",
                               "a" = "#00BFC4")) + 
  labs(fill = "Type")
ggsave(plot = p, filename = 'results/oclr_stemness_index.png', device = 'png', width = 5, height = 5)
