# Function to Predict using PCBC dataset
# Arguments:
  # fnSig - PCBC training data
  # expr - expression matrix with rownames collapse to HUGO genes and colnames are samples
  # fnOut - Test data output

expr <- 'data/input.RData'
fnSig <- 'data/pcbc-stemsig.tsv'

main.predict <- function(fnSig = "data/pcbc-stemsig.tsv", expr, fnOut = "data/mRNA_StemScore.tsv") {
  ## Load the signature
  w <- read.delim(fnSig, header = FALSE, row.names = 1) %>% 
    as.matrix() %>% 
    drop()
  
  ## load test dataset
  X <- get(load(expr))
  X <- X[rownames(X) %in% names(w),] # subset to training set genes
  
  ## Reduces HUGO|POSITION gene IDs to just HUGO
  # f <- function(v) unlist(lapply(strsplit(v, "\\|"), "[[", 1))
  
  # s <- synGet("syn4976369", downloadLocation = "/data/pancan")
  # X <- read.delim(s@filePath, as.is=TRUE, check.names = FALSE) %>%  ## Read the raw values
  #   filter( !grepl( "\\?", gene_id ) ) %>%      ## Drop genes with no mapping to HUGO
  #   mutate( gene_id = f( gene_id ) ) %>%        ## Clip gene ids to HUGO
  #   filter( gene_id %in% names(w) )         ## Reduce to the signature's gene set
  # 
  # ## SLC35E2 has multiple entries with the same HUGO id
  # ## Keep the first entry only
  # j <- grep( "SLC35E2", X[,1] )
  # if( length(j) > 1 )
  #   X <- X[-j[-1],]
  # 
  # ## Convert to a matrix
  # rownames(X) <- NULL
  # X <- X %>% tibble::column_to_rownames( "gene_id" ) %>% as.matrix()
  
  ## Reduce the signature to the common set of genes
  stopifnot( all( rownames(X) %in% names(w) ) )
  w <- w[rownames(X)]
  
  ####### Score via Spearman correlation
  s <- apply(X, 2, function(z) {cor( z, w, method = "sp", use = "complete.obs" )})
  
  ## Scale the scores to be between 0 and 1
  s <- s - min(s)
  s <- s / max(s)
  
  write.table(cbind(s), file = fnOut, sep = "\t", quote = FALSE, col.names = FALSE)
}
