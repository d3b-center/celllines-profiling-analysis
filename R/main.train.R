# Function to Train using PCBC dataset
# Arguments:
  # fnOut - filename of the output signature
  # fnGenes - [optional] filename of the list of entrez ID to consider

main.train <- function(fnOut = "data/pcbc-stemsig.tsv", fnGenes = NULL) {
  # Train (get PCBC dataset and save locally)
  synRNA <- synGet("syn2701943", downloadLocation = "data/PCBC" )
  X <- read.delim(synRNA$path) %>%
    tibble::column_to_rownames( "tracking_id" ) %>% 
    as.matrix
  
  ## Retrieve metadata
  synMeta <- synTableQuery( "SELECT UID, Diffname_short FROM syn3156503" )
  Y <- synMeta %>%
    as.data.frame(synMeta) %>%
    mutate(UID = gsub("-", ".", UID)) %>%
    dplyr::select(Diffname_short, UID) %>%
    tibble::column_to_rownames("UID")
  
  # Subset Y using X
  y <- Y[colnames(X),]
  names(y) <- colnames(X)
  
  ## Fix the missing labels by hand
  y["SC11.014BEB.133.5.6.11"] <- "EB"
  y["SC12.039ECTO.420.436.92.16"] <- "ECTO"
  
  ## Drop the splice form ID from the gene names
  v <- strsplit( rownames(X), "\\." ) %>% lapply( "[[", 1 ) %>% unlist()
  rownames(X) <- v
  
  ## Map Ensembl IDs to HUGO
  V <- genes2hugo(rownames(X))
  
  # change rownames to hugo symbols
  X <- X[V[,1],]
  rownames(X) <- V[,2]
  
  # Reduce gene set to the provides list. (no input yet)
  if(exists('fnGenes') && !is.null(fnGenes)){
    vGenes <- read.delim(fnGenes, header=FALSE) %>% as.matrix() %>% drop()
    VE <- genes2hugo( vGenes, "entrezgene" )
    X <- X[intersect( rownames(X), VE[,2] ),]
  } else {
    print("Do nothing")
  }
  
  # Mean-center the data
  m <- apply(X, 1, mean)
  X <- X - m
  
  # Identify stem cell samples
  # Break up all samples into 2 groups:
  j <- which( y == "SC" )
  X.tr <- X[,j] # stem cell
  X.bk <- X[,-j] # non stem-cell
  
  # Train a one-class model
  mm <- gelnet(X = t(X.tr), y = NULL, l1 = 0, l2 = 1)
  
  ## Store the signature to a file
  write.table(mm$w, file = fnOut, sep = "\t", quote = FALSE, col.names = FALSE)
  
  # Perform leave-one-out cross-validation
  auc <- c()
  for(i in 1:ncol(X.tr)){
    ## Train a model on non-left-out data
    X1 <- X.tr[,-i]
    m1 <- gelnet( t(X1), NULL, 0, 1 )
    
    ## Score the left-out sample against the background
    s.bk <- apply( X.bk, 2, function(z) {cor( m1$w, z, method="sp" )} )
    s1 <- cor( m1$w, X.tr[,i], method="sp" )
    
    ## AUC = P( left-out sample is scored above the background )
    auc[i] <- sum( s1 > s.bk ) / length(s.bk)
    cat( "Current AUC: ", auc[i], "\n" )
    cat( "Average AUC: ", mean(auc), "\n" )
  }
  return(auc)
}
