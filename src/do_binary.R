#######################################################################
# Calculate binarization from MCMC iterations of already saved files  #
# Perform calculations and reading in parallel                        #
# Same nomenclarture from La Manno, 2016                                  #
#######################################################################


binary.to.data.frame <- function(OUTPUT.folder=OUTPUT,basal.label="Sz",cores=4){
  require(doParallel)
  
  #Read METADATA
  load(file=paste0(OUTPUT.folder,"Genes_cell_types.Rdata",collapse = "")) 
  Cell_type <- metadata$Cell_type
  TGTGENE <- metadata$genes
  
  #Iterative conditional binarization of calculations in parallel
    cl <- makeCluster(cores)
    registerDoParallel(cl)
  
  
  BINS <- foreach(i=seq_along(TGTGENE),.packages = c("dplyr","matrixStats"),.combine = rbind) %dopar% {
    #Each gene MCMC permutation
    Gene <- TGTGENE[i] 
    INFILE <- paste0("OUTPUT/MAP_", Gene, ".rds", collapse="")
    MCMC <- (readRDS(INFILE) - 1)[,Cell_type]
    #get basal permutation, by label
    Basal <- MCMC[,basal.label]
    difference_basal <- MCMC[,Cell_type]-Basal+1e-12 #add 1e-12 to avoid zeroes
    
    #Condition for binarized n1
    TEMP.BIN <- data.frame()    
    bigger_condition <- colSums(difference_basal>0) >= 998 #For max strictness do == 1000
    #Check if gene is "basal"
    bigger_condition[basal.label] <- mean(difference_basal > 0) > 0.999
    TEMP.BIN <- rbind(TEMP.BIN,bigger_condition)
    
    #Condition for binarized n2
    median_categories <- colMedians(as.matrix(MCMC[,Cell_type])) 
    median_condition <- median_categories >= 0.35 * max(median_categories)
    TEMP.BIN <- rbind(TEMP.BIN,median_condition)
    
    #Condition for binarized n3
    high_enough <- colMaxs(as.matrix(MCMC[,Cell_type]))
    TEMP.BIN <- rbind(TEMP.BIN,high_enough) > 0.4
    
    #Logic all conditions must be TRUE, taking advantange of colProds default output
    TEMP.BIN <- colProds(as.matrix(TEMP.BIN))
    
    #Stored in binary matrix with gene name
    TEMP.BIN <- c(Gene,as.vector(TEMP.BIN)) %>% noquote()
    names(TEMP.BIN) <- c("Feature",colnames(MCMC[,Cell_type]))
    return(TEMP.BIN)
    
  }
  stopCluster(cl)
  BINS <- BINS %>% as.data.frame() 
  #drop rownames, they are not usefull
  rownames(BINS) <- c()
  return(BINS)
}
