#################################################################
# Calculate medians from MCMC iterations of already saved files #
# Perform calculations and reading in parallel                  #
#################################################################

# This is just wrong
# Medians.to.data.frame <- function(OUTPUT.folder=OUTPUT,cores=4){
#   require(dplyr)
#   require(tibble)
#   require(doParallel)
#   
#   #Read METADATA
#   load(file=paste0(OUTPUT.folder,"Genes_cell_types.Rdata",collapse = "")) 
#   Cell_type <- metadata$Cell_type
#   TGTGENE <- metadata$genes
#   
#   #Iterative median calculation in parallel
#   
#   cl <- makeCluster(cores)
#   
#   registerDoParallel(cl)
#   
#   MED <- data.frame(row.names = Cell_type)
#   MEDS <- foreach(i=seq_along(TGTGENE),.packages = "dplyr",.combine = cbind) %dopar% {
#     Gene <- TGTGENE[i] 
#     INFILE <- paste0("OUTPUT/MAP_", Gene, ".rds", collapse="")
#     BETA <- (readRDS(INFILE) - 1)[,Cell_type]
#     TEMP <- apply(BETA,MARGIN = 2,FUN = median)
#     MED  <- data.frame(MED,TEMP) %>% rename_(., .dots = setNames("TEMP", Gene))
#     return(MED)
#     }
#     
#   stopCluster(cl)
#   
#   MEDS <- MEDS %>% t() %>% as.data.frame() %>% rownames_to_column("Features")
#   return(MEDS)
# }

Medians.to.data.frame <- function(OUTPUT.folder=OUTPUT,cores=8){
  require(dplyr)
  require(tibble)
  require(doParallel)
  #Read METADATA
  load(file=paste0(OUTPUT.folder,"Genes_cell_types.Rdata",collapse = "")) 
  Cell_type <- c("Sz",metadata$Cell_type) #Add basal
  TGTGENE <- as.character(metadata$genes)
  #Iterative median calculation in parallel
  cl <- makeCluster(cores)
  registerDoParallel(cl)
  
  MEDS <- foreach(i=seq_along(TGTGENE),.packages = "dplyr",.combine = cbind) %dopar% {
    Gene <- TGTGENE[i] 
    INFILE <- paste0("OUTPUT/MAP_", Gene, ".rds", collapse="")
    BETA <- (readRDS(INFILE) - 1)[,Cell_type]
    TEMP <- apply(BETA,MARGIN = 2,FUN = median) %>% as.data.frame()
    colnames(TEMP) <- Gene
    return(TEMP)
    }
  
  stopCluster(cl)
  
  MEDS <- MEDS %>% t() %>% as.data.frame() %>% rownames_to_column("Features")
  return(MEDS)
}