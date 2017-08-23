source("src/Ceftools.R")
source("src/BayesianGLM.R")
source("src/Violinplot.R")
source("src/whole_dataset_BayesianGLM.R")
source("src/do_binary.R")
source("src/whole_dataset_median_DF.R")

# Whole dataset calculations

DATASET <- "DATA/Human_Embryo_fulldataset.cef"

# Run MCMC on whole dataset
# Files saved on folder OUTPUT


BayesianGLM_alldataset(DATASET = DATASET,parallel = T)


#Not run until end, test sample with 389 genes

# get genes analyzed
# Only for this example
library(dplyr)
library(stringr)
genes <- dir(path = "OUTPUT/",pattern = "MAP") %>% noquote()
genes <- data.frame(files=genes,stringsAsFactors = F)

genes2 <- data.frame(do.call('rbind', strsplit(as.character(genes$files),'AP_',fixed=TRUE))) 
genes2 <- genes2 %>% select(X2) %>% data.frame(do.call('rbind', strsplit(as.character(genes2$X2),'.',fixed=TRUE))) %>% 
  select(X1) %>% rename(Genes=X1)

TGTGENE <- as.character(genes2$Genes)


# functions for example dataset
Medians.to.data.frame <- function(OUTPUT.folder=OUTPUT,cores=4){
  require(dplyr)
  require(tibble)
  require(doParallel)
  #Read METADATA
  load(file=paste0(OUTPUT.folder,"Genes_cell_types.Rdata",collapse = "")) 
  Cell_type <- c("Sz",metadata$Cell_type) #Add basal
  #TGTGENE <- metadata$genes
  TGTGENE <- as.character(genes2$Genes)
  #Iterative median calculation in parallel
  cl <- makeCluster(cores)
  registerDoParallel(cl)
  #MED <- data.frame(row.names = Cell_type)

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

binary.to.data.frame <- function(OUTPUT.folder=OUTPUT,basal.label="Sz",cores=4){
  require(doParallel)
  
  #Read METADATA
  load(file=paste0(OUTPUT.folder,"Genes_cell_types.Rdata",collapse = "")) 
  Cell_type <- c("Sz",metadata$Cell_type) #Add basal
  #Cell_type <- metadata$Cell_type
  #TGTGENE <- metadata$genes
  #TGTGENE <- metadata$genes
  TGTGENE <- as.character(genes2$Genes)
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
  rownames(BINS) <- c()
  
  return(BINS)
}




#Test Cells
medians.cells <- Medians.to.data.frame(OUTPUT.folder = "OUTPUT/",cores = 8)


binary.cells <- binary.to.data.frame(OUTPUT.folder = "OUTPUT/",basal.label = "Sz",cores = 8)

