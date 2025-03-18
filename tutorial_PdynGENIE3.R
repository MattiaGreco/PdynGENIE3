#Source of the dynGENIE3 code with priors 
source("~/Desktop/PdynGENIE3/PdynGENIE3_sourceCode/PdynGENIE3.R")
#Source of the code with a lot of useful functions
source("~/Desktop/PdynGENIE3/UsefulFunctions.R")

setwd("~/Desktop/PdynGENIE3")

#Name of the files with golds and time series
nomenclature="goldSigned_"

#Paths for IO, easiest thing is to do is just change the number of in-out degree
path <- "~/Desktop/PdynGENIE3/TestData_InOut1/"
path_simu <- paste(path,"Simulation_files_Time_Series/",sep="")
list_files <- list.files(path_simu) # gives us the names of the files under the given path
# we indicate the directory under which we wish to save the different predictions

path_pred <- paste(path,"Predictions/",sep="")
ex_file <- read.table(paste(path_simu,list_files[1],sep=""),header=TRUE,sep="\t")
num <- numeration_of_list_files(list_files)

for (i in 1:5){
  
  print("######################################")
  print(paste("graph number",i))
  print("######################################")
  i_this_listed_file <- which(num == i)
  path_data=paste(path_simu,list_files[i_this_listed_file],sep="")
  data = read.table(path_data,header=TRUE,sep="\t") # rows are samples, columns are genes
  this_graph_error <- FALSE

  df <- read.expr.matrix(path_data,form="rows.are.samples")
  # df has: rows are genes, columns are samples
  
  #Get time series data, in this case the time series have 21 time points and we have 10 expriments 
  TS.data <- list(df[2:nrow(df),1:21], df[2:nrow(df),22:42], df[2:nrow(df),43:63],  df[2:nrow(df),64:84], df[2:nrow(df),85:105], 
                  df[2:nrow(df),106:126], df[2:nrow(df),127:147], df[2:nrow(df),148:168], df[2:nrow(df),169:189], df[2:nrow(df),190:210])
  
  time.points <- list(df[1,1:21], df[1,1:21], df[1,1:21], df[1,1:21], df[1,1:21], df[1,1:21], df[1,1:21], df[1,1:21], df[1,1:21], df[1,1:21])
  regulator_genes <- row.names(df)[2:101]

  #Compute Correlations and decay rates estimates
  print("Compute correlation matrix and decay rates")
  C <- different_time_corr(TS.data, time.points)

  new_alphas <- C$alphas
  names(new_alphas) <- row.names(df)[2:101]
  
  res <- dynGENIE3(TS.data, time.points, alpha = new_alphas, regulators = regulator_genes, ntrees = 100,
                   K = 10, ncores = 1, seed=12345, verbose = TRUE, priors = abs(C$C))
  weights <- res$weight.matrix
  diag(weights) <- 0
  colnames(weights) <- regulator_genes
  rownames(weights) <- regulator_genes

}

