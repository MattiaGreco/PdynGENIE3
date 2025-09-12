# Example for running PdynGENIE3 on a few test examples
# Compared to dynGENIE3, it provides priors based on time-lagged correlation
# functions and estimates better decay constants
# See also https://github.com/MattiaGreco/PdynGENIE3
# To run this on your computer, set these paths to their correct value:
# path_source  path_input  path_pred

# Libraries needed
# library("reshape2")
# library(doRNG) if use multicore

# Prefix of the files with golds and time series
nomenclature="goldSigned_"

# Path for R scripts
path_source <- "~/LOCAL/GRN_Inference/BENCHMARKING/PdynGENIE3/Source_PdynGENIE3"

# Path for input data (time series expression)
path_input <- "~/LOCAL/GRN_Inference/Mattia_TEST/In_out_degree_1_num_nodes_100/Simulation_files_Time_Series/"
list_files <- list.files(path_input) # gives us the names of the files under the given path

# Path to output predictions (importances as a file "from" "to" "strength")
path_pred <- "~/LOCAL/GRN_Inference/Mattia_TEST/In_out_degree_1_num_nodes_100/Predictions_Time_Series/"

setwd(path_source)
# Source "PdynGENIE3.R" to get main functions
# Note: the library object file PdynGENIE3.so should be here too
source("PdynGENIE3.R")

# Source "UsefulFunctions.R" for other useful functions
source("UsefulFunctions.R")

num <- numeration_of_list_files(list_files)
ii <- 1
for (ii in 1:length(num)){
  
  print("######################################")
  print(paste("graph number",num[ii]))
  print("######################################")
  path_data=paste(path_input,list_files[ii],sep="")
  this_graph_error <- FALSE

  df <- read.expr.matrix(path_data,form="rows.are.samples")
  # row 1 is time
  # following rows are genes
  # columns are samples
  
  # Get time series data
  # This example is of the DREAM type: perturbed basal transcription rate for first half of series
  # The time series have 21 time points and we have 10 experiments
  TS.data <- list()
  time.points <- list()
  num_time_points <- 21
  start_col <- 1
  end_col <- num_time_points
  for(i_time_series in 1:10) {
    TS.data <- c(TS.data,list(df[2:nrow(df),start_col:end_col]))
    time.points <- c(time.points,list(df[1,start_col:end_col]))
    start_col <- start_col+num_time_points
    end_col <- end_col+num_time_points
  }
  
  regulator_genes <- row.names(df)[2:101]  # Default: all genes are putative regulators

  # Compute Correlations and decay rates estimates
  print("Compute correlation matrix and decay rates")
  theory_estimates <- different_time_corr(TS.data, time.points)
  str(theory_estimates)
  theory_priors <- abs(theory_estimates$Lambda_target_i_driver_j)
  theory_alphas <- pmax(theory_estimates$alphas,0)  # in case some values are negative
  names(theory_alphas) <- row.names(df)[2:101]
  outfile <- paste0(path_pred,nomenclature,num[ii],"_alphas.tsv")
  # write.table(theory_alphas,file=outfile,row.names=TRUE,col.names=FALSE,quote=FALSE,sep="\t")
  # next # uncomment if only want to get the alphas
  ncores <- 8  # The number of cores of your computer you can use
 
  res <- PdynGENIE3(TS.data, time.points, alpha = theory_alphas, regulators = regulator_genes, ntrees = 100,
                K = "sqrt", ncores = ncores, seed=12345, verbose = TRUE, priors = theory_priors)
  weights <- res$weight.matrix
  
  diag(weights) <- 0  # Remove self interactions. If one is dealing with a rectangular matrix one must use the gene names.
  colnames(weights) <- regulator_genes
  rownames(weights) <- regulator_genes
  
  # write predictions to file
  pred <- get.link.list(weights,threshold=0.0)
  colnames(pred) <- c("from","to","strength") # to standardize the notations
  outfile <- paste0(path_pred,nomenclature,num[ii],"_prediction.tsv")
  write.table(pred,file=outfile,row.names=FALSE,col.names=FALSE,sep="\t")
}
