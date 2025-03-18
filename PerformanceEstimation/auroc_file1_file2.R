# Get AUROC when comparing two files (2 predictions or gold vs prediction)
# usage from a shell: Rscript auroc_file1_file2.R file1 file2
# where the paths are absolute or relative (but do not start with / nor ~, instead ./ should be used)
# calling directly the R from the command line messes up the arguments as they start with . or / or ~

#############################################################################
#                     Recover PARAMETERS passed                             #
#############################################################################
args <- commandArgs(trailingOnly = TRUE)
#  print(args)

if(length(args)!=2){
    print("Usage: Rscript auroc_file1_file2.R file1 file2")
    exit()
}else{
    file1<-args[[1]]
    file2<-args[[2]]
}

file1 <- "./Goldstandards/goldSigned_1.tsv"
file2 <- './Predictions/GENIE3'
print(file1)
print(file2)



#############################################################################
#                     Function for data manipulation                        #
#############################################################################

make_dataframe_for_auroc=function(gold,prediction) {
# Make a global dataframe combining the gold and prediction information
# It is convenient to creates the row names for the underlying dataframes
# Careful: importance lists are no longer ordered!
  these_row_names <- paste(gold[,1],"_",gold[,2],sep="")
  row.names(gold) <- these_row_names
  these_row_names <- paste(prediction[,1],"_",prediction[,2],sep="")
  row.names(prediction) <- these_row_names

  all_names <- unique(c(row.names(prediction),row.names(gold)))
  from_and_to <- stringr::str_split(all_names,pattern="_")
  num_arcs <- length(all_names)
  from <- rep(NA,num_arcs)
  to <- from
  for(ii in 1:num_arcs) {
    from[ii] <- from_and_to[[ii]][1]
    to[ii] <- from_and_to[[ii]][2]
  }
  prediction_gold <- data.frame(from=from,to=to,importance_prediction=rep(0,num_arcs),importance_gold=rep(0,num_arcs))
  row.names(prediction_gold) <- all_names
  prediction_gold[row.names(prediction),]$importance_prediction <- prediction[row.names(prediction),3]
  prediction_gold[row.names(gold),]$importance_gold <- gold[row.names(gold),3]
  return(prediction_gold)
}


#############################################################################
#                 AUROC and AUPRC between the 2 files                       #
# Just interpret the gold numbers as the probability of being TRUE          #
# AUROC: order on the x axis according to decreasing weight(prediction)     #
#        x=fraction_of_total_(1-weights)(predictions)                       #
#        y=fraction_of_total_weights(gold)                                  #
# AUPRC: 
#############################################################################


  this_gold <- read.table(file=file1,header=FALSE) # real numbers or only 0s and 1s
  this_prediction <- read.table(file=file2,sep="",header=FALSE)
  my_df <- make_dataframe_for_auroc(this_gold,this_prediction)
  # that is sorted according to importance_prediction
  print(paste("spearman correlation coefficient",cor(my_df$importance_prediction,my_df$importance_gold,method="spearman")))
  my_test <- cor.test(my_df$importance_prediction,my_df$importance_gold,method="spearman")
  print(paste("p.value",my_test$p.value))
  quit()


  FPR <- cumsum(my_df$importance_prediction)
  my_df <- cbind(my_df,fpr=FPR)
  # limma::auROC(my_df$importance_gold,my_df$importance_prediction) 
  # Can also get different curves with the package library(ROCR)
  pred <- ROCR::prediction(my_df$importance_prediction,my_df$importance_gold)
  my_auroc <- ROCR::performance(pred,"auc",fpr.stop=1.0)
  print(paste("auroc = ",my_auroc@y.values[[1]]))
  my_auprc <- ROCR::performance(pred,"aucpr",fpr.stop=1.0)
  print(paste("auprc = ",my_auprc@y.values[[1]]))



