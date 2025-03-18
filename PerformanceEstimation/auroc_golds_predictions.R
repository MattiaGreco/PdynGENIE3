# Get AUROC averaged over all the graphs in their 2 directories (goldSigned and predicted)
# usage from a shell: Rscript auroc_golds_predictions.R golds_directory predictions_directory
# where the paths are absolute or relative (but do not start with / nor ~, instead ./ should be used)
# and do not end with / as that is added by the code.
# calling directly the R from the command line messes up the arguments as they start with . or / or ~

#Put your paths here
golds_directory <-  "~/Desktop/ComputerVecchio/PaperCodes/GNWData/In_out_degree_4_num_nodes_100/Goldstandards"  
predictions_directory <- "~/Desktop/ComputerVecchio/PaperCodes/GNWData/In_out_degree_4_num_nodes_100/Predictions/Temperature/beta_5"      


#############################################################################
#                     Info for the files                                    #
#############################################################################
if(!("stringr" %in% installed.packages())) {
  install.packages("stringr")
}
library(stringr)

numbers_for_list_files <- function(list_files){
  # retrieves the numbers (graph number) from the files in the
  # given directory (e.g., 1  10 100  11  12  13...)
  #
  # Parameters(required): 
  #   list_files = output of list_files: the names of the files in the directory
  #
  # Returns:
  #   a dataframe containing
  #    - the entry in list_files of the file (e.g., 1   2   3   4   5   6...)
  #    - the numbers in the name of the file (e.g., 1  10 100  11  12  13...)
  #
  l <- length(list_files)
  num <- data.frame(order_appearance=seq(1,l,by=1),graph_number=rep(0,l))
  for (i in 1:l){
    split1 <- str_split(list_files[i],pattern="\\.")
    split2 <- str_split(split1[[1]][1],pattern="_")
    num[i,]$graph_number <- as.numeric(split2[[1]][2])
  }
  return(num)
}

golds_list_files <- list.files(golds_directory)
golds_numbers <- numbers_for_list_files(golds_list_files)
predictions_list_files <- list.files(predictions_directory)
predictions_numbers <- numbers_for_list_files(predictions_list_files)
valid_graphs_lines <- golds_numbers$graph_number %in% predictions_numbers$graph_number
valid_numbers <- golds_numbers[valid_graphs_lines,]

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
  prediction_gold <- data.frame(from=from,to=to,importance_prediction=rep(0,num_arcs),gold=rep(0,num_arcs))
  row.names(prediction_gold) <- all_names
  prediction_gold[row.names(prediction),]$importance_prediction <- prediction[row.names(prediction),3]
  
  lines <- which((gold[,3] == "+") | (gold[,3] == "-"))
  gold[lines, 3] <- 1
  
  prediction_gold[row.names(gold),]$gold <- gold[row.names(gold),3]
  return(prediction_gold)
}


#############################################################################
#                  Mean, sd and se of AUROC across the graphs               #
#############################################################################
# library(limma)
library(ROCR)

all_AUROCs_AUPRCs <- data.frame(graph_num=valid_numbers$graph_number,
                    auroc=rep(NA,dim(valid_numbers)[1]),auprc=rep(NA,dim(valid_numbers)[1]))

pr1 <- matrix(0, nrow = 1, ncol=400)
ii <-1

#Computation of average AUROC and AUPRC across all graphs in a folder
for(ii in 1:dim(valid_numbers)[1]) {
  this_graph_num <- valid_numbers[ii,]$graph_number
  all_AUROCs_AUPRCs[ii,]$graph_num <- this_graph_num
  jj <-  valid_numbers[ii,]$order_appearance
  this_gold <- read.table(file=paste(golds_directory,golds_list_files[jj],sep="/"),header=FALSE) # Only 0s and 1s
  kk <- predictions_numbers[which(predictions_numbers$graph_number == this_graph_num),]$order_appearance
  this_prediction <- read.table(file=paste(predictions_directory,predictions_list_files[kk],sep="/"),header=FALSE)
  my_df <- make_dataframe_for_auroc(this_gold,this_prediction)
  # all_AUROCs_AUPRCs[ii,]$auroc <- limma::auROC(my_df$gold,my_df$importance_prediction) 
  # Can also get different curves with the package library(ROCR)
  pred <- ROCR::prediction(my_df$importance_prediction,my_df$gold)
  my_auroc <- ROCR::performance(pred,"auc",fpr.stop=1.0)
  all_AUROCs_AUPRCs[ii,]$auroc <- my_auroc@y.values[[1]]
  my_auprc <- ROCR::performance(pred,"aucpr",fpr.stop=1.0)
  all_AUROCs_AUPRCs[ii,]$auprc <- my_auprc@y.values[[1]]
  
  my_auprc <- ROCR::performance(pred,"prec","tpr")
  pr1 <- pr1 + my_auprc@y.values[[1]][2:401]
  #plot(my_auprc)
}

print(paste("Mean_AUROC = ",mean(all_AUROCs_AUPRCs$auroc),"sd AUROC = ",sd(all_AUROCs_AUPRCs$auroc),"se AUROC = ",sd(all_AUROCs_AUPRCs$auroc)/sqrt(dim(valid_numbers)[1])))
print(paste("Mean_AUPRC = ",mean(all_AUROCs_AUPRCs$auprc),"sd AUPRC = ",sd(all_AUROCs_AUPRCs$auprc),"se AUPRC = ",sd(all_AUROCs_AUPRCs$auprc)/sqrt(dim(valid_numbers)[1])))

### FOR FIGURE 10
b <- c(0, 1, 2, 3, 4, 5)

####Results for In-out=1
R <- c(0.890610897959178, 0.920611239795917, 0.920237489795911,  0.920962280612236, 0.914202448979586, 0.909787658163264 )
P <- c(0.323081448705781, 0.621060097398035, 0.652697781352994,  0.711954023335683, 0.664338333807382, 0.660918903349862 )

plot(b, R, main='AUROC in-out=1', xlab='1/T', ylab='AUROC', cex.axis = 1.2, cex.lab =1.5, cex.main=2, lwd = 5, col='red' , type = 'b', pch = 5)
plot(b, P, main='AUPRC in-out=1', xlab='1/T', ylab='AUPRC', cex.axis = 1.2, cex.lab =1.5, cex.main=2, lwd = 5, col='red', type = 'b', pch = 5)

####Results for In-out=2
R <- c(0.742373453608247, 0.784128121134022, 0.78482411597937 , 0.78, 0.773427335051545, 0.76920655154639)
P <- c(0.135505169524548, 0.29689189707259, 0.313190624044092 , 0.36, 0.325240852131572, 0.324266664386999)

plot(b, R, main='AUROC in-out=2', xlab='1/T', ylab='AUROC', cex.axis = 1.2, cex.lab =1.5, cex.main=2, lwd = 5, col='red' , type = 'b', pch = 5)
plot(b, P, main='AUPRC in-out=2', xlab='1/T', ylab='AUPRC', cex.axis = 1.2, cex.lab =1.5, cex.main=2, lwd = 5, col='red', type = 'b', pch = 5)

####Results for In-out=3
R <- c(0.69663090451389, 0.70258682986111, 0.70099951388889, 0.695736654513889, 0.690265803819444, 0.685774276041667 )
P <- c(0.16726985961397, 0.19649762790480, 0.21561767818141, 0.220344854788354, 0.220231777357555, 0.218937162939258 )

plot(b, R, main='AUROC in-out=3', xlab='1/T', ylab='AUROC', cex.axis = 1.2, cex.lab =1.5, cex.main=2, lwd = 5,col='red' , type = 'b', pch = 5)
plot(b, P, main='AUPRC in-out=3', xlab='1/T', ylab='AUPRC', cex.axis = 1.2, cex.lab =1.5, cex.main=2, lwd = 5, col='red', type = 'b', pch = 5)

####Results for In-out=4
R <- c(0.65124685394737, 0.654927953947368, 0.652913567105264, 0.647331585526316, 0.642654881578947, 0.639761592105263 )
P <- c(0.14327251619766, 0.164297846730227, 0.180164557328066, 0.182345043192405, 0.177427293387241, 0.176764613686765 )

plot(b, R, main='AUROC in-out=4', xlab='1/T', ylab='AUROC',col='red' , cex.axis = 1.2, cex.lab =1.5, cex.main=2, lwd = 5, type = 'b', pch = 5)
plot(b, P, main='AUPRC in-out=4', xlab='1/T', ylab='AUPRC', col='red', cex.axis = 1.2, cex.lab =1.5, cex.main=2, lwd = 5, type = 'b', pch = 5)








plot(as.vector(prN[2:101])/100, col = rgb(1, 0, 0, alpha = 0.75), ylim=c(0.4, 1), pch = 19, xlab="Interactions", ylab="PR")
points(as.vector(prA[2:101])/100, col = rgb(0, 1, 0, alpha = 0.75), pch = 19)
points(as.vector(prOriginal[2:101])/100, col = rgb(0, 0, 1, alpha = 0.75), pch = 19)

legend("topright", legend = c("dynGENIE3", "Priors", "Priors + alphas"),
       col = c(rgb(0, 0, 1, alpha = 0.75), rgb(0, 1, 0, alpha = 0.75),rgb(1, 0, 0, alpha = 0.75)) , pch = 19)
title(main='Average PR curves over 100 graphs In-out=1')


prN <- pr1
prA <- pr1
prO <- pr1

plot(1:200, my_auprc@y.values[[1]][1:200])
plot(1:8000, my_auprc@y.values[[1]][1:8000])
plot(1:length(my_auprc@y.values[[1]]), my_auprc@y.values[[1]])

PR <- cbind(as.vector(prN),as.vector(prA),as.vector(prO))

write.table( PR,paste("~/Desktop/PaperCodes/GNWData/In_out_degree_4_num_nodes_100/Predictions/",
                              "PRavg.txt",sep=""),row.names=FALSE,col.names=FALSE,sep="\t")





ii <-1
for(ii in 1:dim(valid_numbers)[1]) {
  this_graph_num <- valid_numbers[ii,]$graph_number
  all_AUROCs_AUPRCs[ii,]$graph_num <- this_graph_num
  jj <- which(golds_list_files=="goldSigned_1.tsv")#valid_numbers[ii,]$order_appearance
  this_gold <- read.table(file=paste(golds_directory,golds_list_files[jj],sep="/"),header=FALSE) # Only 0s and 1s
  kk <- predictions_numbers[which(predictions_numbers$graph_number == this_graph_num),]$order_appearance
  this_prediction <- read.table(file=paste(predictions_directory,predictions_list_files[kk],sep="/"),header=FALSE)
  my_df <- make_dataframe_for_auroc(this_gold,this_prediction)
  # all_AUROCs_AUPRCs[ii,]$auroc <- limma::auROC(my_df$gold,my_df$importance_prediction) 
  # Can also get different curves with the package library(ROCR)
  pred <- ROCR::prediction(my_df$importance_prediction,my_df$gold)
  my_auroc <- ROCR::performance(pred,"auc",fpr.stop=1.0)
  all_AUROCs_AUPRCs[ii,]$auroc <- my_auroc@y.values[[1]]
  my_auprc <- ROCR::performance(pred,"aucpr",fpr.stop=1.0)
  all_AUROCs_AUPRCs[ii,]$auprc <- my_auprc@y.values[[1]]
}
print(paste("Mean_AUROC = ",mean(all_AUROCs_AUPRCs$auroc),"sd AUROC = ",sd(all_AUROCs_AUPRCs$auroc),"se AUROC = ",sd(all_AUROCs_AUPRCs$auroc)/sqrt(dim(valid_numbers)[1])))
print(paste("Mean_AUPRC = ",mean(all_AUROCs_AUPRCs$auprc),"sd AUPRC = ",sd(all_AUROCs_AUPRCs$auprc),"se AUPRC = ",sd(all_AUROCs_AUPRCs$auprc)/sqrt(dim(valid_numbers)[1])))

v <- all_AUROCs_AUPRCs[order(all_AUROCs_AUPRCs$graph_num),]

plot(beta, v$auroc, xlab = "1/T", ylab = "AUROC" )
abline(h = v$auroc[3])

plot(beta, v$auprc, xlab = "1/T", ylab = "AUPRC" )
abline(h = v$auprc[3])

plot(v$auroc, e1)






beta <- c(0.5, 1, 2, 3, 4, 5, 6, 7)
a00 <- all_AUROCs_AUPRCs

AUROC <-cbind(au0$auroc, au1$auroc, au2$auroc, au3$auroc, au4$auroc, au5$auroc, au6$auroc, au7$auroc)
AUPRC <- cbind(au0$auprc, au1$auprc, au2$auprc, au3$auprc, au4$auprc, au5$auprc, au6$auprc, au7$auprc)

aa <- c(mean(au0$auroc/au1$auroc), 1, mean(au2$auroc/au1$auroc), mean(au3$auroc/au1$auroc),
        mean(au4$auroc/au1$auroc), mean(au5$auroc/au1$auroc), mean(au6$auroc/au1$auroc), mean(au7$auroc/au1$auroc))

plot(beta, aa, xlab = "1/T", ylab = "AUROC(1/T) / AUROC(1)", main = "AUROC in-out=4", pch =17, col="red")

aa <- c(mean(au0$auprc/au1$auprc), 1, mean(au2$auprc/au1$auprc), mean(au3$auprc/au1$auprc),
        mean(au4$auprc/au1$auprc), mean(au5$auprc/au1$auprc), mean(au6$auprc/au1$auprc), mean(au7$auprc/au1$auprc))

plot(beta, aa, xlab = "1/T", ylab = "AUPRC(1/T) / AUPRC(1)", main = "AUPRC in-out=3", pch =17, col="green")



sroc <- c(sd(au0$auroc), sd(au1$auroc), sd(au2$auroc), sd(au3$auroc), sd(au4$auroc), sd(au5$auroc))
beta <- c(0.5, 1, 2, 3, 4, 5)


plot(beta, colSums(AUROC)/100,  xlab = '1/T', ylab="AUROC", 
     main = "AUROC vs 1/T, in-out=4", pch = 17, col ="red")
abline(h = mean(auf$auroc))

plot(beta, colSums(AUPRC)/100,ylim = c(0, 0.2),  xlab = '1/T', ylab="AUPRC",
     main = "AUPRC vs 1/T, in-out=4", pch = 17, col ="red")

abline(h = mean(auf$auprc))



plot(au0$auprc, au$auprc, xlim = c(0., 0.2), ylim = c(0., 0.25), xlab = "No priors", ylab = "priors",
     main = "AUPRC In-Out=4")
abline(0,1)


