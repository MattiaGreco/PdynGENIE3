# Tutorial on the package ROCR for ROC characteristics
library(ROCR)
data(ROCR.hiv)

predictions <- ROCR.hiv$hiv.svm$predictions  # a list, one can also have just one vector
labels <- ROCR.hiv$hiv.svm$labels    # a list, one can also have just one vector

pred <- prediction(predictions, labels)
pred
#> A prediction instance
#>   with 10 cross validation runs (equal lengths)

perf <- ROCR::performance(pred, "tpr", "fpr")
perf
plot(perf,
     avg= "threshold",
     colorize=TRUE,
     lwd= 3,
     main= "AUROC curve")
plot(perf,
     lty=3,
     col="grey78",
     add=TRUE)


perf <- ROCR::performance(pred, "prec", "rec")
plot(perf,
     avg= "threshold",
     colorize=TRUE,
     lwd= 3,
     main= "... Precision/Recall plot ...")
plot(perf,
     lty=3,
     col="grey78",
     add=TRUE)

demo(ROCR) # starts a demonstration of further graphical capabilities of ROCR. 
help(package=ROCR)


library(ROCR)
data(ROCR.simple)
reorder <- order(ROCR.simple$predictions)
ordered_data <- ROCR.simple$predictions[reorder]
ordered_labels <- ROCR.simple$labels[reorder]
pred <- prediction(ordered_data, ordered_labels)
pred

num_data <- length(reorder)
x <- 1:num_data
y <- pred@tp[[1]][2:(num_data+1)]/(pred@tp[[1]][2:(num_data+1)]+pred@fp[[1]][2:(num_data+1)])
plot(x,y)

perf <- ROCR::performance(pred, "prec", "rec")
perf
plot(perf,
     avg= "threshold",
     colorize=TRUE,
     lwd= 3,
     main= "AUROC curve")


