
#Function that performs LTR, needs two matrices weights and edges, plus a numerical threshold.

edge_sign_estimate <- function(corr, new_weights){

  num.genes <- nrow(corr)
  edges <- matrix(0, nrow = num.genes, ncol = num.genes)
  
  edges <- corr*new_weights
  
  edges[edges > 0] <- 1
  edges[edges < 0] <- -1
  
  
  return(edges)
  
}

library('lmtest')
corr_mat <- function(TS.data, time.points){
  
  num.exp <- length(TS.data)
  gene.names <- rownames(TS.data[[1]])
  num.genes <- length(gene.names)
  
  edges <- matrix(0, nrow = num.genes, ncol = num.genes)
  colnames(edges) <- gene.names
  row.names(edges) <- gene.names

  for(k in 1:num.exp){
    
    current.timeseries <- TS.data[[k]]
    current.time.points <- time.points[[k]]
    
    t <- length(time.points[[k]])
    print(k)
    for(j in 1:num.genes){
      
      for(i in 1:num.genes){
        edges[i, j] <- edges[i, j] + cor(current.timeseries[j,], current.timeseries[i,])
        
       
      }
    }
  }
  
  #set diagonal to zero and normalize
  diag(edges) <- 0
  
  edges <- edges/t
  
  return(edges)
  
}


LTR <- function(edges, new_weights, corr, beta=1){
  
  #Size of perturbation graph
  #n is the number of drivers thus number of lines in edges
  #n <- length(edges[,1])
  
  colnames(edges) <- colnames(new_weights)
  row.names(edges) <- row.names(new_weights)
  #Find edges present from i to j
  E_itoj <- abs(edges)   #E_ij =1 means interaction i to j
  
  #Find non edges i to j
  E_not_itoj <- abs(!(edges))
  
  #Initialize adjacency matrix of removable edges
  RE <- matrix(0, nrow = nrow(edges), ncol = ncol(edges))

  start <- Sys.time()  

  count <-0
  
  corr <- abs(corr)
  order_to_test <- order(corr)
  
  nrows <- nrow(corr)
  ncols <- ncol(corr)
  
  # --------for(i in 1:nrows){----------------
    
  #  j <- which(E_itoj[i,] != 0)    # & E_not_itoj[,i])
    
  #  l1 <- which(E_itoj[i,] & E_itoj[j,])
  #  l2 <- which((edges[i, j] == edges[i,]*edges[j,]))
  #  l3 <- which((corr[i, j] >= a*(corr[i,])/(corr[j,])))
    
  #  K <- intersect(l1, l2)
  #  K <- intersect(K, l3)
    
   
  #  RE[K] <-1
  #}
  
  #n <- length(which(RE == 1))
  #print("Edges removed:")
  #-----print(n)----
  flag <-0
  
    for (index in 1:length(order_to_test)){
  
      driver <- (order_to_test[index]-1)%%nrows +1
      target <- floor((order_to_test[index]-1)/nrows)+1
    
    if(corr[driver, target] == 0){
      next
    }
      #all possible js connecting i and k making a triangle ijk
      l1 <- which(E_itoj[driver,] & E_itoj[,target])   
      l2 <- which(edges[driver,target] == edges[driver,]*edges[,target])
      #l3 <- which(E_not_itoj[target, driver] & E_itoj[driver, target])
      
      l <- intersect(l1, l2)
      
      if(length(l) == 0){
        next
      }
      
      weights_triangle_edges <- matrix(0, nrow = length(l), ncol = 2)
      weights_triangle_edges[,1] <- corr[driver,l]**beta
      weights_triangle_edges[,2] <- corr[l,target]**beta
      
      #threshold <- max(weights_triangle_edges[,1]*weights_triangle_edges[,2])
      
      
      
      for(i in 1:length(l)){
        threshold <- min(weights_triangle_edges[i,1],weights_triangle_edges[i,2])
        
        if((corr[driver, target] < threshold) ){
          new_weights[driver, target] <- new_weights[driver, target]*corr[driver, target]**beta/max(weights_triangle_edges[i,1],weights_triangle_edges[i,2])
          #*corr[driver,target]/threshold
          #*threshold/(weights_triangle_edges[i,1]+weights_triangle_edges[i,2])
          flag <- flag+1
        }#else{
          #new_weights[]
          #new_weights[which(new_weights == threshold, arr.ind = TRUE)] <- threshold**2/(weights_triangle_edges[i,1]+weights_triangle_edges[i,2])
         #   }
        
      
      
      #if(flag == TRUE){
      #  count <- count+1
      #}
      }
    }
  #    if((new_weights[driver, target] < threshold) ){
        
  #      new_weights[driver, target] <- 0#new_weights[driver, target]**2/threshold
  #      flag <- flag+1
  #    }
      
    #}
  
    print("Edges touched:")
    print(flag)
  
  end <- Sys.time()
  elapsed <- end - start
  print(elapsed)
  
  #new_weights <- transitive_reduction(edges, RE, new_weights)
 
  return(new_weights)
}

transitive_reduction <- function(edges, RE, new_weights){
  
  #M <- abs(edges)
  
  indices <- which(RE != 0, arr.ind = TRUE)
  RE1 <- indices[, 1]
  RE2 <- indices[, 2]
  
  #M[RE1, RE2] <-  0
  
  #for (i in 1:length(RE1)) {
  #  if (!any(M[RE1[i],] & M[RE2[i], ])) {
  #    # Reinsert edge if distance is not explained by a triplet (i,j,k)
  #    M[RE1[i],RE2[i]] <- 1
  #  }
  #}
  
  
  new_weights[RE1, RE2] <- new_weights[RE1, RE2]/10
  
  
  return(new_weights)
}

find_edges <- function(weights){
  
  #thresholds <- find_all_thresholds(weights)
  
  new_weights <- matrix(0, nrow = nrow(weights), ncol = ncol(weights))
  
  for(j in 1:ncol(weights)){
    q <- quantile(weights[,j], prob=0.)
    for(i in 1:nrow(weights)){
      if(weights[i, j] > q){#thresholds[i]){
        new_weights[i, j] <- weights[i,j] #+ 1 
      }else{
        new_weights[i, j] <- 0#weights[i,j]**2/q
      }
    }
    
    new_weights[,j] <- new_weights[,j]/sum(new_weights[,j])
  }
  
  row.names(new_weights) <- row.names(weights)
  colnames(new_weights) <- colnames(weights)
  
  return(new_weights)
}

#Shannon entropy
Shannon_h <- function(x){
  x <- x/sum(x, na.rm = T)
  return(-sum(log(x**x), na.rm = T))
}

#Tsallis entropy 
Tsallis_h <- function(x, q){
  if(q == 1){
    return(Shannon_h(x))
  }else{
    x <- x/sum(x, na.rm=T)
    return(1/(q-1)*(1 - sum(x**q, na.rm=T)))
  }
  
}

#Function that finds a threshold based on Tsallis entropy 
#returns a single number, ie. the threshold

find_one_threshold <- function(x){
  x <- sort(x)
  l <- length(x)
  h <- c()
  for(j in 1:l){
    d <- x[1:j]
    h <- append(h, Tsallis_h(d, 1)) 
  }
  
  return(x[which(h == max(h, na.rm = TRUE))])
}

#Function that returns all the thresholds for a given weight matrix, based on Tsallis entropy 
find_all_thresholds <- function(x){
  thresholds <- list()
  for(i in 1:(ncol(x))){
    thresholds <- append(thresholds, find_one_threshold(x[,i]))
  }
  return(as.vector(thresholds))
}

#Function to produce a final network from a weight matrix. 
#The final network is still weighted one can remove the weights by adding weights[weights>0] = 1

final_network<- function(weights){
  
  thr <- find_all_thresholds(weights)
  
  for(i in 1:(ncol(weights))){
    for(j in 1:nrow(weights)){
      if(weights[j,i] < thr[i]){
        weights[j,i] <- 0 
      }
    }
  }
  
  return(weights)
}


odiag <- function(x) x[col(x) != row(x)]


#Function that calculates correlations
my_corr <- function(x, y){
  #ssum <- sum((x -mean(x))*(y -mean(y) ))/sqrt(sum( (x -mean(x))**2) *sum( (y -mean(y))**2))
  ssum <- sum((x -mean(x))*(y-mean(y) ))/sqrt(sum( (x -mean(x))**2) *sum( (y -mean(y))**2))
  
  return(ssum)
}

#Function to include the protein step
gamma_mat <- function(TS.data, time.points, alphas){
  
  num.exp <- length(TS.data)
  gene.names <- rownames(TS.data[[1]])
  num.genes <- length(gene.names)
  
  Gammas <- matrix(0, nrow = num.genes, ncol = num.genes)
  correlations <- matrix(0, nrow = num.genes, ncol = num.genes)
  
  for(k in 1:num.exp){
    TS.data[[k]] <- t(TS.data[[k]])
  }
  
  
  for(k in 1:num.exp){
    
    current.timeseries <- (TS.data[[k]])
    current.time.points <- time.points[[k]]
    num.points <- length(current.time.points)
    time.diff.current <- current.time.points[(2):num.points] - current.time.points[1:(num.points-1)]
    
    current.timeseries.output <- (current.timeseries[(2):num.points,] - current.timeseries[1:(num.points-1),]) / replicate(num.genes,time.diff.current) 
    + t(replicate((num.points-1), alphas)) * (current.timeseries[(2):num.points,] + current.timeseries[1:(num.points-1),]) /2
    
  
    
    t <- length(time.points[[k]])
    
    print(k)
      
    #loop over targets
    for(j in 1:num.genes){
      
      GeneName <- gene.names[j]
      
      #loop over drivers
      for(i in 1:num.genes){
        
        
        f <- function(g){
          return( abs( cor( as.vector(current.timeseries.output[, j]) ,as.vector(current.timeseries[2:(num.points),i]) + g*(as.vector(current.timeseries[1:(num.points-1),i])) )     ))
        }
        
        res <- optimize(f, c(-10,10), maximum = TRUE)
        
        Gammas[j, i] <- Gammas[j,i] + res$maximum
        correlations[j,i] <- correlations[j, i] +res$objective
      }
    }
    
    
  }
  
  #The matrices have the usual convention where rows are targets and columns are drivers
  
  Gammas <- Gammas/10
  correlations <- correlations/10
  diag(correlations) <- NA
  results <- list("gammas" = Gammas, "correlations" = correlations)
  
  return(results)
  
}

#Function that calculates time-lagged correlations and effective decay rates
#can probably be improved 
different_time_corr <- function(TS.data, time.points){
 
  num.exp <- length(TS.data)
  gene.names <- rownames(TS.data[[1]])
  num.genes <- length(gene.names)
  
  #CAVEAT the correlations have to be normalized well using the variances
  C <- matrix(0, nrow = num.genes, ncol = num.genes)      #different time correlation cor(x_i(t), x_j(t+1))
  D <- matrix(0, nrow = num.genes, ncol = num.genes)      #correlation of variations cor(x_i(t), x_j(t+1)-x_j(t))
  
  A <- matrix(0, nrow = num.genes, ncol = num.exp)
  
  for(k in 1:num.exp){
    
    current.timeseries <- (TS.data[[k]])
    current.time.points <- time.points[[k]]
    num.points <- length(current.time.points)
    time.diff.current <- current.time.points[(2):num.points] - current.time.points[1:(num.points-1)]
  
    t <- length(time.points[[k]])
    
    print(paste("Time series:  ",k,"/",num.exp,sep=""))
    
    #loop over targets
    for(j in 1:num.genes){
      #print(paste("Gene:  ",j,"/",num.genes,sep=""))

      cy <- c()
      #loop over drivers
      for(i in 1:num.genes){
        
       C[i, j] <- C[i, j] + my_corr(current.timeseries[j,1:(t-1)], current.timeseries[i,2:t])/num.exp
       
       dif <- (current.timeseries[j,2:t] - current.timeseries[j,1:(t-1)] )
       D[i, j] <- D[i, j]+ my_corr(current.timeseries[i,2:t], dif)/num.exp
       
       cy <- c(cy, my_corr(current.timeseries[j,1:(t-1)], current.timeseries[i,1:(t-1)])*sqrt(var(current.timeseries[j,1:(t-1)])* var(current.timeseries[i,1:(t-1)]  ))   )
      }
      
      N <- sqrt(var(current.timeseries[j,1:(t-1)]) * var(current.timeseries[j,2:t]))
      
      A[j, k] <- 1/50*(1-N*C[j, j]/var(current.timeseries[j,1:(t-1)])) 
  
        }
  }
  
  alphas <- rowSums(A)/num.exp

  results <- list("C" = C, "D" = D, "alphas" = alphas)
  
  return(results)
  
}

#Function that calculates TS averages 
TS_averages <- function(TS.data, time.points){
  
  num.exp <- length(TS.data)
  gene.names <- rownames(TS.data[[1]])
  num.genes <- length(gene.names)
  
  avg <- matrix(0, nrow = num.genes, ncol = 1)
  full_avg <- matrix(0,  nrow = num.genes, ncol = 1)
  VAR <- matrix(0,  nrow = num.genes, ncol = 1)
  
  for(k in 1:num.exp){
    
    current.timeseries <- (TS.data[[k]])
    current.time.points <- time.points[[k]]
    num.points <- length(current.time.points)
    time.diff.current <- current.time.points[(2):num.points] - current.time.points[1:(num.points-1)]
    
    t <- length(time.points[[k]])
    
    print(k)
    
    #loop over targets
    for(j in 1:num.genes){
      
    full_avg[j] <- full_avg[j] + mean(current.timeseries[j,])/num.exp
    avg[j] <- avg[j] + mean(current.timeseries[j,1:num.points-1])/num.exp
    VAR[j] <- VAR[j] + var(current.timeseries[j,])/num.exp
      }
  
  }
  
  results <- list("avg" = avg,"full_avg" = full_avg, "var" = VAR)
  return(results)
}

#Function that normalizes the priors by row
rowNormalize <- function(M){
  
  for(i in 1:length(M[1,])){ #rows
    Z <- sum(M[i,])
    for(j in 1:length(M[,1])){ #cols
      M[i,j] <- M[i,j]/Z
    }
  }
  
  return(M)
}


average_list_importances=function(list_pred.random) {
  # Takes a list of importance predictions and returns the average importances
  # One has to consider that an absent importance is equal to 0
  num_reps <- length(list_pred.random)
  # get all possible interactions
  all_gene_names <- c()
  for(irep in 1:num_reps) {
    all_gene_names <- c(all_gene_names,unique(list_pred.random[[irep]][,1]))
    all_gene_names <- c(all_gene_names,unique(list_pred.random[[irep]][,2]))
  } 
  all_gene_names <- sort(unique(all_gene_names))
  num_genes <- length(all_gene_names)
  num_interactions <- num_genes*(num_genes-1)
  all_interactions <- c()  # of size num_genes*(num_genes-1)
  all_from <- c()
  all_to   <- c()
  for(ifrom in 1:num_genes) {
    all_from <- c(all_from,rep(all_gene_names[ifrom],num_genes-1))
    all_to   <- c(all_to,all_gene_names[-ifrom])
    all_interactions <- c(all_interactions,paste(all_gene_names[ifrom],"_",all_gene_names[-ifrom],sep=""))
  }
  all_importances <- data.frame(from=all_from,to=all_to)  # we'll add the columns for each rep
  row.names(all_importances) <- all_interactions
  # have each line name be that of the interaction
  for(irep in 1:num_reps) {
    gene_from <- list_pred.random[[irep]][,1]
    gene_to   <- list_pred.random[[irep]][,2]
    these.row.names <- paste(gene_from,"_",gene_to,sep="")
    row.names(list_pred.random[[irep]]) <- these.row.names
    this_col <- rep(0,num_interactions)
    all_importances <- cbind(all_importances,this_col)
    all_importances[these.row.names,irep+2] <- list_pred.random[[irep]][,3]
  }
  mean_importances <- rowMeans(all_importances[,-c(1,2)])
  all_importances <- all_importances[,1:3]
  all_importances[,3] <- mean_importances
  colnames(all_importances)[3] <- "importance"
  all_importances <- all_importances[order(all_importances[,3],decreasing=TRUE),]
  return(all_importances)
}

true_strengths.bs <- function(prediction.bs){
  # true_strengths.bs: retrieve the true value of strength for each interaction
  # in a boot.strength prediction
  #
  # Parameters(required): 
  #   -- prediction.bs: a boot.strength prediction that is a dataframe of 4 columns
  #      "from", "to", "strength" and "direction" corresponding to the interactions
  #      inferred by the method, the respective confidence levels of their presence
  #      and their direction. The problem here is that the strength of each interaction
  #      corresponds to the frequency of appearance of the interaction in both senses
  #      in the bootstraps. The real frequency of each interaction in its given sense
  #      corresponds actually to its level of confidence for its direction multiplied by
  #      its strength.
  #
  # Returns:
  #   a prediction as a dataframe of 3 columns: "from", "to" and "strength" where the
  #   interactions are ordered in the descending order regarding their strength
  #
  n <- length(prediction.bs$from)
  for (i in 1:n){
    prediction.bs$strength[i] <- prediction.bs$strength[i]*prediction.bs$direction[i]
  }
  final.pred <- prediction.bs[,1:3]# we don't keep the "direction" column
  # we order the arcs in descending order regarding their strength
  final.pred <- final.pred[order(final.pred$strength,decreasing=TRUE),]
  return(final.pred)
}

if(!"fs" %in% installed.packages()){
  install.packages("fs")
}
library(fs)
if(!("stringr" %in% installed.packages())) {
  install.packages("stringr")
}
library(stringr)

numeration_of_list_files <- function(list_files){
  # retrieves the numbers (graph number) from the files in the
  # given directory (e.g., 1  10 100  11  12  13...)
  #
  # Parameters(required): 
  #   list_files = output of list_files: the names of the files in the directory
  #
  # Returns:
  #   a vector containing the numbers of those files (e.g., 1  10 100  11  12  13...)
  #
  l <- length(list_files)
  num <- rep(0,l)
  for (i in 1:l){
    split1 <- str_split(list_files[i],pattern="\\.")
    split2 <- str_split(split1[[1]][1],pattern="_")
    num[i] <- split2[[1]][2]
  }
  num <- as.numeric(num)
  return(num)
}




sign_inference <- function(TS.data, time.points){
  
  num.exp <- length(TS.data)
  gene.names <- rownames(TS.data[[1]])
  num.genes <- length(gene.names)
  
  #CAVEAT the correlations have to be normalized well using the variances
  C <- matrix(0, nrow = num.genes, ncol = num.genes)      #different time correlation cor(x_i(t), x_j(t+1))
  D <- matrix(0, nrow = num.genes, ncol = num.genes)      #correlation of variations cor(x_i(t), x_j(t+1)-x_j(t)
  S <- matrix(0, nrow = num.genes, ncol = num.genes) 
  
  for(k in 1:num.exp){
    
    current.timeseries <- (TS.data[[k]])
    current.time.points <- time.points[[k]]
    num.points <- length(current.time.points)
    time.diff.current <- current.time.points[(2):num.points] - current.time.points[1:(num.points-1)]
    
    t <- length(time.points[[k]])
    
    print(paste("Time series:  ",k,"/",num.exp,sep=""))
    
    #loop over targets
    for(j in 1:num.genes){
      
      N <- sqrt(var(current.timeseries[j,1:(t-1)]) * var(current.timeseries[j,2:t]))
      a <- 1/50*(1-N*C[j, j]/var(current.timeseries[j,1:(t-1)])) 
      
      #loop over drivers
      for(i in 1:num.genes){
        
        C[i, j] <- C[i, j] + my_corr(current.timeseries[j,1:(t-1)], current.timeseries[i,2:t])/num.exp
        
        dif <- (current.timeseries[j,2:t] - current.timeseries[j,1:(t-1)] )
        D[i, j] <- D[i, j]+ my_corr(current.timeseries[i,2:t], dif)/num.exp
      
        S[i, j] <- S[i, j] + D[i, j]*sd(dif)/sd(current.timeseries[i,])/50 +
          C[i, j]*sd(current.timeseries[j,])/sd(current.timeseries[i,])*a
        
        }
      
      
    }
  }
  
  results <- list("C" = C, "D" = D, "S" = S)
  
  return(results)
  
}







