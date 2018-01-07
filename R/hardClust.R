#' Hard clustering using k-means
#'
#' This function expects output from custom minimap test dataset that contains original locations of mapped reads in the genome.
#'
#' @param counts.l A \code{list} of directional read counts per PB read per library.
#' @inheritParams SaaRclust
#' @return A \code{list} of estimated theta values for every cluster and cell.
#' @author David Porubsky
#' @export

hardClust <- function(counts.l=NULL, num.clusters=NULL, alpha=0.1, nstart=10, iter.max=10) {

  message("Hard clustering")
  ptm <- startTimedMessage("    Kmeans clustering for ",num.clusters," clusters ...") 
  
  ratios.l <- list()
  for (j in 1:length(counts.l)) {
    #lib.name <- names(tab.l[j])
    #message("\tWorking on ",lib.name)
    counts <- counts.l[[j]]
   
    ratios <- (counts[,2]-counts[,1])/(counts[,2]+counts[,1]) #calculate ratio of WW reads
    ratios[is.nan(ratios)] <- 0
    ratios.l[[j]] <- ratios
  }
   
  ratios.m <- do.call(cbind, ratios.l)
  ratios.m[ratios.m<0] <- -1
  ratios.m[ratios.m>0] <- 1
  km <- suppressWarnings( kmeans(ratios.m, centers = num.clusters, nstart = nstart, iter.max = iter.max) )
  ord <- km$cluster
  #ratios.m.ord <- ratios.m[order(ord),]
  stopTimedMessage(ptm)

  ptm <- startTimedMessage("    Estimate theta values ...") 
  theta.estim <- list()
  for (j in 1:length(counts.l)) {
    minus.c <- split(counts.l[[j]][,1], ord)
    plus.c <- split(counts.l[[j]][,2], ord)
    
    #minus.counts <- sapply(minus.c, sum)
    #plus.counts <- sapply(plus.c, sum)
    #probs <- countProb2(minusCounts = minus.counts, plusCounts = plus.counts)
  
    clust.prob <- mapply(function(X,Y) { countProb(X,Y) }, X=minus.c, Y=plus.c)
    clust.prob.norm <- lapply(clust.prob, function(x) colSums(log(x)))
    estimates <- sapply(clust.prob.norm, which.max)
    
    #Assign cell type probs based on the majority type in each cluster
    probs <- list()
    for (i in 1:length(clust.prob)) {
      estim <- estimates[i]
      if (estim == 1) {
        theta <- c(1-alpha, alpha/2, alpha/2)
      } else if (estim == 2) {
        theta <- c(alpha/2, 1-alpha, alpha/2)
      } else {
        theta <- c(alpha/2, alpha/2, 1-alpha)
      }
      probs[[i]] <- theta
    }
    probs <- do.call(rbind, probs)
    theta.estim[[j]] <- probs
  }
  stopTimedMessage(ptm)
  
  #return(list(theta.estim=theta.estim, clust.id=ord, raw.counts=counts.l))
  return(list(theta.estim=theta.estim, clust.id=ord))
}
                              
#' Hierarchical clustering for merging the kmeans clusters
#'
#' This function takes as input the kmeans hard clustering output and the initialized thetas and merges the kmeans clusters based on thetas
#'
#' @param theta.l A \code{list} of estimated theta values for each cluster and cell.
#' @param kmeans.clust The kmeans hard clustering.
#' @inheritParams SaaRclust
#' @return A new hard clustering with the correct number of clusters
#' @author Maryam Ghareghani
#' @export


mergeClusters <- function(kmeans.clust, theta.l)
{
  theta.all <- do.call(cbind, theta.l)
  hc <- hclust(dist(theta.all))
  hc.clust <- cutree(hc, k=46)
  
  return(list(clust.id = sapply(kmeans.clust$clust.id, function(i) hc.clust[i])))
}

#' Get the feature vector based on the WmiunsC ratios
#'
#' This function expects output from custom minimap test dataset that contains original locations of mapped reads in the genome.
#'
#' @param counts.l A \code{list} of directional read counts per PB read per library.
#' @inheritParams SaaRclust
#' @return A \code{list} of estimated theta values for every cluster and cell.
#' @author Maryam Ghareghani
#' @export

WminusCratiosFeatures <- function(counts.l=NULL) {
  ratios.l <- list()
  for (j in 1:length(counts.l)) {
    #lib.name <- names(tab.l[j])
    #message("\tWorking on ",lib.name)
    counts <- counts.l[[j]]
    
    ratios <- (counts[,2]-counts[,1])/(counts[,2]+counts[,1]) #calculate ratio of WW reads
    ratios[is.nan(ratios)] <- 0
    ratios.l[[j]] <- ratios
  }
  
  ratios.m <- do.call(cbind, ratios.l)
  #ratios.m[ratios.m<0] <- -1
  #ratios.m[ratios.m>0] <- 1
  
  ratios.m
}

#' Get the feature vector based on the max{W,C} ratios
#'
#' This function expects output from custom minimap test dataset that contains original locations of mapped reads in the genome.
#'
#' @param counts.l A \code{list} of directional read counts per PB read per library.
#' @inheritParams SaaRclust
#' @return A \code{list} of estimated theta values for every cluster and cell.
#' @author Maryam Ghareghani
#' @export

maxWandCratiosFeatures <- function(counts.l=NULL) {
  ratios.l <- list()
  for (j in 1:length(counts.l)) {
    #lib.name <- names(tab.l[j])
    #message("\tWorking on ",lib.name)
    counts <- counts.l[[j]]
    
    ratios <- counts[,2]/(counts[,2]+counts[,1])  #calculate ratio of WW reads
    ratios <- sapply(ratios, function(x) max(x, 1-x))
    ratios[is.nan(ratios)] <- 0
    ratios.l[[j]] <- ratios
  }
  
  ratios.m <- do.call(cbind, ratios.l)
  #ratios.m[ratios.m<0] <- -1
  #ratios.m[ratios.m>0] <- 1
  
  ratios.m
}

#' Get the feature vector based on both W and C ratios
#'
#' This function expects output from custom minimap test dataset that contains original locations of mapped reads in the genome.
#'
#' @param counts.l A \code{list} of directional read counts per PB read per library.
#' @inheritParams SaaRclust
#' @return A \code{list} of estimated theta values for every cluster and cell.
#' @author Maryam Ghareghani
#' @export

WandCratiosFeatures <- function(counts.l=NULL) {
  ratios.l <- list()
  for (j in 1:length(counts.l)) {
    #lib.name <- names(tab.l[j])
    #message("\tWorking on ",lib.name)
    counts <- counts.l[[j]]
    
    ratios <- counts/(counts[,2]+counts[,1]) #calculate ratio of WW reads
    ratios[is.nan(ratios)] <- 0
    ratios.l[[j]] <- ratios
  }
  
  ratios.m <- do.call(cbind, ratios.l)
  #ratios.m[ratios.m<0] <- -1
  #ratios.m[ratios.m>0] <- 1
  
  ratios.m
}
