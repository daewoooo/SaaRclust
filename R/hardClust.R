#' Hard clustering using k-means
#'
#' This function expects output from custom minimap test dataset that contains original locations of mapped reads in the genome.
#'
#' @param counts.l A \code{list} of directional read counts per long read/contig per library.
#' @param num.clusters Desired number of clusters or cluster centers (k-means clustering).
#' @param nstart Number of random sets to choose cluster centers (k-means clustering).
#' @param iter.max A maximum number of allowed interations (k-means clustering).
#' @return A \code{list} of estimated theta values for every cluster and cell.
#' @importFrom stats kmeans
#' @author David Porubsky
#' @export
#' @examples
#'## Get an example file
#'exampleFile <- system.file("extdata", "rawCounts_5e+06bp_dynamic.RData", package = "SaaRclust")
#'## Load BAM count table
#'counts.l <- get(load(exampleFile))
#'## Get hard clustering results
#'hardClust.ord <- hardClust(counts.l, num.clusters=100, nstart = 100)
#'
hardClust <- function(counts.l=NULL, num.clusters=NULL, nstart=10, iter.max=10) {

  message("Hard clustering")
  ptm <- startTimedMessage("    Kmeans clustering for ",num.clusters," clusters") 
  
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
  #ratios.m[ratios.m < 0.5] <- -1
  #ratios.m[ratios.m > 0.5] <- 1
  ratios.m[] <- as.numeric(findInterval(ratios.m, c(-0.5, 0.5)))
  
  #hard clustering using kmeans
  km <- suppressWarnings( stats::kmeans(ratios.m, centers = num.clusters, nstart = nstart, iter.max = iter.max) )
  ord <- km$cluster
  #ratios.m.ord <- ratios.m[order(ord),]
  stopTimedMessage(ptm)
  return(ord)
}  
  
#' Estimate theta values based on hard clustering
#'
#' This function takes results of hard clustering and estimates majority cell types for each Strand-seq library
#'
#' @param hard.clust A \code{integer} of cluster assignments for each long read or contig. 
#' @inheritParams hardClust
#' @inheritParams countProb
#' @return A \code{list} of estimated theta values for every cluster and cell.
#' @author David Porubsky
#' @export
#' @examples
#'## Get an example file
#'exampleFile <- system.file("extdata", "rawCounts_5e+06bp_dynamic.RData", package = "SaaRclust")
#'## Load BAM count table
#'counts.l <- get(load(exampleFile))
#'## Get hard clustering results
#'hardClust.ord <- hardClust(counts.l, num.clusters=100, nstart = 100)
#'## Estimate theta parameter
#'theta.param <- estimateTheta(counts.l, hard.clust=hardClust.ord, alpha=0.1)
#'
estimateTheta <- function(counts.l=NULL, hard.clust=NULL, alpha=0.1) {
  ptm <- startTimedMessage("Estimate theta values") 
  theta.estim <- list()
  for (j in 1:length(counts.l)) {
    minus.c <- split(counts.l[[j]][,1], hard.clust)
    plus.c <- split(counts.l[[j]][,2], hard.clust)
    
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
  
  return(theta.estim)
}
 
                             
#' Hierarchical clustering for merging the kmeans clusters.
#'
#' This function takes as input the kmeans hard clustering output and the initialized thetas and merges the kmeans clusters based on thetas
#'
#' @param theta.param A \code{list} of estimated theta values for each cluster and cell.
#' @param k Desired number of clusters after merging.
#' @inheritParams estimateTheta
#' @return A new hard clustering with the correct number of clusters
#' @importFrom stats dist hclust cutree
#' @author Maryam Ghareghani, David Porubsky
#' 
mergeClusters <- function(hard.clust, theta.param, k=46)
{
  ptm <- startTimedMessage("Merging clusters")  
  theta.all <- do.call(cbind, theta.param)
  hc <- stats::hclust(stats::dist(theta.all))
  hc.clust <- stats::cutree(hc, k=k)
  
  stopTimedMessage(ptm)
  return(sapply(hard.clust, function(i) hc.clust[i]))
}
