#' Hard clustering using k-means
#'
#' This function expects output from custom minimap test dataset that contains original locations of mapped reads in the genome.
#'
#' @param tab.l A \code{list} of PB alignmetns separated per cell.
#' @param clusters Desired number of clusters.
#' @return A \code{list} of estimated theta values for every cluster and cell.
#' @author David Porubsky, Maryam Ghareghani
#' @export

hardClust <- function(tab.l=NULL, clusters=NULL) {

  message("Hard clustering")
  ptm <- startTimedMessage("    Counting directional reads ...") 
  ratios.l <- list()
  counts.l <- list()
  for (j in 1:length(tab.l)) {
    lib.name <- names(tab.l[j])
    #message("\tWorking on ",lib.name)
    lib.aligns <- tab.l[[j]]
    aligns.per.read <- split(lib.aligns$strand, lib.aligns$PBreadNames)
  
    counts <- t(sapply(aligns.per.read, function(x) table(x))) #count directional reads per PB read
    ratios <- (counts[,2]-counts[,1])/(counts[,2]+counts[,1]) #calculate ratio of WW reads
    ratios[is.nan(ratios)] <- 0
    ratios.l[[j]] <- ratios
    counts.l[[j]] <- counts
  }
  stopTimedMessage(ptm)

  ptm <- startTimedMessage("    Kmeans clustering for ",clusters," clusters ...") 
  ratios.m <- do.call(cbind, ratios.l)
  km <- suppressWarnings( kmeans(ratios.m, centers = clusters, nstart = 100) )
  ord <- km$cluster
  ratios.m.ord <- ratios.m[order(ord),]
  stopTimedMessage(ptm)

  ptm <- startTimedMessage("    Estimate theta values ...") 
  theta.estim <- list()
  for (j in 1:length(counts.l)) {
    minus.c <- split(counts.l[[j]][,1], ord)
    plus.c <- split(counts.l[[j]][,2], ord)
  
    clust.prob <- mapply(function(X,Y) { countProb(X,Y) }, X=minus.c, Y=plus.c)
    clust.prob.norm <- lapply(clust.prob, function(x) colSums(log(x)))
    estimates <- sapply(clust.prob.norm, which.max)
    
    #Assign cell type probs based on the majority type in each cluster
    probs <- list()
    for (i in 1:length(clust.prob)) {
      estim <- estimates[i]
      if (estim == 1) {
        theta <- c(0.9, 0.05, 0.05)
      } else if (estim == 2) {
        theta <- c(0.05, 0.9, 0.05)
      } else {
        theta <- c(0.05, 0.05, 0.9)
      }
      probs[[i]] <- theta
    }
    probs <- do.call(rbind, probs)
    theta.estim[[j]] <- probs
  }
  stopTimedMessage(ptm)
  
  return(list(theta.estim=theta.estim, clust.id=ord, raw.counts=counts.l))
}  