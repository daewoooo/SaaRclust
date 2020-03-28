#' Get pairs of clusters coming from the same chromosome that differs in directionality.
#'
#' This function solves problem of finding the best clustering partners using linear programing function.
#' 
#' @param theta.param A \code{list} of estimated cell types for each cluster and each cell.
#' @return A \code{matrix} of pairs of clusters IDs that belong to the same chromosome.
#' @importFrom lpSolve lp.assign
#' @author David Porubsky
#' 
findClusterPartners <- function(theta.param=NULL) {
  
  ## If there is an uneven number of clusters remove the one with the most WC states
  num.clusters <- nrow(theta.param[[1]])
  if (num.clusters %% 2 != 0) {
    #Find cluster with WC state in majority of cells
    theta.sums <- Reduce("+", theta.param)
    remove.clust <- which.max(theta.sums[,3])
    message("    Removed cluster ", remove.clust, " to ensure even number of clusters!!!")
    theta.param <- lapply(theta.param, function(x) x[-remove.clust,])
  }
  
  ## get only wc thetas
  theta.param.wc <- lapply(theta.param, function(x) x[,3])
  ## cbind wc thetas for all single cells
  all.theta.param.wc <- do.call(cbind, theta.param.wc)
  ## compute the pairwise distance of all clusters wc thetas
  d <- as.matrix(dist(all.theta.param.wc))
  ## convert distance to a similarity measure
  d <- max(d) - d
  ## set diagonal values to zero
  diag(d) <- 0
  ## Find pairs of clusters with the highest similarity
  max.partners <- lpSolve::lp.assign(d, "max")
  max.partners.m <- max.partners$solution #matrix with pairs of clusters with maximal similarity
  ## Extract indices of pair of clusters
  max.partners.idx <- which(max.partners.m > 0, arr.ind = TRUE)
  max.partners.idx <- max.partners.idx[max.partners.idx[,1] < max.partners.idx[,2],] #remove duplicate cluster partners
  colnames(max.partners.idx) <- c('Cluster1', 'Cluster2')
  return(max.partners.idx)
}