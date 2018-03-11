
#' Get pairs of clusters coming from the same chromosome that differs in directionality.
#'
#' This function solves the Maximum matching problem for all possible pairs of clusters and reports pairs with highest similarity.
#'
#' @param theta.param A \code{list} of estimated cell types for each cluster and each cell.
#' @return A \code{vector} of pairs of clusters IDs that belong to the same chromosome.
#' @importFrom maxmatching maxmatching
#' @importFrom igraph graph E
#' @author David Porubsky, Maryam Ghareghani
#' @export

findClusterPartners <- function(theta.param=NULL) {
  
  euc.dist.v <- function(v) sqrt(sum((v[1] - v[2]) ^ 2)) #calculate euclidean distance for pair of datapoints
  
  pairwise.dist <- list()
  for (i in 1:length(theta.param)) {
    cell.theta <- theta.param[[i]]
    pairs <- t(combn(nrow(cell.theta), 2))
    pairs.m <- cbind( cell.theta[pairs[,1],1], cell.theta[pairs[,2],2] )
    dist <- apply(pairs.m, 1, euc.dist.v)
    pairwise.dist[[i]] <- dist
  }
  pairwise.dist.m <- do.call(cbind, pairwise.dist)
  pairwise.simil.m <- max(pairwise.dist.m)-pairwise.dist.m
  simil.sum <- rowSums(pairwise.simil.m)
  G <- igraph::graph(c(rbind(pairs[,1], pairs[,2])), directed = FALSE)
  igraph::E(G)$weight <- simil.sum
  max.match <- maxmatching::maxmatching(G, weighted = TRUE)
  
  return(max.match$matching)
}


#' Get pairs of clusters coming from the same chromosome
#'
#' This function finds cluster coming from the same chromosome and having the same directionality
#'
#' @param theta.param A \code{list} of estimated cell types for each cluster and each cell.
#' @return A \code{vector} of pairs of clusters IDs that belong to the same chromosome.
#' @importFrom igraph graph max_cliques
#' @author David Porubsky
#' @export

findSplitedClusters <- function(theta.param=NULL) {
  
  euc.dist.v <- function(v) sqrt(sum((v[1] - v[2]) ^ 2))
  
  pairwise.dist <- list()
  for (i in 1:length(theta.param)) {
    cell.theta <- theta.param[[i]]
    pairs <- t(combn(nrow(cell.theta), 2))
    pairs.w <- cbind( cell.theta[pairs[,1],1], cell.theta[pairs[,2],1] )
    dist <- apply(pairs.w, 1, euc.dist.v)
    pairwise.dist[[i]] <- dist
  }
  pairwise.dist.m <- do.call(cbind, pairwise.dist)
  pairwise.simil.m <- max(pairwise.dist.m)-pairwise.dist.m
  simil.sum <- rowSums(pairwise.simil.m)
  zscores <- (simil.sum - mean(simil.sum)) / sd(simil.sum)
  idx <- zscores > 3.291 #99.9 confidence interval
  G <- igraph::graph(c(rbind(pairs[idx,1], pairs[idx,2])), directed = FALSE)
  sub.G <- igraph::max_cliques(G, min = 2)
  
  return(sub.G)
}
