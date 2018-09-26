#' Get pairs of clusters coming from the same chromosome that differs in directionality.
#'
#' This function solves problem of finding the best clustering partners using linear programing function.
#' 
#' @param theta.param A \code{list} of estimated cell types for each cluster and each cell.
#' @return A \code{matrix} of pairs of clusters IDs that belong to the same chromosome.
#' @importFrom lpSolve lp.assign
#' @author David Porubsky
#' @export

findClusterPartners <- function(theta.param=NULL) {
    # get only wc thetas
    theta.param.wc <- lapply(theta.param, function(x) x[,3])
    # cbid wc thetas for all single cells
    all.theta.param.wc <- do.call(cbind, theta.param.wc)
    # compute the pairwise distance of all clusters wc thetas
    d <- as.matrix(dist(all.theta.param.wc))
    # convert distance to a similarity measure
    d <- max(d) - d
    # set diagonal values to zero
    diag(d) <- 0
    # Find pairs of clusters with the highest similarity
    max.partners <- lp.assign(d, "max")
    max.partners.m <- max.partners$solution #matrix with pairs of clusters with maximal similarity
    # Extract indices of pair of clusters
    max.partners.idx <- which(max.partners.m > 0, arr.ind = TRUE)
    max.partners.idx <- max.partners.idx[max.partners.idx[,1] < max.partners.idx[,2],] #remove duplicate cluster partners
    colnames(max.partners.idx) <- c('Cluster1', 'Cluster2')
    return(max.partners.idx)
}
    
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

findClusterPartners_maxMatch <- function(theta.param=NULL) {
  
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
    
    

#' Get pairs of clusters coming from the same chromosome with different directionalities
#'
#' This function finds the pairs by computing the pairwise distances between clusters theta_wc parameters, and pais each cluster with the cluster with min dist
#' Note that the first 23 minimum values should form a matching in the complete graph of clusters (this condition is observed to be held in practice)
#' TODO: We should replace this simple function by a proper max matching algorithm or some other solution that properly solves this problem and finds non-overlappinh pairs of clusters
#' 
#' @param theta.param A \code{list} containing estimated cell types per cluster for each cell.
#' @return A \code{vector} of pairs of clusters IDs that belong to the same chromosome.
#' @importFrom maxmatching maxmatching
#' @importFrom igraph graph E
#' @author Maryam Ghareghani
#' @export

findClusterPartners_simple <- function(theta.param=NULL) {
  
  # get only wc thetas
  theta.param.wc <- lapply(theta.param, function(x) x[,3])
  # cbid wc thetas for all single cells
  all.theta.param.wc <- do.call(cbind, theta.param.wc)
  # compute the pairwise distance of all wc thetas (all pairs of clusters)
  d <- as.matrix(dist(all.theta.param.wc))
  # compute the min index in each row (excluding the diogonal element)
  min.cost.clust.pair <- sapply(1:nrow(d), function(i) return(which.min(d[i,-i])))
  min.cost.clust.pair <- data.table(first_clust=1:nrow(d), second_clust=as.numeric(names(min.cost.clust.pair)))
  # put the min rank cluster in the first column
  min.cost.clust.pair[, min_id_clust:=min(first_clust, second_clust), by=1:nrow(min.cost.clust.pair)]
  min.cost.clust.pair[, second_clust:=max(first_clust, second_clust), by=1:nrow(min.cost.clust.pair)]
  min.cost.clust.pair[, `:=`(first_clust=min_id_clust, min_id_clust=NULL)]
  
  # compute the frequencies of each cluster in this pairing
  fr <- table(c(min.cost.clust.pair$first_clust, min.cost.clust.pair$second_clust))
  
  # check whether the pairing is a matching
  assert_that(all(sort(unique(fr))==1:3) & length(which(duplicated(min.cost.clust.pair)))==floor(nrow(min.cost.clust.pair)/2)) %>% invisible()
  
  # keep only the duplicated rows (keep the matching and kick out the garbage cluster)
  min.cost.clust.pair <- min.cost.clust.pair[duplicated(min.cost.clust.pair)]
  
  return(min.cost.clust.pair)
}
