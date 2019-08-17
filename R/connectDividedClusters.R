#' Get pairs of clusters likely coming from the same chromosome
#'
#' This function finds cluster likely coming from the same chromosome by looking for clusters that share the same
#' directionality across multiple Strand-seq libraries.
#'
#' @param theta.param A \code{list} of estimated cell types for each cluster and each cell.
#' @param z.limit Connect clusters with z-score equal or above this limit.
#' @param remove.always.WC Set to \code{TRUE} if the cluster with majority of WC states should be removed.
#' @return A \code{matrix} of pairs of clusters IDs that belong to the same chromosome.
#' @author David Porubsky
#' @export
#' 
connectDividedClusters <- function(theta.param=NULL, z.limit=3.29, remove.always.WC=FALSE) {
  
  ## Helper function ##
  ## Calculate euclidean distance for pair of datapoints
  euc.dist.v <- function(v) sqrt(sum((v[1] - v[2]) ^ 2))
  
  ## Remove cluster with the most WC states
  if (remove.always.WC) {
    #Find cluster with WC state in majority of cells
    theta.sums <- Reduce("+", theta.param)
    remove.clust <- which.max(theta.sums[,3])
    message("    Removed cluster ", remove.clust, " to ensure even number of clusters!!!")
    theta.param <- lapply(theta.param, function(x) x[-remove.clust,])
  }
  
  dist.wc <- list()
  dist.ww <- list()
  dist.cc <- list()
  dist.het <- list()
  for (i in 1:length(theta.param)) {
    cell.theta <- theta.param[[i]]
    ## Get all possible cluster pairs
    pairs <- t(combn(nrow(cell.theta), 2))
    ## Calculate similarity for WC pairs
    pairs.wc <- cbind( cell.theta[pairs[,1],3], cell.theta[pairs[,2],3] )
    dist.wc[[i]] <- apply(pairs.wc, 1, euc.dist.v)
    ## Calculate similarity for WW pairs
    pairs.ww <- cbind( cell.theta[pairs[,1],1], cell.theta[pairs[,2],1] )
    dist.ww[[i]] <- apply(pairs.ww, 1, euc.dist.v)
    ## Calculate similarity for CC pairs
    pairs.cc <- cbind( cell.theta[pairs[,1],2], cell.theta[pairs[,2],2] )
    dist.cc[[i]] <- apply(pairs.cc, 1, euc.dist.v)
    ## Check cluster connections for HET inv
    ## HET inversion appears as region which is always WC while the rest of the chromsome is CC or WW (or vice versa)
    WWorCC <- pmax(cell.theta[pairs[,1],1], cell.theta[pairs[,1],2]) # Get prob. for WW or CC state
    pairs.het <- cbind( WWorCC, cell.theta[pairs[,2],3] )
    dist.het[[i]] <- apply(pairs.het, 1, euc.dist.v)
    
  }
  ## Get the most significant connections for WC state
  dist.wc.m <- do.call(cbind, dist.wc)
  simil.wc.m <- max(dist.wc.m) - dist.wc.m
  simil.wc.sum <- rowSums(simil.wc.m)
  zscores <- (simil.wc.sum - mean(simil.wc.sum)) / sd(simil.wc.sum)
  idx <- zscores > z.limit
  vertices.wc <- c(rbind(pairs[idx,1], pairs[idx,2]))
  ## Get the most significant connections for WW state
  dist.ww.m <- do.call(cbind, dist.ww)
  simil.ww.m <- max(dist.ww.m) - dist.ww.m
  simil.ww.sum <- rowSums(simil.ww.m)
  zscores <- (simil.ww.sum - mean(simil.ww.sum)) / sd(simil.ww.sum)
  idx <- zscores > z.limit
  vertices.ww <- c(rbind(pairs[idx,1], pairs[idx,2]))
  ## Get the most significant connections for CC state
  dist.cc.m <- do.call(cbind, dist.cc)
  simil.cc.m <- max(dist.cc.m) - dist.cc.m
  simil.cc.sum <- rowSums(simil.cc.m)
  zscores <- (simil.cc.sum - mean(simil.cc.sum)) / sd(simil.cc.sum)
  idx <- zscores > z.limit
  vertices.cc <- c(rbind(pairs[idx,1], pairs[idx,2]))
  ## Get the most significant connections for HET inversion cases
  dist.het.m <- do.call(cbind, dist.het)
  simil.het.m <- max(dist.het.m) - dist.het.m
  simil.het.sum <- rowSums(simil.het.m)
  zscores <- (simil.het.sum - mean(simil.het.sum)) / sd(simil.het.sum)
  idx <- zscores > z.limit
  vertices.het <- c(rbind(pairs[idx,1], pairs[idx,2]))
  ## Get cluster ID of putative het INVs
  putative.HETs <- unique(vertices.het[duplicated(vertices.het)])
  
  ## Merge all vertices
  #vertices <- c(vertices.ww, vertices.cc)
  vertices <- c(vertices.ww, vertices.cc, vertices.wc, vertices.het)
  ## Find strongly connected clusters
  G <- igraph::graph(vertices, directed = FALSE)
  clusters <- groups(components(G, mode = 'strong'))
  
  return(list(clusters=clusters, putative.HETs=putative.HETs))
}
