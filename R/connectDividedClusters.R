#' Merge divided clusters likely coming from the same chromosome
#'
#' This function finds cluster likely coming from the same chromosome by looking for clusters that share the same
#' directionality across multiple Strand-seq libraries.
#'
#' @param theta.param A \code{list} of estimated cell types for each cluster and each cell.
#' @param z.limit Connect clusters with z-score equal or above this limit.
#' @param remove.always.WC Set to \code{TRUE} if the cluster with majority of WC states should be removed.
#' @param desired.num.clusters Desired number of clusters after cluster merging step.
#' @return A \code{matrix} of pairs of clusters IDs that belong to the same chromosome.
#' @importFrom igraph graph groups components
#' @author David Porubsky
#' @export
#'
connectDividedClusters <- function(theta.param=NULL, z.limit=3.29, remove.always.WC=FALSE, desired.num.clusters=NULL) {
  
  ## Helper function ##
  ## Calculate euclidean distance for pair of datapoints
  euc.dist.v <- function(v) sqrt(sum((v[1] - v[2]) ^ 2))
  
  ptm <- startTimedMessage("Detecting divided clusters")
  
  ## Find clusters with WC state in majority of cells
  theta.sums <- Reduce("+", theta.param)
  theta.zscore <- (theta.sums[,3] - mean(theta.sums[,3])) / sd(theta.sums[,3])
  wc.clust.idx <- which(theta.zscore > 2.5)
  ## Remove cluster with the most WC states
  if (remove.always.WC) {
    theta.param <- lapply(theta.param, function(x) x[-wc.clust.idx,])
  }
  
  ## Get all possible cluster pairs
  pairs <- t(combn(nrow(theta.param[[1]]), 2))
  
  dist.wc <- list()
  dist.ww <- list()
  dist.cc <- list()
  dist.het <- list()
  for (i in 1:length(theta.param)) {
    cell.theta <- theta.param[[i]]
    
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
  vertices.wc <- cbind(pairs, zscores)
  vertices.wc <- vertices.wc[order(vertices.wc[,3], decreasing = TRUE),]
  ## Get the most significant connections for WW state
  dist.ww.m <- do.call(cbind, dist.ww)
  simil.ww.m <- max(dist.ww.m) - dist.ww.m
  simil.ww.sum <- rowSums(simil.ww.m)
  zscores <- (simil.ww.sum - mean(simil.ww.sum)) / sd(simil.ww.sum)
  vertices.ww <- cbind(pairs, zscores)
  vertices.ww <- vertices.ww[order(vertices.ww[,3], decreasing = TRUE),]
  ## Get the most significant connections for CC state
  dist.cc.m <- do.call(cbind, dist.cc)
  simil.cc.m <- max(dist.cc.m) - dist.cc.m
  simil.cc.sum <- rowSums(simil.cc.m)
  zscores <- (simil.cc.sum - mean(simil.cc.sum)) / sd(simil.cc.sum)
  vertices.cc <- cbind(pairs, zscores)
  vertices.cc <- vertices.cc[order(vertices.cc[,3], decreasing = TRUE),]
  ## Get the most significant connections for HET inversion cases
  dist.het.m <- do.call(cbind, dist.het)
  simil.het.m <- max(dist.het.m) - dist.het.m
  simil.het.sum <- rowSums(simil.het.m)
  zscores <- (simil.het.sum - mean(simil.het.sum)) / sd(simil.het.sum)
  vertices.het <- cbind(pairs, zscores)
  vertices.het <- vertices.het[order(vertices.het[,3], decreasing = TRUE),]
  ## Get cluster ID of putative het INVs
  putative.HETs <- c(vertices.het[,1][vertices.het[,3] > z.limit], vertices.het[,2][vertices.het[,3] > z.limit])
  putative.HETs <- unique(putative.HETs[duplicated(putative.HETs)])
  
  ## Merge z-scores across all pairwise comparisons and order them by decreasing z-score
  vertices <- rbind(vertices.wc, vertices.ww, vertices.cc, vertices.het)  
  vertices <- vertices[order(vertices[,3], decreasing = TRUE),]
  
  ## There is a user defined desired.num.clusters keep decresing z.limit to reach user defined
  ## number of clusters
  if (!is.null(desired.num.clusters)) {
    cl.num <- nrow(theta.param[[1]])
    clusters <- NULL
    clusters.prev <- NULL
    while (desired.num.clusters < cl.num) {
      vertices.sub <- vertices[vertices[,3] >= z.limit,]
      vertices.sub <- c(rbind(vertices.sub[,1], vertices.sub[,2]))
      ## Keep previous iteration 
      if (!is.null(clusters)) {
        clusters.prev <- clusters
      }
      ## Find strongly connected clusters
      G <- igraph::graph(vertices.sub, directed = FALSE)
      clusters <- igraph::groups(igraph::components(G, mode = 'strong'))
      cl.num <- length(clusters)
      z.limit <- z.limit - 0.1
    }
    ## Report previous iteration in case current number of clusters is smaller than 'desired.num.clusters'
    if (length(clusters) < desired.num.clusters & !is.null(clusters.prev)) {
      clusters <- clusters.prev
    }
  } else {
    vertices.sub <- vertices[vertices[,3] >= z.limit,]
    vertices.sub <- c(rbind(vertices.sub[,1], vertices.sub[,2]))
    ## Find strongly connected clusters
    G <- igraph::graph(vertices.sub, directed = FALSE)
    clusters <- igraph::groups(igraph::components(G, mode = 'strong'))
  }
  
  stopTimedMessage(ptm)
  return(list(clusters=clusters, putative.HETs=putative.HETs))
}

## Older version of cluster merging function ##
#' #' Get pairs of clusters likely coming from the same chromosome
#' #'
#' #' This function finds cluster likely coming from the same chromosome by looking for clusters that share the same
#' #' directionality across multiple Strand-seq libraries.
#' #'
#' #' @param theta.param A \code{list} of estimated cell types for each cluster and each cell.
#' #' @param z.limit Connect clusters with z-score equal or above this limit.
#' #' @param remove.always.WC Set to \code{TRUE} if the cluster with majority of WC states should be removed.
#' #' @return A \code{matrix} of pairs of clusters IDs that belong to the same chromosome.
#' #' @importFrom igraph graph groups components
#' #' @author David Porubsky
#' #' @export
#' #'
#' connectDividedClusters <- function(theta.param=NULL, z.limit=3.29, remove.always.WC=FALSE) {
#' 
#'   ## Helper function ##
#'   ## Calculate euclidean distance for pair of datapoints
#'   euc.dist.v <- function(v) sqrt(sum((v[1] - v[2]) ^ 2))
#' 
#'   ptm <- startTimedMessage("Detecting divided clusters")
#' 
#'   ## Find cluster with WC state in majority of cells
#'   theta.sums <- Reduce("+", theta.param)
#'   max.WC <- which.max(theta.sums[,3])
#'   ## Remove cluster with the most WC states
#'   if (remove.always.WC) {
#'     message("\n    Removed cluster ", max.WC, " with the most WC states!!!")
#'     theta.param <- lapply(theta.param, function(x) x[-max.WC,])
#'   }
#' 
#'   dist.wc <- list()
#'   dist.ww <- list()
#'   dist.cc <- list()
#'   dist.het <- list()
#'   for (i in 1:length(theta.param)) {
#'     cell.theta <- theta.param[[i]]
#'     ## Get all possible cluster pairs
#'     pairs <- t(combn(nrow(cell.theta), 2))
#'     ## Calculate similarity for WC pairs
#'     pairs.wc <- cbind( cell.theta[pairs[,1],3], cell.theta[pairs[,2],3] )
#'     dist.wc[[i]] <- apply(pairs.wc, 1, euc.dist.v)
#'     ## Calculate similarity for WW pairs
#'     pairs.ww <- cbind( cell.theta[pairs[,1],1], cell.theta[pairs[,2],1] )
#'     dist.ww[[i]] <- apply(pairs.ww, 1, euc.dist.v)
#'     ## Calculate similarity for CC pairs
#'     pairs.cc <- cbind( cell.theta[pairs[,1],2], cell.theta[pairs[,2],2] )
#'     dist.cc[[i]] <- apply(pairs.cc, 1, euc.dist.v)
#'     ## Check cluster connections for HET inv
#'     ## HET inversion appears as region which is always WC while the rest of the chromsome is CC or WW (or vice versa)
#'     WWorCC <- pmax(cell.theta[pairs[,1],1], cell.theta[pairs[,1],2]) # Get prob. for WW or CC state
#'     pairs.het <- cbind( WWorCC, cell.theta[pairs[,2],3] )
#'     dist.het[[i]] <- apply(pairs.het, 1, euc.dist.v)
#' 
#'   }
#'   ## Get the most significant connections for WC state
#'   dist.wc.m <- do.call(cbind, dist.wc)
#'   simil.wc.m <- max(dist.wc.m) - dist.wc.m
#'   simil.wc.sum <- rowSums(simil.wc.m)
#'   zscores <- (simil.wc.sum - mean(simil.wc.sum)) / sd(simil.wc.sum)
#'   idx <- zscores > z.limit
#'   vertices.wc <- c(rbind(pairs[idx,1], pairs[idx,2]))
#'   ## Get the most significant connections for WW state
#'   dist.ww.m <- do.call(cbind, dist.ww)
#'   simil.ww.m <- max(dist.ww.m) - dist.ww.m
#'   simil.ww.sum <- rowSums(simil.ww.m)
#'   zscores <- (simil.ww.sum - mean(simil.ww.sum)) / sd(simil.ww.sum)
#'   idx <- zscores > z.limit
#'   vertices.ww <- c(rbind(pairs[idx,1], pairs[idx,2]))
#'   ## Get the most significant connections for CC state
#'   dist.cc.m <- do.call(cbind, dist.cc)
#'   simil.cc.m <- max(dist.cc.m) - dist.cc.m
#'   simil.cc.sum <- rowSums(simil.cc.m)
#'   zscores <- (simil.cc.sum - mean(simil.cc.sum)) / sd(simil.cc.sum)
#'   idx <- zscores > z.limit
#'   vertices.cc <- c(rbind(pairs[idx,1], pairs[idx,2]))
#'   ## Get the most significant connections for HET inversion cases
#'   dist.het.m <- do.call(cbind, dist.het)
#'   simil.het.m <- max(dist.het.m) - dist.het.m
#'   simil.het.sum <- rowSums(simil.het.m)
#'   zscores <- (simil.het.sum - mean(simil.het.sum)) / sd(simil.het.sum)
#'   idx <- zscores > z.limit
#'   vertices.het <- c(rbind(pairs[idx,1], pairs[idx,2]))
#'   ## Get cluster ID of putative het INVs
#'   putative.HETs <- unique(vertices.het[duplicated(vertices.het)])
#' 
#'   ## Merge all vertices
#'   #vertices <- c(vertices.ww, vertices.cc)
#'   vertices <- c(vertices.ww, vertices.cc, vertices.wc, vertices.het)
#'   ## Find strongly connected clusters
#'   G <- igraph::graph(vertices, directed = FALSE)
#'   clusters <- igraph::groups(igraph::components(G, mode = 'strong'))
#' 
#'   stopTimedMessage(ptm)
#'   return(list(clusters=clusters, putative.HETs=putative.HETs, max.WC=max.WC))
#' }