#' Re-orient contigs based on shared strand states
#' 
#' This function tries to find set of contigs inverted in respect to each other in every single cell
#' and flips directionality of such 'misoriented' contigs.
#'
#' @param contig.states A \code{data.frame} of strand states per contig and per cell.
#' @author David Porubsky
#' @export
#' 
syncClusterDir <- function(contig.states) {
  ## Remove cell with wc states
  mask <- apply(contig.states, 2, function(x) all(x != 3))
  contig.states.sub <- as.matrix(contig.states[,mask])
  ## Remove cell with 'pure' ww or cc state (non-informative)
  mask <- apply(contig.states.sub, 2, function(x) length(unique(x)) > 1)
  contig.states.sub <- as.matrix(contig.states.sub[,mask])
  if (ncol(contig.states.sub) > 2 & nrow(contig.states.sub) > 2) {
    ## Divide antiparallel set of contigs by hierarchical clustering
    hc.obj <- hclust(dist(contig.states.sub))
    hc.clusters <- cutree(hc.obj, k = 2)
    misorients <- split(names(hc.clusters), hc.clusters)
    ## Flip WW and CC states in the smaller group of misorients
    contigsToFlip <- misorients[[which.min(lengths(misorients))]]
    contigsToFlip.newID <- paste0(contigsToFlip, "_revcomp")
    idx.toFlip <- which(rownames(contig.states) %in% contigsToFlip)
    contig.states.recoded <- cluster.m[idx.toFlip,]
    contig.states.recoded[cluster.m[idx.toFlip,] == 1] <- 2
    contig.states.recoded[cluster.m[idx.toFlip,] == 2] <- 1
    contig.states[idx.toFlip,] <- contig.states.recoded
    ## Rename flipped contigs
    rownames(contig.states)[idx.toFlip] <- contigsToFlip.newID
    rownames(contig.states)[-idx.toFlip] <- paste0(rownames(contig.states)[-idx.toFlip], "_dir")
    ## Returnet re-oriented matrix  
    return(contig.states)
  } else {
    rownames(contig.states) <- paste0(rownames(contig.states), "_dir")
    return(contig.states)
  }
}
