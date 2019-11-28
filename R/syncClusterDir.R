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
  #mask <- apply(contig.states, 2, function(x) all(x != 3))
  #contig.states.sub <- as.matrix(contig.states[,mask])
  ## Set wc states to NA
  contig.states.sub <- contig.states
  contig.states.sub[contig.states.sub == 3] <- NA
  ## Check if non-WC cells have any switch in directionality
  mask <- apply(contig.states.sub, 2, function(x) all(!is.na(x)))
  non.WC <- contig.states.sub[,mask, drop=FALSE]
  if (ncol(non.WC) > 2) {
    switch <- any(apply(non.WC, 2, function(x) any(diff(x) != 0)))
  } else {
    ## Attempt to orient contigs in case switch in directionality cannot be realiably determined
    switch <- TRUE
  }  
  
  if (switch) {
    ## Remove cell with 'pure' ww or cc state (non-informative)
    mask <- apply(contig.states.sub, 2, function(x) length(unique(x)) > 1)
    contig.states.sub <- as.matrix(contig.states.sub[,mask])
    ## Check if strand-state matrix can be clustered using HC
    hc.obj <- tryCatch({
      ## Cluster contigs by hierarchical clustering
      hclust(dist(contig.states.sub))
      }, error = function(e) {return('error')}
    )
    
    if (ncol(contig.states.sub) > 2 & nrow(contig.states.sub) > 2 & class(hc.obj) == 'hclust') {
      ## Divide antiparallel set of contigs
      hc.clusters <- cutree(hc.obj, k = 2)
      misorients <- split(names(hc.clusters), hc.clusters)
      ## Flip WW and CC states in the smaller group of misorients
      contigsToFlip <- misorients[[which.min(lengths(misorients))]]
      contigsToFlip.newID <- paste0(contigsToFlip, "_revcomp")
      idx.toFlip <- which(rownames(contig.states) %in% contigsToFlip)
      contig.states.recoded <- contig.states[idx.toFlip,]
      contig.states.recoded[contig.states[idx.toFlip,] == 1] <- 2
      contig.states.recoded[contig.states[idx.toFlip,] == 2] <- 1
      contig.states[idx.toFlip,] <- contig.states.recoded
      ## Rename flipped contigs
      rownames(contig.states)[idx.toFlip] <- contigsToFlip.newID
      rownames(contig.states)[-idx.toFlip] <- paste0(rownames(contig.states)[-idx.toFlip], "_dir")
      ## Return re-oriented matrix  
      return(contig.states)
    } else {
      message("    Contig orienting failed!!!")
      rownames(contig.states) <- paste0(rownames(contig.states), "_dir")
      return(contig.states)
    }
  } else {
    message("    No switch detected, skipping!!!")
    rownames(contig.states) <- paste0(rownames(contig.states), "_dir")
    return(contig.states)
  }  
}
