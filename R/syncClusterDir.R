#' Re-orient contigs based on shared strand states
#' 
#' This function tries to find set of contigs inverted in respect to each other in every single cell
#' and flips directionality of such 'misoriented' contigs.
#'
#' @param contig.states A \code{data.frame} of strand states per contig and per cell.
#' @importFrom stats hclust dist cutree
#' @author David Porubsky
#' @export
#' 
syncClusterDir <- function(contig.states) {
  ## Remove cell with wc states
  #mask <- apply(contig.states, 2, function(x) all(x != 3))
  #contig.states.sub <- as.matrix(contig.states[,mask])
  contig.states.sub <- contig.states
  ## Set wc states to NA
  contig.states.sub[contig.states.sub == 3] <- NA
  ## Check if non-WC cells have any switch in directionality
  mask <- apply(contig.states.sub, 2, function(x) all(!is.na(x)))
  #mask <- apply(contig.states.sub, 2, function(x) all(x != 3))
  non.WC <- contig.states.sub[,mask, drop=FALSE]
  if (ncol(non.WC) > 2) {
    switches <- apply(non.WC, 2, function(x) any(diff(x) != 0))
    ## Set switch to TRUE if more than half of the non.WC cells contains switch in directionality
    switch <- (length(switches[switches == 'TRUE']) / length(switches)) > 0.5
  } else {
    ## Attempt to orient contigs in case switch in directionality cannot be realiably determined
    switch <- TRUE
  }  
  
  if (switch) {
    ## Remove cell with 'pure' ww or cc state (non-informative)
    mask <- apply(contig.states.sub, 2, function(x) length(unique(x)) > 1)
    contig.states.sub <- as.matrix(contig.states.sub[,mask])
    
    hc.obj <- NULL
    while (class(hc.obj) == 'NULL') {
      ## Check if strand-state matrix can be clustered using Hierarchical clustering
      hc.obj <- tryCatch({
        ## Cluster contigs by hierarchical clustering
        stats::hclust(stats::dist(contig.states.sub))
        }, error = function(e) {return(NULL)}
      )
      ## Remove contig with the maximum number of missing values (NA's)
      if (class(hc.obj) == 'NULL') {
        if (nrow(contig.states.sub) > 2) {
          contig.states.sub <- contig.states.sub[-which.min(rowSums(contig.states.sub, na.rm = TRUE)),]
        } else {
          break
        }  
      }
    }  
    
    if (ncol(contig.states.sub) > 2 & nrow(contig.states.sub) >= 2 & class(hc.obj) == 'hclust') {
      ## Divide antiparallel set of contigs
      hc.clusters <- stats::cutree(hc.obj, k = 2)
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
