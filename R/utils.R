#' Collapse consecutive bins with the same ID value
#'
#' Collapse consecutive bins with the same value defined in 'id.field'. 
#'
#' @param gr A \code{\link{GRanges-class}} object.
#' @param id.field A number of metadata column to use for region merging.
#' @return A \code{\link{GRanges-class}} object.
#' @author David Porubsky
#' @export
#' 
collapseBins <- function(gr, id.field=0, measure.field=NULL) {
  ## Include seqnames into the unique ID
  unique.ID <- paste0(seqnames(gr), '_', mcols(gr)[,id.field])
  ## Get continous runs of the same unique ID
  unique.ID.runLength <- runLength(Rle(unique.ID))
  ind.last <- cumsum(unique.ID.runLength) ##get indices of last range in a consecutive(RLE) run of the same value
  ind.first <- c(1,cumsum(unique.ID.runLength) + 1) ##get indices of first range in a consecutive(RLE) run of the same value
  ind.first <- ind.first[-length(ind.first)]  ##erase last index from first range indices 
  collapsed.gr <- GenomicRanges::GRanges(seqnames=seqnames(gr[ind.first]), ranges=IRanges(start=start(gr[ind.first]), end=end(gr[ind.last])), mcols=mcols(gr[ind.first]))
  names(mcols(collapsed.gr)) <- names(mcols(gr[ind.first]))
  ## Sum values in user defined measure fields(s)
  if (!is.null(measure.field)) {
    run.ID  <- factor(rep(seq_along(unique.ID.runLength), unique.ID.runLength))
    for (field in measure.field) {
      mcols(collapsed.gr)[,field] <-  tapply(mcols(gr)[,field], run.ID, sum)
      #collapsed.gr$C <-  tapply(gr$C, run.ID, sum)
      #collapsed.gr$W <-  tapply(gr$W, run.ID, sum)
    }
  }
  return(collapsed.gr)
}


#' Reverse orientation of directional read counts
#'
#' This function flips the orientation of directional reads counts for each genomic position.
#'
#' @param counts.l A \code{list} of plus and minus alignments per genomic region.
#' @return A \code{list} of matrices (columns: minus (W) and plus (C) counts; rows: genomic regions).
#' @export
#' 
flipCounts <- function(counts.l) {
  flipped.counts.l <- list()
  for (i in 1:length(counts.l)) {
    toFlip <- counts.l[[i]]
    toFlip[,1] <- counts.l[[i]][,2]
    toFlip[,2] <- counts.l[[i]][,1]
    rownames(toFlip) <- paste0(rownames(toFlip),"_rev")
    flipped.counts.l[[i]] <- toFlip
  }
  return(flipped.counts.l)
}


#' Function to obtain co-inheritance matrix based on shared strand states
#' 
#' This function takes \code{data.frame} of strand states per contig [rows] and per cell [columns]
#' and converts it into a matrix of contig similarities based on co-inheritance patterns.
#'
#' @param contig.states A \code{data.frame} of strand states per contig and per cell.
#' @param scale If set to \code{TRUE} returned values will be scaled to 1.
#' @return A symetric \code{matrix} of similarities for each contig pair.
#' @author David Porubsky
#' @export
#' 
getCoinheritanceMatrix <- function(contig.states = NULL, scale=FALSE) {
  ## Initialize empty matrix
  n.contigs <- nrow(contig.states)
  coinherit.m <- matrix(nrow = n.contigs, ncol = n.contigs)
  pairs <- as.matrix(expand.grid(1:n.contigs, 1:n.contigs))
  for (i in 1:nrow(pairs)) {
    pair <- pairs[i,]
    comp <- contig.states[pair[1],] == contig.states[pair[2],]
    coinherit.score <- length(comp[comp == TRUE])
    coinherit.m[pair[1], pair[2]] <- coinherit.score
  }
  ## Convert coinheritance to similarity
  coinherit.m <- max(coinherit.m) - coinherit.m
  ## Scale similarities to 1
  if (scale) {
    coinherit.m <- coinherit.m / sum(coinherit.m)
  }
  ## Add row and column names and return
  rownames(coinherit.m) <- rownames(contig.states)
  colnames(coinherit.m) <- rownames(contig.states)
  return(coinherit.m)
}


#' Assign group.ID for similar cluster.IDs
#' 
#' This function will write extra metacolumn 'group.ID' into an input \code{\link{GRanges-class}} object
#' based on cluster.groups parameter.
#'
#' @param cluster.gr A \code{\link{GRanges-class}} object with contig position and their cluster assignment in 'clust.ID' metacolumn.
#' @param cluster.groups A \code{list} of grouped 'clust.IDs'.
#' @return A \code{\link{GRanges-class}} object with an extra metacolumn group.ID.
#' @author David Porubsky
#' @export
#' 
addClusterGroup <- function(cluster.gr=NULL, cluster.groups=NULL) {
  ## Initialize cluster group with clust.ID
  cluster.gr$group.ID <- cluster.gr$clust.ID
  ## Add cluster group for pair of clusters specified in cluster.pairs
  for (i in seq_along(cluster.groups)) {
    clust.group <- cluster.groups[[i]]
    ## Assign common cluster group for a pair of clusters
    cluster.gr.sub <- cluster.gr[cluster.gr$clust.ID %in% clust.group]
    cluster.gr.sub$group.ID <- paste0('cluster',i)
    cluster.gr[cluster.gr$clust.ID %in% clust.group] <- cluster.gr.sub
  }
  return(cluster.gr)
} 