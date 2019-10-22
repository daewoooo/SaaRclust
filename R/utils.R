#' Collapse consecutive bins with the same ID value
#'
#' Collapse consecutive bins with the same value defined in 'id.field'. 
#'
#' @param gr A \code{\link{GRanges-class}} object.
#' @param id.field A number of metadata column to use for region merging.
#' @param measure.field A field column that contains measured value to be sum up.
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
getCoinheritanceMatrix <- function(contig.states = NULL, scale=FALSE, WWtoCC.penalty=TRUE) {
  ## Initialize empty matrix
  n.contigs <- nrow(contig.states)
  coinherit.m <- matrix(nrow = n.contigs, ncol = n.contigs)
  pairs <- as.matrix(expand.grid(1:n.contigs, 1:n.contigs))
  for (i in 1:nrow(pairs)) {
    pair <- pairs[i,]
    comp <- contig.states[pair[1],] == contig.states[pair[2],]
    coinherit.score <- length(comp[comp == TRUE])
    if (WWtoCC.penalty) {
      ## Penalize cases where there is W and C state in a single cell
      WWtoCC <- abs(contig.states[pair[1],] - contig.states[pair[2],]) == 1
      WWtoCC.score <- length(WWtoCC[WWtoCC == TRUE]) #/ ncol(contig.states)
      ## Assign coinherit.score
      coinherit.m[pair[1], pair[2]] <- coinherit.score - WWtoCC.score
    } else {
      coinherit.m[pair[1], pair[2]] <- coinherit.score
    }  
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

#' Convert strings in a data.frame object into a \code{\link{GRanges-class}} object.
#'
#' @param region.string A \code{vector} object 
#' @return A \code{\link{GRanges-class}} object with an extra metacolumn group.ID.
#' @importFrom tidyr separate
#' @author David Porubsky
#' @export
#' 
string2GRanges <- function(region.string=NULL) {
  regions.df <- data.frame(regions=region.string)
  regions.df <- tidyr::separate(data = regions.df, col = regions, into = c('chr','posANDdir'), sep = ":")
  regions.df <- tidyr::separate(data = regions.df, col = posANDdir, into = c('start','end','dir'), sep = "-|_")
  regions.gr <- GenomicRanges::GRanges(seqnames=regions.df$chr, 
                                       ranges=IRanges(start=as.numeric(regions.df$start), 
                                                      end=as.numeric(regions.df$end)), 
                                       dir=regions.df$dir)
  return(regions.gr)
}


#' Fill in gaps in between a set of genomic ranges.
#' 
#' This function takes a \code{\link{GRanges-class}} object and extend gaps between
#' subsequent ranges by merging with the preceeding or following genomic range.
#' 
#' @param gr A \code{\link{GRanges-class}} object.
#' @return A \code{\link{GRanges-class}} object.
#' @author David Porubsky
#' @export

expandGaps <- function(gr) {
  ## Check if seqlengths are defined
  if (any(is.na(seqlengths(gr)))) {
    stop("Undefined seqlengths for parameter gr ...")
  }
  
  ## Extend gaps between ranges
  new.gr <- gr
  gaps.gr <- GenomicRanges::gaps(GenomicRanges::sort(gr))
  gaps.gr <- gaps.gr[strand(gaps.gr) == '*']
  if (length(gaps.gr) > 0) {
    ## Merge with preceeding range
    preceed.idx <- GenomicRanges::start(gaps.gr) > 1
    GenomicRanges::start(gaps.gr)[preceed.idx] <- GenomicRanges::start(gaps.gr)[preceed.idx] - 1
    preceed.hits <- GenomicRanges::findOverlaps(gaps.gr[preceed.idx], new.gr)
    GenomicRanges::end(new.gr[subjectHits(preceed.hits)]) <- GenomicRanges::end(gaps.gr[preceed.idx])
    ## Merge with following range
    follow.idx <- GenomicRanges::start(gaps.gr) == 1
    GenomicRanges::end(gaps.gr)[follow.idx] <- GenomicRanges::end(gaps.gr)[follow.idx] + 1
    follow.hits <- GenomicRanges::findOverlaps(gaps.gr[follow.idx], new.gr)
    GenomicRanges::start(new.gr[subjectHits(follow.hits)]) <- GenomicRanges::start(gaps.gr[follow.idx])
    return(new.gr)
  } else {
    return(gr)
  }
}
