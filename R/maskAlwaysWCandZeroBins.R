#' This function takes list of BAM files and enumerates user defined genomic bins 
#' for recurrent regions of WC strand-state and low read counts. Such bins are reported
#' as a custom \code{\link[GenomicRanges]{GRanges}} object or regions to be blacklisted.
#' 
#' @param genomic.bins A \code{\link[GenomicRanges]{GRanges}} object with defined genomic bins.
#' @inheritParams importBams
#' @importFrom bamsignals bamCount
#' @author David Porubsky
#' @export
#' 
maskAlwaysWCandZeroBins <- function(bamfolder=bamfolder, genomic.bins=NULL, min.mapq=10, pairedEndReads=TRUE, max.frag=1000) {
  ptm <- startTimedMessage("Detecting always WC and low read count regions")
  
  ## List bams present in a directory
  bamfiles <- list.files(bamfolder, pattern = '.bam$', full.names = T)
  
  ## Set parameter for bamsignals counts
  paired.end <- 'ignore'
  if (pairedEndReads) {
    paired.end <- 'filter'
  }
  
  ## Go over all Strand-seq libraries and report strand-state for each bin
  wc.counts <- list()
  total.read.counts <- list()
  for (i in 1:length(bamfiles)) {
    bam <- bamfiles[i]
    counts <- suppressMessages( bamsignals::bamCount(bam, genomic.bins, 
                                                     mapq=min.mapq, 
                                                     filteredFlag=1024, 
                                                     paired.end=paired.end, 
                                                     tlenFilter=c(0, max.frag), 
                                                     verbose=FALSE, 
                                                     ss=TRUE) 
    )
    genoT <- countProb(minusCounts = counts[2,], plusCounts = counts[1,], log = TRUE)
    genoT <- apply(genoT, 1, which.max)
    wc.genoT <- rep(0, length(genoT))
    wc.genoT[genoT == 3] <- 1
    wc.counts[[i]] <- wc.genoT
    
    total.read.counts[[i]] <- colSums(counts)
  }
  wc.counts.m <- do.call(cbind, wc.counts)
  wc.counts.sums <- rowSums(wc.counts.m)
  total.read.counts.m <- do.call(cbind, total.read.counts)
  total.read.sums <- rowSums(total.read.counts.m)
  
  ## Find bins that are always WC in a majority of cells
  ## Remove bins that are WC in more than 90% of all cells (bamfiles)
  z.score <- (wc.counts.sums - median(wc.counts.sums)) / sd(wc.counts.sums)
  mask.bins <- which(z.score >= 3)
  #thresh <- length(bamfiles) * 0.90
  #mask.bins <- which(wc.counts.sums >= thresh)
  if (length(mask.bins) > 0) {
    alwaysWC <- genomic.bins[mask.bins]
  } else {
    alwaysWC <- NULL
  } 
  
  ## Find bins that have close to zero counts
  z.score <- (total.read.sums - median(total.read.sums)) / sd(total.read.sums)
  zero.bins <- which(z.score <= -3)
  if (length(zero.bins) > 0) {
    alwaysZero <- genomic.bins[zero.bins]
  } else {
    alwaysZero <- NULL
  } 
  
  stopTimedMessage(ptm)
  return(list(alwaysWC=alwaysWC, alwaysZero=alwaysZero))
}
