#' Import BAM(s) and count reads
#'
#' Import aligned reads from a multiple BAM files and counts directional reads in specified genomic locations.
#' Results are stored in a \code{list} of matrices with each element of a \code{list} representing counts for single BAM file.
#'
#' @param bamfolder A folder containing BAM files with Strand-seq reads aligned to denovo assembly.
#' @param bin.size A length of a bin to count reads in.
#' @param reads.per.bin An approximate number of desired reads per bin. The bin size will be selected accordingly. Forces 'bin.method' to be 'fixed'.
#' @param max.frag A maximum fragment length to import from the BAM file.
#' @param bin.method One of the 'fixed' or 'dynamic' binning method.
#' @param blacklist A \code{\link[GenomicRanges]{GRanges}} object containing genonic regions that should be excluded from the analysis.
#' @return A \code{list} of matrices (columns: minus (W) and plus (C) counts; rows: genomic regions).
#' @importFrom bamsignals bamCount
#' @importFrom Rsamtools BamFile
#' @inheritParams readBamFileAsGRanges
#' @inheritParams makeFixedBins
#' @author David Porubsky
#' @export
#' 
importBams <- function(bamfolder=bamfolder, chromosomes=NULL, pairedEndReads=TRUE, min.mapq=10, bin.size=100000, step.size=NULL, reads.per.bin=NULL, max.frag=1000, bin.method='fixed', blacklist=NULL) {
  ## Get total processing time
  ptm <- proc.time()
  message("Preparing BAM count table ...")
  
  ## Helper functions
  # extendZeroBins <- function(gr, reads.per.bin=50) {
  #   read.cs <- cumsum(gr$total.reads)
  #   min.reads.cs <- seq(from=reads.per.bin, to=max(sum(gr$total.reads), 2*reads.per.bin), by=reads.per.bin)
  #   intervals <- findInterval(read.cs, min.reads.cs, all.inside = TRUE)
  #   gr$group <- intervals
  #   gr.concat <- collapseBins(gr, id.field = 4, measure.field = c(1,2,3))
  #   return(gr.concat[,-4])
  # } 
  
  ## List bams present in a directory
  bamfiles <- list.files(bamfolder, pattern = '.bam$', full.names = T)
  
  ### Get chromosome lengths ###
  bamfile <- bamfiles[1] ## process a first bamfile in a list
  if (!is.null(bamfile)) {
    chrom.lengths <- GenomeInfoDb::seqlengths(Rsamtools::BamFile(bamfile))
  }
  chroms.in.data <- names(chrom.lengths)
  if (is.null(chromosomes)) {
    chromosomes <- chroms.in.data
  }
  chroms2use <- intersect(chromosomes, chroms.in.data)
  ## Stop if none of the specified chromosomes exist
  if (length(chroms2use) == 0) {
    chrstring <- paste0(chromosomes, collapse=', ')
    stop('Could not find length information for any of the specified chromosomes: ', chrstring, '. Pay attention to the naming convention in your data, e.g. "chr1" or "1".')
  }
  ## Issue warning for non-existent chromosomes
  diff <- setdiff(chromosomes, chroms.in.data)
  if (length(diff) > 0) {
    diffs <- paste0(diff, collapse=', ')
    warning('Could not find length information for the following chromosomes: ', diffs)
  }
  chrom.lengths <- chrom.lengths[chroms2use]
  
  ## Make chromosome/contig ranges
  chr.gr <- GenomicRanges::GRanges(seqnames=chroms2use, ranges=IRanges(start=1, end=chrom.lengths))
  
  ## Make genome bins
  if (!is.null(bin.size) & bin.method == 'fixed') {
    bins.gr <- makeFixedBins(bamfile = bamfile, bin.size = bin.size, step.size = step.size, chromosomes = chroms2use, keep.small.chr = TRUE)
  } else if (!is.null(bin.size) & bin.method == 'dynamic') {
    bins.gr <- makeDynamicBins(bamfiles = bamfiles, bin.size = bin.size, step.size = step.size, chromosomes = chroms2use, keep.small.chr = TRUE)
  } else {
    warning("Unsupported binning method!!!, Set 'bin.method' to 'fixed' or 'dynamic'")
  }
  
  ## Set parameter for bamsignals counts
  paired.end <- 'ignore'
  if (pairedEndReads) {
    paired.end <- 'filter'
  }
  
  ## Mask regions with and excess of read coverage that appears always WC
  ## and bins that have very low read counts
  if (!is.null(blacklist)) {
    #removed.gr <- subsetByOverlaps(bins.gr, blacklist)
    bins.gr <- subsetByOverlaps(bins.gr, blacklist, invert = TRUE)
    bins.gr <- keepSeqlevels(bins.gr, value = as.character(unique(seqnames(bins.gr))), pruning.mode = 'coarse')
  }
  
  counts.l <- list()
  for (i in 1:length(bamfiles)) {
    bam <- bamfiles[i]
    bam.name <- basename(bam)
    message("    Processing ", bam.name)
    
    ## Scale bin size to the required minimum number of reads in a bin
    if (!is.null(reads.per.bin)) {
      bin.method <- 'fixed'
      chr.counts <- bamsignals::bamCount(bam, chr.gr, mapq=min.mapq, filteredFlag=1024, paired.end=paired.end, tlenFilter=c(0, max.frag), verbose=FALSE)
      n.reads <- as.numeric(sum(chr.counts))
      num.counts.perbp <- n.reads / sum(as.numeric(chrom.lengths))
      bin.size <- round(reads.per.bin / num.counts.perbp, -2)
      ## Make genome bins
      bins.gr <- makeFixedBins(bamfile = bam, bin.size = bin.size, step.size = step.size, chromosomes = chroms2use)
    }
    
    if (length(bins.gr) > 0) {
      ## Get read counts per bin
      counts <- bamsignals::bamCount(bam, bins.gr, mapq=min.mapq, filteredFlag=1024, paired.end=paired.end, tlenFilter=c(0, max.frag), verbose=FALSE, ss=TRUE)
      mcols(bins.gr) <- t(counts)
      bins.gr$total.reads <- colSums(counts)
    } else {
      ## Get read counts per chromosome/contig
      counts <- bamsignals::bamCount(bam, chr.gr, mapq=min.mapq, filteredFlag=1024, paired.end=paired.end, tlenFilter=c(0, max.frag), verbose=FALSE, ss=TRUE)
      mcols(chr.gr) <- t(counts)
      chr.gr$total.reads <- colSums(counts)
    }
    
    ## Extend bins with zero read counts
    # if(!is.null(reads.per.bin)) {
    #   bins.grl <- split(bins.gr, seqnames(bins.gr))
    #   bins.grl <- endoapply(bins.grl, function(x) extendZeroBins(gr=x, reads.per.bin=config[['reads.per.bin']]))
    #   bins.gr <- unlist(bins.grl, use.names = FALSE)
    # }
    
    ## Get continous groups of zero counts
    # zero.bins <- which(bins.gr$total.reads < min.reads)
    # groups <- cumsum(c(1, abs(zero.bins[-length(zero.bins)] - zero.bins[-1]) > 1))

    ## Create a count matrix
    bins.gr.new <- bins.gr
    counts.m <- cbind(bins.gr.new$antisense, bins.gr.new$sense)
    ## Extend gaps between ranges
    bins.gr.ext <- expandGaps(bins.gr)
    rownames(counts.m) <- as.character(bins.gr.ext)
    
    counts.l[[bam.name]] <- counts.m
  }
  time <- proc.time() - ptm
  message("\nTime spent: ", round(time[3],2), "s")
  
  return(counts.l)
}
