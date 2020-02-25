#' Make fixed-width bins
#'
#' Make fixed-width bins based on given bin and step size.
#'
#' @param bamfile A BAM file from which the header is read to determine the chromosome lengths. If a \code{bamfile} is specified, option \code{assembly} is ignored.
#' @param bin.size A size of the genomic bin to split genome into.
#' @param step.size A size of the genomic interval to move each bin.
#' @param chromosomes A user defined set of chromosomes for binning (eg. 'chr1')
#' @param keep.small.chr Set to \code{TRUE} if chromosome/contigs smaller than the 'bin.size' should be kept.
#' @return A \code{\link{GRanges-class}} object with fixed-width bins.
#' @author David Porubsky
#' @importFrom Rsamtools BamFile
#' @export
#'
makeFixedBins <- function(bamfile=NULL, bin.size=100000, step.size=NULL, chromosomes=NULL, keep.small.chr=FALSE) {
  
  ptm <- startTimedMessage("Preparing fixed-width bins for bin size ", bin.size)
  
  ### Check user input ###
  if (length(bin.size) == 0) {
    return(list())
  }
  if (is.null(bamfile)) {
    stop("Please specify a 'bamfile'!!!")
  }
  
  ### Get chromosome lengths ###
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
  
  chrom.lengths.floor <- floor(chrom.lengths / bin.size) * bin.size
  bins <- unlist(GenomicRanges::tileGenome(chrom.lengths.floor[chroms2use], tilewidth=bin.size), use.names=FALSE)
  ## remove chromosomes that are smaller than bin.size?
  bins <- bins[end(bins) > 0]
  
  if (any(width(bins) != bin.size)) {
    stop("tileGenome failed")
  }
  
  ## Add sequence lengths
  seqlengths(bins) <- chrom.lengths[chroms2use]
  
  ## Add step size if defined
  if (!is.null(step.size)) {
    shift.bp <- 0
    bins.list.step <- GenomicRanges::GRangesList()
    while (shift.bp < bin.size) {
      bins.list.step[[as.character(shift.bp)]] <- suppressWarnings( trim(GenomicRanges::shift(bins, shift.bp)) )
      shift.bp <- step.size + shift.bp
    }
    bins <- sort(unlist(bins.list.step, use.names = FALSE))  
  }
  
  ## Keep chromosomes/contigs smaller than the bin.size if desired
  if (keep.small.chr) {
    chroms2keep <- chroms2use[!chroms2use %in% unique(seqnames(bins))]
    if (length(chroms2keep) > 0) {
      chroms2keep.gr <- GenomicRanges::GRanges(seqnames = chroms2keep, ranges=IRanges(start = 1, end = chrom.lengths[chroms2keep]))
      bins <- GenomicRanges::sort(c(bins, chroms2keep.gr))
    }  
  }
  
  ## Report chromosome/contigs smaller than the bin.size
  skipped.chroms <- setdiff(seqlevels(bins), as.character(unique(seqnames(bins))))
  if (length(skipped.chroms) > 0) {
    warning("The following chromosomes/contigs were skipped because they are smaller than bin.size ", bin.size, ": ", paste0(skipped.chroms, collapse=', '))
  }
  
  stopTimedMessage(ptm)
  return(bins)
}


#' Make dynamic-width bins
#'
#' Make dynamic-width bins based on given bin and step size.
#' NOTE: Here 'bin.size' represents number of mappable positions in each bin !!!
#'
#' @inheritParams makeFixedBins
#' @return A \code{\link{GRanges-class}} object with fixed-width bins.
#' @author David Porubsky
#' @importFrom Rsamtools BamFile
#' @importFrom bamsignals bamCoverage
#' @export
#'
makeDynamicBins <- function(bamfiles=NULL, bin.size=100000, step.size=NULL, chromosomes=NULL, keep.small.chr=FALSE) {
  
  ## Helper function
  reformat <- function(x) {
    out_list <- list() 
    for ( i in seq(1, length(x), 2) ) {
      out_list[[i]] <- c(x[i], x[i+1])
    }
    mt <- do.call("rbind",out_list)
    df <- data.frame(mt)
    colnames(df) <- c("start", "end")
    return(df)
  } 
  
  ptm <- startTimedMessage("Preparing dynamic-width bins for bin size ", bin.size)
  
  ### Check user input ###
  if (length(bin.size) == 0) {
    return(list())
  }
  if (is.null(bamfiles)) {
    stop("Please specify list of 'bamfiles'!!!")
  }
  
  ### Get chromosome lengths ###
  bamfile <- bamfiles[1]
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
  ## Make chromosome/contig ranges
  chr.gr <- GenomicRanges::GRanges(seqnames=chroms2use, ranges=IRanges(start=1, end=chrom.lengths[chroms2use]))
  
  ## Make variable genome bins
  bins.l <- GenomicRanges::GRangesList()
  for (i in seq_along(chr.gr)) {
    chrom.gr <- chr.gr[i]
    chr.cov <- rep(0, width(chrom.gr))
    chr <- as.character(seqnames(chrom.gr))
    #message("Processing :", chr)
    for (j in 1:length(bamfiles)) {
      bam <- bamfiles[j]
      bam.name <- basename(bam)
      
      ## Check if the bamfile is indexed
      bamindex <- paste0(bam,'.bai')
      if (!file.exists(bamindex)) {
        warning("Couldn't find BAM index-file, indexing ...")
        bamindex.own <- Rsamtools::indexBam(bam)
      }
      
      #message("    Processing ", bam.name)
      chr.counts <- bamsignals::bamCoverage(bam, gr = chrom.gr, mapq=10, filteredFlag=1024, paired.end='ignore', verbose=FALSE)
      chr.cov <- chr.cov + chr.counts[1]
    }
    #chr.cov.binned <- zoo::rollapply(chr.cov, width=bin.size, sum, by=step.size)
    chr.cov[chr.cov > 0] <- 1
    cs <- cumsum(chr.cov)
    if (bin.size >= max(cs)) {
      bin.starts <- 1
      bin.ends <- max(cs)
      bins <- c(rbind(bin.starts, bin.ends)) 
    } else {
      bin.starts <- seq(from = 1, to = max(cs) - bin.size, by = step.size)
      bin.ends <- seq(from = bin.size, to = max(cs), by = step.size)
      bins <- c(rbind(bin.starts, bin.ends)) 
    }  
    
    gen.pos <- findInterval(bins, cs, rightmost.closed = T)
    gen.bins <- reformat(gen.pos)
    gen.bins.gr <- GenomicRanges::GRanges(seqnames = chr, ranges=IRanges(start=gen.bins$start, end=gen.bins$end))
    
    bins.l[[i]] <- gen.bins.gr
  }  
  bins <- unlist(bins.l, use.names = FALSE)
  
  ## Add sequence lengths
  seqlengths(bins) <- chrom.lengths[chroms2use]
  
  ## Keep chromosomes/contigs smaller than the bin.size if desired
  if (keep.small.chr) {
    chroms2keep <- chroms2use[!chroms2use %in% unique(seqnames(bins))]
    if (length(chroms2keep) > 0) {
      chroms2keep.gr <- GRanges(seqnames = chroms2keep, ranges=IRanges(start = 1, end = chrom.lengths[chroms2keep]))
      bins <- sort(c(bins, chroms2keep.gr))
    }  
  }
  
  ## Report chromosome/contigs smaller than the bin.size
  skipped.chroms <- setdiff(GenomeInfoDb::seqlevels(bins), as.character(unique(seqnames(bins))))
  if (length(skipped.chroms) > 0) {
    warning("The following chromosomes/contigs were skipped because they are smaller than bin.size ", bin.size, ": ", paste0(skipped.chroms, collapse=', '))
  }
  
  stopTimedMessage(ptm)
  return(bins)
}