#' Make fixed-width bins
#'
#' Make fixed-width bins based on given bin and step size.
#'
#' @param bamfile A BAM file from which the header is read to determine the chromosome lengths. If a \code{bamfile} is specified, option \code{assembly} is ignored.
#' @param bin.size A size of the genomic bin to split genome into.
#' @param step.size A size of the genomic interval to move each bin.
#' @param chromosomes A user defined set of chromosomes for binning (eg. 'chr1')
#' @return A \code{\link{GRanges-class}} object with fixed-width bins.
#' @author David Porubsky
#' @importFrom Rsamtools BamFile
#' @export
#'
makeBins <- function(bamfile=NULL, bin.size=100000, step.size=NULL, chromosomes=NULL) {
  
  ptm <- startTimedMessage("Preparing fixed-width bins for bin size ", bin.size, " ...")
  
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
  bins <- bins[end(bins) > 0] # remove chromosomes that are smaller than bin.size
  if (any(width(bins) != bin.size)) {
    stop("tileGenome failed")
  }
  ## Add sequence lengths
  seqlengths(bins) <- chrom.lengths[chroms2use]
  ## Add step size if defined
  if (!is.null(step.size)) {
    shift.bp <- 0
    bins.list.step <- GRangesList()
    while (shift.bp < bin.size) {
      bins.list.step[[as.character(shift.bp)]] <- suppressWarnings( trim(shift(bins, shift.bp)) )
      shift.bp <- step.size + shift.bp
    }
    bins <- sort(unlist(bins.list.step, use.names = FALSE))  
  }
  ## Report chromosome/contigs smaller than the bin.size
  skipped.chroms <- setdiff(seqlevels(bins), as.character(unique(seqnames(bins))))
  if (length(skipped.chroms) > 0) {
    warning("The following chromosomes/contigs were skipped because they are smaller than bin.size ", bin.size, ": ", paste0(skipped.chroms, collapse=', '))
  }
  
  stopTimedMessage(ptm)
  return(bins)
}
