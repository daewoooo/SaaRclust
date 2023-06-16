#' Export FASTA sequences from a set of Genomic Ranges.
#'
#' This function takes a \code{\link{GRanges-class}} object and extracts a genomic sequence 
#' from these regions either from an original range or a range expanded on each side by
#' defined number of bases.
#' 
#' @param gr A \code{\link{GRanges-class}} object of genomic regions to extract genomic sequence from. 
#' @param bsgenome A \code{\link{BSgenome-class}} object of reference genome to get the genomic sequence from. 
#' @param asm.fasta An assembly FASTA file to extract DNA sequence determined by 'gr' parameter.
#' @param index.field A user defined column number to be used as a header id for the exported FASTA.
#' @param expand Expand ends of each genomic region by this length (in bp).
#' @param fasta.save A path to a filename where to store final FASTA file.
#' @importFrom Rsamtools indexFa FaFile scanFa scanFaIndex
#' @importFrom BSgenome getSeq
#' @importFrom Biostrings writeXStringSet
#' @author David Porubsky
#' @export
#'
regions2FASTA <- function(gr, bsgenome=NULL, asm.fasta=NULL, index.field=NULL, expand=0, fasta.save=NULL) {
  ## Load BSgenome object
  if (class(bsgenome) != 'BSgenome') {
    if (is.character(bsgenome)) {
      bsgenome <- tryCatch({
        suppressPackageStartupMessages(library(bsgenome, character.only=TRUE))
        bsgenome <- eval(parse(text=bsgenome)) ## replacing string by object
      }, error = function(e) {return(NULL)})
    }
  }
  
  if (!is.null(asm.fasta)) {	
    ## Check if submitted fasta file is indexed
    asm.fasta.idx <- paste0(asm.fasta, ".fai")
    if (!file.exists(asm.fasta.idx)) {
      fa.idx <- Rsamtools::indexFa(file = asm.fasta)
    }
  }
  
  ## Expand regions of interest by certain size (downstream and upstream)
  if (expand > 0) {
    gr.seqLen <- seqlengths(gr)[as.character(seqnames(gr))]
    if (!is.na(gr.seqLen)) {
      start(gr) <- pmax(1, start(gr) - expand)
      end(gr) <- pmin(gr.seqLen, end(gr) + expand)
    } else {
      warning("Could not expand the exported regions, missing 'seqlengths' in submitted 'gr' object!")
    }  
  }
  
  if (!is.null(bsgenome)) {
    ## Extract FASTA from BSgenome object
    gr.seq <- BSgenome::getSeq(bsgenome, gr)
    names(gr.seq) <- as.character(gr)
  } else if (is.character(asm.fasta)) {
    ## Extract FASTA from user defined FASTA file
    fa.file <- open(Rsamtools::FaFile(asm.fasta))
    ## Remove sequences not present in the FASTA index
    fa.idx <- Rsamtools::scanFaIndex(fa.file)
    gr <- suppressWarnings( subsetByOverlaps(gr, fa.idx) )
    ## Read in contigs for a given cluster
    gr.seq <- Rsamtools::scanFa(file = fa.file, param = gr, as = "DNAStringSet")
    names(gr.seq) <- gsub(as.character(gr), pattern = ':', replacement = '_')
  } else {
    stop("Please set a 'bsgenome' or 'asm.fasta' parameter!!!")
  }
  
  ## Used user defined column to name FASTA sequences
  if (!is.null(index.field)) {
    if (index.field <= length(GenomicRanges::mcols(gr))) {
      names(gr.seq) <- as.character(GenomicRanges::mcols(gr)[[index.field]])
    }
  }
  
  ## Write final FASTA
  if (is.character(fasta.save)) {
    Biostrings::writeXStringSet(x = gr.seq, filepath = fasta.save, format = 'fasta')
  } else {
    warning("Please speficify 'fasta.save' if you want to export FASTA into a file!!!")
  }  
  return(gr.seq)
}
