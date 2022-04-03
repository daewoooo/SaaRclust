#' A function to split FASTA file into a separate files one per FASTA sequence.
#' 
#' @param fasta.file A FASTA file to be divided into individual FASTA sequences.
#' @param outputfolder A path to a folder where the output will be written.
#' @importFrom Rsamtools FaFile scanFa
#' @importFrom Biostrings writeXStringSet
#' @author David Porubsky
#' @export
#'
splitFASTA <- function(fasta.file=NULL, outputfolder=NULL) {
  ## Extract FASTA from user defined FASTA file
  fa.file <- open(Rsamtools::FaFile(fasta.file))
  ## Read in contigs for a given cluster
  gr.seq <- Rsamtools::scanFa(file = fa.file, as = "DNAStringSet")
  
  if (is.null(outputfolder)) {
    return(gr.seq)
  } else {
    if (!dir.exists(outputfolder)) {
      dir.create(outputfolder)
    }
    ## Loop over each FASTA and print into a separate file
    for (i in seq_along(gr.seq)) {
      fasta.seq <- gr.seq[i]
      ## Print out fasta
      fasta.name <- paste0(names(fasta.seq), '.fasta')
      fasta.save <- file.path(outputfolder, fasta.name)
      Biostrings::writeXStringSet(x = fasta.seq, filepath = fasta.save, format = 'fasta')
    }
  }
}