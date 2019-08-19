#' Concatenate FASTA sequences given a set of Genomic Ranges and cluster identity.
#'
#' This function takes a \code{\link{GRanges-class}} object and cluster assignment stored 
#' in 'ID' field and export concated FASTA file for each cluster and genomic region/contig 
#' assigned to a given cluster.
#' 
#' @param clustered.gr A \code{\link{GRanges-class}} object with given fields: dir - contig directionality
#' order - contig order and ID - contig assignment to clusters.
#' @param assembly.fasta An assembly FASTA file to extract DNA sequence determined by clustered.gr.
#' @param outputfolder A path to a folder to export concatenated FASTA sequences.
#' @importFrom Rsamtools indexFa FaFile scanFa
#' @importFrom Biostrings reverseComplement DNAStringSet writeXStringSet
#' @author David Porubsky
#' @export
exportPseudoChromosomalScaffolds <- function(clustered.gr=NULL, assembly.fasta=NULL, outputfolder=NULL) {
  ## Check if submitted fasta file is indexed
  assembly.fasta.idx <- paste0(assembly.fasta, ".fai")
  if (!file.exists(assembly.fasta.idx)) {
    message("Fasta file is not indexed, indexing ...")
    fa.idx <- Rsamtools::indexFa(file = assembly.fasta)
  }
  ## Open FASTA file instance
  fa.file <- open(Rsamtools::FaFile(assembly.fasta))
  
  ## Go through all clusters and exports pseudo-chromosomal scaffolds
  clustered.grl <- split(clustered.gr, clustered.gr$ID)
  for (i in seq_along(clustered.grl)) {
    cluster.regions <- clustered.grl[[i]]
    cluster.ID <- unique(cluster.regions$ID)
    message("Exporting FASTA for: ", cluster.ID)
    ## Read in contigs for a given cluster
    cluster.seq <- Rsamtools::scanFa(file = fa.file, param = cluster.regions, as = "DNAStringSet")
    ## Reverse complement based on 'dir' field
    revcomp <- which(cluster.regions$dir == 'revcomp')
    cluster.seq[revcomp] <- Biostrings::reverseComplement(cluster.seq[revcomp])
    ## Concatenate all sequences into a single FASTA separated by 100 N's.
    delim <- paste(rep('N', 100), collapse = '')
    cluster.seq.collapsed <- Biostrings::DNAStringSet(paste(cluster.seq, collapse = delim))
    ## Create unique sequence name
    seq.name <- paste0(cluster.ID, '_', paste(names(cluster.seq), collapse = "."))
    names(cluster.seq.collapsed) <- seq.name
    ## Write final FASTA
    destination <- file.path(outputfolder, paste0(cluster.ID, '.fasta')) 
    Biostrings::writeXStringSet(x = cluster.seq.collapsed, filepath = destination, format = 'fasta')
  }
}
