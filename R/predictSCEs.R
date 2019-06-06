predictSCEs <- function(contig.states=NULL, matrix.ord=NULL, filt.cols=TRUE) {
  ## Remove columns that have the same strand state across all contigs ('uninformative cells')
  if (filt.cols) {
    mask <- apply(contig.states, 2, function(x) length(unique(x)) > 1)
    if (length(mask[mask == TRUE]) > 1) {
      contig.states <- contig.states[,mask]
    } else {
      message("Parameter 'filt.cols' would leave only one cell, skipping ...")
    }  
  }
  if (!is.null(matrix.ord) & is.numeric(matrix.ord)) {
    contig.states <- contig.states[matrix.ord,]
  }
  ## Extract positions from the input matrix
  regions.df <- data.frame(regions=rownames(contig.states))
  regions.df <- tidyr::separate(data = regions.df, col = regions, into = c('chr','start','end'), sep = ":|-")
  regions.gr <- GenomicRanges::GRanges(seqnames=regions.df$chr, ranges=IRanges(start=as.numeric(regions.df$start), end=as.numeric(regions.df$end)))
  ## Rescale coinheritance matrix based on SCE information
  mcols(regions.gr) <- contig.states
  SCEs.grl <- GRangesList()
  for (i in 1:ncol(cluster.m.new)) {
    regions.cell.gr <- regions.gr[,i]
    rle.obj <- Rle(mcols(regions.cell.gr)[,1])
    split.v <- rep(1:length(runLength(rle.obj)), runLength(rle.obj))
    regions.cell.grl <- split(regions.cell.gr, split.v)
    for (j in 1:length(regions.cell.grl)) {
      region.sub <- regions.cell.grl[[j]]
      region.sub <- region.sub[length(region.sub)]
      SCE.gr <- GRanges(seqnames=seqnames(region.sub), ranges=IRanges(start=end(region.sub)+2, end=end(region.sub)+3))
      SCEs.grl[[length(SCEs.grl) + 1]] <- SCE.gr
    }
  }
  SCEs.predict <- unlist(SCEs.grl, use.names = FALSE)
  return(SCEs.predict)
}  
