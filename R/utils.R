#' Collapse consecutive bins with the same ID value
#'
#' Collapse consecutive bins with the same value defined in 'id.field'. 
#'
#' @param gr A \code{\link{GRanges-class}} object.
#' @param id.field A number of metadata column to use for region merging.
#' @param measure.field A field column that contains measured value to be sum up.
#' @importFrom S4Vectors runLength Rle
#' @return A \code{\link{GRanges-class}} object.
#' @author David Porubsky
#' @export
#' 
collapseBins <- function(gr, id.field=0, measure.field=NULL) {
  ## Include seqnames into the unique ID
  unique.ID <- paste0(seqnames(gr), '_', mcols(gr)[,id.field])
  ## Get continous runs of the same unique ID
  unique.ID.runLength <- S4Vectors::runLength(S4Vectors::Rle(unique.ID))
  ind.last <- cumsum(unique.ID.runLength) ##get indices of last range in a consecutive(RLE) run of the same value
  ind.first <- c(1, cumsum(unique.ID.runLength) + 1) ##get indices of first range in a consecutive(RLE) run of the same value
  ind.first <- ind.first[-length(ind.first)]  ##erase last index from first range indices 
  collapsed.gr <- GenomicRanges::GRanges(seqnames=GenomeInfoDb::seqnames(gr[ind.first]), 
                                         ranges=IRanges::IRanges(start=start(gr[ind.first]), end=end(gr[ind.last])), 
                                         mcols=GenomicRanges::mcols(gr[ind.first]))
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
#' @param WWtoCC.penalty Penalize cases where there is W and C state in a single cell
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
  cluster.gr$group.ID <- paste0('cluster', cluster.gr$clust.ID)
  ## Add cluster group for pair of clusters specified in cluster.pairs
  for (i in seq_along(cluster.groups)) {
    clust.group <- cluster.groups[[i]]
    ## Assign common cluster group for a pair of clusters
    if (any(clust.group %in% cluster.gr$clust.ID)) {
      cluster.gr.sub <- cluster.gr[cluster.gr$clust.ID %in% clust.group]
      cluster.gr.sub$group.ID <- paste0('cluster',i)
      cluster.gr[cluster.gr$clust.ID %in% clust.group] <- cluster.gr.sub
    }  
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
#' @examples 
#' ## An example string to convert into \code{\link{GRanges-class}} object.
#' genomic.location.str <- 'chr4:128840654-143665255'
#' gr <- string2GRanges(region.string = genomic.location.str)
#' 
string2GRanges <- function(region.string=NULL) {
  regions.df <- data.frame(regions=region.string)
  regions.df <- tidyr::separate(data = regions.df, col = 'regions', into = c('chr','posANDdir'), sep = ":")
  regions.df <- tidyr::separate(data = regions.df, col = 'posANDdir', into = c('start','end','dir'), sep = "-|_", fill = 'right')
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
#' @importFrom S4Vectors subjectHits
#' @author David Porubsky
#' @export
#' @examples
#'## Get an example file
#'ctg.file <- system.file("extdata/clustered_assembly", "ordered&oriented_5e+06bp_chunks.RData", package = "SaaRclust")
#'## Load ordered and oriented contigs
#'ctg.gr <- get(load(ctg.file))
#'## Expand gaps inside contigs
#'ctg.gr.expanded <- expandGaps(gr = ctg.gr)

expandGaps <- function(gr) {
  ## Check if seqlengths are defined
  if (any(is.na(GenomeInfoDb::seqlengths(gr)))) {
    stop("Undefined seqlengths for parameter gr ...")
  }
  
  ## Extend gaps between ranges
  new.gr <- gr
  gaps.gr <- GenomicRanges::gaps(GenomicRanges::sort(gr))
  gaps.gr <- gaps.gr[GenomicRanges::strand(gaps.gr) == '*']
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
    GenomicRanges::start(new.gr[S4Vectors::subjectHits(follow.hits)]) <- GenomicRanges::start(gaps.gr[follow.idx])
    return(new.gr)
  } else {
    return(gr)
  }
}


#' Report summary of misassembled contigs/scaffolds
#' 
#' This function takes a \code{\link{GRanges-class}} object that contains putatively
#' misassembled contigs/scaffolds and report if this misassembly is likely caused by
#' an misorient, chimerism or both.
#' 
#' @param gr A \code{\link{GRanges-class}} object.
#' @return A \code{data.frame} object.
#' @importFrom BiocGenerics as.data.frame
#' @author David Porubsky
#' @export
#' @examples 
#'## Get an example file
#'putative.errors.file <- system.file("extdata/data", "putativeAsmErrors_5e+06bp_dynamic.RData", package = "SaaRclust")
#'## Load putative errors 
#'putative.errors.gr <- get(load(putative.errors.file))
#'## Export putative errors 
#'putative.errors.report <- reportMisAsmCTGs(gr = putative.errors.gr)
#'
reportMisAsmCTGs <- function(gr) {
  grl <- GenomicRanges::split(gr, seqnames(gr))
  
  report.l <- list()
  for (i in seq_along(grl)) {
    gr.sub <- grl[[i]]
    
    ## Get number of misassembled bases
    max.ctg <- gr.sub[which.max(width(gr.sub))]
    misasm.bases <- sum(width(gr.sub[!gr.sub$collapse.ID %in% max.ctg$collapse.ID]))
    
    ## Split contig chunks per cluster ID
    gr.sub.perClust <- GenomicRanges::split(gr.sub, gr.sub$ID)
    
    ## Classify observed assembly errors as misorient and/or chimerism
    if (any(sapply(gr.sub.perClust, function(x) length(unique(x$dir)) > 1))) {
      is.misorient <- TRUE
    } else {
      is.misorient <- FALSE
    }
    if (length(unique(gr.sub$ID)) > 1) {
      is.chimerism <- TRUE
    } else {
      is.chimerism <- FALSE
    }
    
    if (is.misorient & !is.chimerism) {
      asm.error <- 'misorient'
    } else if (is.chimerism & !is.misorient) {
      asm.error <- 'chimerism'
    } else if (is.chimerism & is.misorient) {
      asm.error <- 'misorient&chimerism'
    } else {
      asm.error <- 'unknown'
    }
    
    ctg.gr <- range(gr.sub)
    ctg.gr$asm.error <- asm.error
    ctg.gr$misasm.bases <- misasm.bases
    
    df <- BiocGenerics::as.data.frame(ctg.gr)
    report.l[[i]] <- df
  }
  report.df <- do.call(rbind, report.l)
  return(report.df)
}


#' Report clustering accuracy of de novo assembled contigs
#' 
#' This function takes alignments of de novo assembled contigs to the referene genome in BED format and reports
#' accuracy of their assignment to unique clusters.
#' 
#' @param bed.data A \code{data.frame} object containing alignments of de novo assembly to the reference genome in BED format.
#' @param cluster.ID A column name the uniquely defines assignment of contigs to unique clusters.
#' @param contig.ID A column name the uniquely defines contigs aligned to the reference genome.
#' @return A \code{data.frame} object.
#' @author David Porubsky
#' @export
#' 
getClusteringAcc <- function(bed.data, cluster.ID='ID', contig.ID='ctg') {
  required.cols <- c('seqnames', 'start', 'end')
  ## Split data per chromosome/scaffold
  if (all(required.cols %in% colnames(bed.data))) {
    data.l <- split(bed.data, bed.data$seqnames)
  } else {
    stop("Reuired columns 'seqnames', 'start' and 'end' not defined in the input data!!!")
  }  
  ## Make sure that contig.ID is defined in columns
  if (!contig.ID %in% colnames(bed.data)) {
    stop("Defined 'contig.ID' not present in the input data!!!")
  }  
  ## Make sure that contig.ID is defined in columns
  if (!cluster.ID %in% colnames(bed.data)) {
    stop("Defined 'cluster.ID' not present in the input data!!!")
  }  
  
  results.l <- list()
  wrong.ctg.grl <- GRangesList()
  for (i in seq_along(data.l)) {
    data <- data.l[[i]]
    seqname <- unique(data$seqnames)
    if(nrow(data) > 0) {
      data.gr <- GenomicRanges::makeGRangesFromDataFrame(data, keep.extra.columns = TRUE)
      ## Get ID of the cluster ID covering the most of the chromosome
      width.per.cluster <- split(width(data.gr), mcols(data.gr)[,cluster.ID])
      width.per.cluster <- sapply(width.per.cluster, sum)
      max.cluster.ID <- names(which.max(width.per.cluster))
      
      correct.ctgs <- unique(data[,contig.ID][data[,cluster.ID] == max.cluster.ID])
      wrong.ctgs <- unique(data[,contig.ID][!data[,contig.ID] %in% correct.ctgs])
      correct.gr <- data.gr[mcols(data.gr)[,cluster.ID] == max.cluster.ID]
      wrong.gr <- data.gr[mcols(data.gr)[,contig.ID] %in% wrong.ctgs]
      correct.ctgs.size <- sum(width(reduce(correct.gr)))
      wrong.ctgs.size <- sum(width(reduce(wrong.gr)))
      df <- data.frame(seqnames=seqname,
                       correct.ctgs=length(correct.ctgs), 
                       wrong.ctgs=length(wrong.ctgs), 
                       correct.ctgs.size=as.numeric(correct.ctgs.size), 
                       wrong.ctgs.size=as.numeric(wrong.ctgs.size), 
                       row.names = NULL)
      results.l[[i]] <- df
      wrong.ctg.grl[[length(wrong.ctg.grl) + 1]] <- wrong.gr
    }  
  }
  results.df <- do.call(rbind, results.l)
  return(results.df)
}


#' Report orienting accuracy of de novo assembled contigs
#' 
#' This function takes alignments of de novo assembled contigs to the referene genome in BED format and reports
#' accuracy of orienting within each unique cluster.
#' 
#' @param dir.ID A column name the uniquely defines mapping directionality of each contigs to the reference genome.
#' @inheritParams getClusteringAcc
#' @return A \code{data.frame} object.
#' @author David Porubsky
#' @export
#' 
getOrientingAcc <- function(bed.data, cluster.ID='ID', dir.ID='dir') {
  required.cols <- c('seqnames', 'start', 'end')
  ## Split data per chromosome/scaffold
  if (all(required.cols %in% colnames(bed.data))) {
    data.l <- split(bed.data, bed.data$seqnames)
  } else {
    stop("Reuired columns 'seqnames', 'start' and 'end' not defined in the input data!!!")
  }  
  ## Make sure that dir.ID is defined in columns
  if (!dir.ID %in% colnames(bed.data)) {
    stop("Defined 'dir.ID' not present in the input data!!!")
  }  
  ## Make sure that contig.ID is defined in columns
  if (!cluster.ID %in% colnames(bed.data)) {
    stop("Defined 'cluster.ID' not present in the input data!!!")
  }  
  
  results.l <- list()
  for (i in seq_along(data.l)) {
    data <- data.l[[i]]
    gr <- GenomicRanges::makeGRangesFromDataFrame(data, keep.extra.columns = TRUE)
    if (length(gr) > 0) {
      seqname <- as.character(unique(seqnames(gr)))
      cluster.id <- names(which.max(BiocGenerics::table(mcols(gr)[,cluster.ID])))
      dir.len <- sapply(split(width(gr), mcols(gr)[,dir.ID]), sum)
      if (length(dir.len) == 2) {
        major.dir <- max(dir.len)
        minor.dir <- min(dir.len)
      } else {
        major.dir <- max(dir.len)
        minor.dir <- 0
      }
      df <- data.frame(seqnames=seqname,
                       cluster.ID=cluster.id,
                       major.dir=major.dir, 
                       minor.dir=minor.dir)
      results.l[[i]] <- df
    }  
  }  
  results.df <- do.call(rbind, results.l)
  return(results.df)
} 


#' Report ordering accuracy of de novo assembled contigs
#' 
#' This function takes alignments of de novo assembled contigs to the referene genome in BED format and reports
#' accuracy of ordering within each unique cluster.
#' 
#' @param order.ID A column name the uniquely defines order of each contigs within a cluster.
#' @inheritParams getClusteringAcc
#' @return A \code{data.frame} object.
#' @author David Porubsky
#' @export
#' 
getOrderingAcc <- function(bed.data, cluster.ID='ID', order.ID='order') {
  required.cols <- c('seqnames', 'start', 'end')
  ## Split data per chromosome/scaffold
  if (all(required.cols %in% colnames(bed.data))) {
    data.l <- split(bed.data, bed.data$seqnames)
  } else {
    stop("Reuired columns 'seqnames', 'start' and 'end' not defined in the input data!!!")
  }  
  ## Make sure that order.ID is defined in columns
  if (!order.ID %in% colnames(bed.data)) {
    stop("Defined 'order.ID' not present in the input data!!!")
  }  
  ## Make sure that contig.ID is defined in columns
  if (!cluster.ID %in% colnames(bed.data)) {
    stop("Defined 'cluster.ID' not present in the input data!!!")
  }  
  
  results.l <- list()
  for (i in seq_along(data.l)) {
    data <- data.l[[i]]
    gr <- GenomicRanges::makeGRangesFromDataFrame(data, keep.extra.columns = TRUE)
    mcols(gr)[,order.ID] <- as.numeric(mcols(gr)[,order.ID])
    
    if (length(gr) > 0) {
      seqname <- as.character(unique(seqnames(gr)))
      cluster.id <- names(which.max(BiocGenerics::table(mcols(gr)[,cluster.ID])))
      total.ctg <- length(gr) 
      
      ctg.ord <- mcols(gr)[,order.ID]
      err <- length(which(diff(ctg.ord) < 0))
      corel1 <- cor(ctg.ord, sort(ctg.ord), method = "pearson")
      corel2 <- cor(rev(ctg.ord), sort(ctg.ord), method = "pearson")
      corel <- pmax(corel1, corel2)
        
      if (is.na(corel) & length(unique(ctg.ord)) == 1) {
        corel <- 1
      }
      
      df <- data.frame(seqnames=seqname, 
                       cluster.ID=cluster.id, 
                       err.ords=err, 
                       total.ctg=total.ctg, 
                       ord.cor=corel)
      results.l[[i]] <- df
    }  
  }
  results.df <- do.call(rbind, results.l)
  return(results.df)
}


#' Subtract and label set of genomic ranges in other set of genomic ranges.
#' 
#' This function takes a \code{\link{GRanges-class}} object of genomic regions to be labeled by user defined label
#'  in other \code{\link{GRanges-class}} object.
#' 
#' @param gr A \code{\link{GRanges-class}} object to be labeled based on the overlap with 'label.gr' object.
#' @param label.gr A \code{\link{GRanges-class}} object of genomic regions to be labeled inside 'gr' object.
#' @param label.gr.ID A user defined label to be used to mark 'label.gr' regions inside 'gr'.
#' @return A \code{\link{GRanges-class}} object.
#' @author David Porubsky
#' @export
#' 
labelGenomicRegions <- function(gr=NULL, label.gr=NULL, label.gr.ID=NULL) {
  ## Helper function
  markPartial <- function(gr=NULL, label.gr=NULL, label.gr.ID=NULL) {
    label.gr.sub <- IRanges::subsetByOverlaps(label.gr, gr)
    gr.new <- GenomicRanges::disjoin(c(gr[,0], label.gr.sub))
    mcols(gr.new) <- mcols(gr)
    hits <- IRanges::findOverlaps(label.gr.sub, gr.new)
    mcols(gr.new[subjectHits(hits)])$ploidy <- '1n'
    return(gr.new)
  }
  ## Set ploidy for contigs completely contained in 'label.gr' ranges
  gr$ploidy <- '2n'
  hits <- suppressWarnings( IRanges::findOverlaps(label.gr, gr, type = 'equal') )
  if (length(hits) > 0) {
    mcols(gr[subjectHits(hits)])$ploidy <- '1n'
  }  
  ## Label regions overlapping only partialy with 'label.gr' ranges
  label.gr.partial <- label.gr[-queryHits(hits)]
  gr.partial <- suppressWarnings( IRanges::subsetByOverlaps(gr, label.gr.partial) )
  gr <- suppressWarnings( IRanges::subsetByOverlaps(gr, label.gr.partial, invert = TRUE) )
  gr.extras <- GenomicRanges::GRangesList()
  for (i in seq_along(gr.partial)) {
    gr.extras[[i]] <- markPartial(gr = gr.partial[i], label.gr = label.gr.partial, label.gr.ID = label.gr.ID)
  }
  gr.extras <- unique(unlist(gr.extras))
  gr <- sort(c(gr, gr.extras))  
  return(gr)
}
