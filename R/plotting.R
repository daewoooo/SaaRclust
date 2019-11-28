#' Plot theta estimates resulting from EM algorithm.
#'
#' @param theta.param A \code{list} of estimated cell types for each cluster and each cell.
#' @param title A \code{character} to use as a title of the plot.
#' @importFrom reshape2 melt
#' @author David Porubsky
#' @export

plotThetaEstimates <- function(theta.param=NULL, title=NULL) {
  
  ptm <- startTimedMessage("Plotting theta estimates")
  plt.data <- list()
  for (j in 1:length(theta.param)) {
    df <- as.data.frame(theta.param[[j]])
    df$clustID <- rownames(df)
    df.plt <- suppressMessages( reshape2::melt(df) )
    df.plt$cell <- j
    plt.data[[j]] <- df.plt
  }
  plt.data.df <- do.call(rbind, plt.data)

  my_theme <-  theme(panel.spacing = unit(0, "lines"), 
                   strip.text.y = element_text(angle = 0),
                   axis.title.y=element_blank(),
                   axis.text.y=element_blank(),
                   axis.ticks.y=element_blank())
  if (is.null(title)) {
    plt <- ggplot(plt.data.df , aes(x=clustID, y=value, fill=variable)) + geom_bar(stat='identity', width=1) + facet_grid(cell ~ .) + scale_fill_manual(values = c('prob.cc'="paleturquoise4", 'prob.mix'="olivedrab",'prob.ww'="sandybrown")) + my_theme
  } else {
    plt <- ggplot(plt.data.df , aes(x=clustID, y=value, fill=variable)) + geom_bar(stat='identity', width=1) + facet_grid(cell ~ .) + scale_fill_manual(values = c('prob.cc'="paleturquoise4", 'prob.mix'="olivedrab",'prob.ww'="sandybrown")) + ggtitle(title) + my_theme
  }
  stopTimedMessage(ptm)
  return(plt)
}


#' Plot distribution of short reads mapped on top of PB reads
#'
#' @param count.list A \code{list} of short read mappings per library.
#' @author David Porubsky
#' @export

plotReadMappingDist <- function(count.list=NULL) {
  
  SSperPB <- list()
  for (j in 1:length(count.list)) {

      lib.aligns <- count.list[[j]]
      counts <- table(lib.aligns$PBreadNames)
      SSperPB[[j]] <- counts
  }
  all.counts <- do.call(rbind, SSperPB)
  plt.df1 <- as.data.frame(table(all.counts))
  plt1 <- ggplot(plt.df1) + geom_bar(aes(x=as.numeric(all.counts), y=Freq), stat='identity', fill="red") + xlab("# of ShortReads per PBread per Library") + ylab("Frequency") + scale_x_continuous(breaks = as.numeric(plt.df1$all.counts), labels = plt.df1$all.counts)
  
  count.list.collapsed <- do.call(rbind, count.list)
  counts <- table(count.list.collapsed$PBreadNames)
  plt.df2 <- as.data.frame(table(counts))
  
  is.odd <- function(x) x %% 2 != 0
  breaks <- as.numeric(plt.df2$counts)[ is.odd(as.numeric(plt.df2$counts)) ]
  plt2 <- ggplot(plt.df2) + geom_bar(aes(x=as.numeric(counts), y=Freq), stat='identity', fill="red") + xlab("# of ShortReads per PBread") + ylab("Frequency") + scale_x_continuous(breaks = breaks, labels = breaks)
  
  plt <- plot_grid(plt1, plt2, nrow = 1, rel_widths = c(1,2))
  return(plt)
}


#' Plot coverage of short reads mapped on top of PB reads
#'
#' @param minimap.tab A \code{data.frame} of short read mappings per PacBio read in maf.
#' @author David Porubsky
#' @export

plotReadAlignments <- function(minimap.tab=NULL) {
  #Convert table of alignments into GRanges object and then split into GRangesList by StrandS library ID
  minimap.tab.gr <- GenomicRanges::GRanges(seqnames=minimap.tab$PBchrom, strand=minimap.tab$strand, ranges=IRanges(start=minimap.tab$TargetCoordStart, end=minimap.tab$TargetCoordend), PBreadLen=minimap.tab$PBreadLen, SSlibNames=minimap.tab$SSlibNames)
  minimap.tab.grl <- GenomicRanges::split(minimap.tab.gr, minimap.tab.gr$SSlibNames)
  
  #get the name of PB read
  readID <- as.character(unique(minimap.tab$PBreadNames))
  
  all.libs <- list()
  #probs.l <- list()
  for (i in 1:length(minimap.tab.grl)) {
    gr <- minimap.tab.grl[[i]]
    gr$level <- GenomicRanges::disjointBins(gr)
    gr$level[which(GenomicRanges::strand(gr) == '-')] <- gr$level[which(GenomicRanges::strand(gr) == '-')] * -1
    
    #Get probabilities for StrandS read distribution
    dirRead.counts <- table(GenomicRanges::strand(gr))
    probs <- countProb(minusCounts = dirRead.counts['-'], plusCounts = dirRead.counts["+"], alpha = 0.1)
    probs.norm <- probs/sum(probs) #normalize prob values to 1
    probs.string <- paste(probs.norm, collapse = ", ")
    gr$probs <- probs.string
    
    #probs.df <- data.frame(minus=dirRead.counts['-'], plus=dirRead.counts["+"], ww=probs[,1], cc=probs[,2] ,wc=probs[,3], max=which.max(probs))
    #probs.l[[i]] <- probs.df
    
    plt.df <- as.data.frame(gr)
    all.libs[[i]] <- plt.df
  }
  all.libs.df <- do.call(rbind, all.libs)
  #all.probs.df <- do.call(rbind, probs.l)
  
  readLen <- data.frame(start=0, end=unique(all.libs.df$PBreadLen))
  plt <- ggplot(all.libs.df) + geom_linerange(data=readLen, aes(x=0, ymin=start, ymax=end), color="black") + geom_linerange(aes(x=level, ymin=start, ymax=end, color=strand)) + coord_flip() + scale_color_manual(values = c("paleturquoise4","sandybrown")) + xlab("") + facet_grid(SSlibNames ~ ., scales = 'free') + geom_text(aes(x=Inf,y=0, vjust=1, hjust=0), label=all.libs.df$probs) + ggtitle(readID) + theme(strip.text.y = element_text(angle = 360))
  return(plt)
}


#' Plot contig strand states per cell 
#' 
#' This function takes \code{data.frame} of strand states per contig [rows] and per cell [columns]
#' and plots heatmap of order or unordered strand states.
#'
#' @param contig.states A \code{data.frame} of strand states per contig and per cell.
#' @param cluster.rows If set to \code{TRUE}, will order rows by hierarchical clustering.
#' @param cluster.cols If set to \code{TRUE}, will order columns by hierarchical clustering.
#' @param filt.cols If set to \code{TRUE}, will remove columns with the same strand-state across all contigs.
#' @author David Porubsky
#' @export
#' 
plotContigStrandStates <- function(contig.states = NULL, cluster.rows=FALSE, cluster.cols=FALSE, filt.cols=FALSE) {
  ## Make sure that submitted object is a data.frame
  if (class(contig.states) != 'data.frame') {
    contig.states <- as.data.frame(contig.states)
  }
  ## Remove columns that have the same strand state across all contigs ('uninformative cells')
  if (filt.cols) {
    mask <- apply(contig.states, 2, function(x) length(unique(x)) > 1)
    if (length(mask[mask == TRUE]) > 1) {
      contig.states <- contig.states[,mask]
    } else {
      message("Parameter 'filt.cols' would leave only one cell, skipping ...")
    }  
  }
  ## Order rows by hierarchical clustering
  plt.df <- contig.states
  if (cluster.rows) {
    contig.dist <- stats::dist(contig.states)
    hc.clust <- stats::hclust(contig.dist)
    contig.order <- hc.clust$order
  }
  ## Order rows by user defined order
  if (is.numeric(cluster.rows)) {
    contig.order <- cluster.rows
  }
  ## Order columns by hierarchical clustering
  if (cluster.cols) {
    cell.dist <- stats::dist(t(contig.states))
    hc.clust <- stats::hclust(cell.dist)
    cell.order <- hc.clust$order
  }
  ## Order columns by hierarchical clustering
  if (is.numeric(cluster.cols)) {
    contig.order <- cluster.cols
  }
  ## Prepare data for plotting
  if (cluster.rows && cluster.cols) {
    plt.df <- plt.df[,cell.order]
    plt.df$contig <- factor(rownames(plt.df), levels = rownames(plt.df)[contig.order])
  } else if (cluster.rows && !cluster.cols) {
    plt.df$contig <- factor(rownames(plt.df), levels = rownames(plt.df)[contig.order])
  } else if (cluster.cols && !cluster.rows) {
    plt.df <- plt.df[,cell.order]
    plt.df$contig <- factor(rownames(plt.df), levels = rownames(plt.df))
  } else {
    plt.df$contig <- factor(rownames(plt.df), levels = rownames(plt.df))
  }
  plt.df <- reshape2::melt(plt.df, id.vars = 'contig')
  ## Plot contigs
  plt <- ggplot2::ggplot(plt.df) + 
    geom_tile(aes(x=variable, y=contig, fill=factor(value))) +
    scale_fill_manual(values = brewer.pal(n=4, name = 'Set1'), name='States') +
    xlab("Cell number") +
    ylab("Contig ID") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
  ## Return final plot
  return(plt)
} 


#' Plot genome-wide positions of clustered contigs
#'
#' @param bedfile An aligned contigs to the reference sequence in bed format.
#' @param min.mapq Minimum mapping quality of a contig to the refrence sequence.
#' @param bsgenome A \code{\link{GBSgenome-class}} object to provide chromosome lengths for plotting.
#' @param blacklist A \code{\link{GRanges-class}} object of regions to be removed.
#' @param report Plot either 'clustering', 'ordering' or 'orienting' of the contigs. Default: 'clustering'.
#' @param min.contig.size Minimal contigs size to plot.
#' @param reverse.x Set to \code{TRUE} if x-axis should be horizontaly reversed.
#' @return A \code{ggplot} object.
#' @importFrom RColorBrewer brewer.pal.info brewer.pal
#' @importFrom tidyr separate
#' @author David Porubsky
#' @export
#' 
plotClusteredContigs <- function(bedfile, min.mapq=10, bsgenome=NULL, blacklist=NULL, report='clustering', min.contig.size=NULL, reverse.x=FALSE) {
  ## Use standard chromosomes only
  chroms <- paste0('chr', c(1:22, 'X','Y'))
  ## Read-in mapped congtigs to the human reference genome
  data <- read.table(bedfile, stringsAsFactors = FALSE)
  colnames(data) <- c('seqnames', 'start', 'end', 'info', 'mapq', 'dir')
  ## Filter contigs by mapping quality
  if (min.mapq > 0) {
    data <- data[data$mapq >= min.mapq,]
  }
  ## Filter by contig size
  if (!is.null(min.contig.size)) {
    if (min.contig.size > 0) {
      mask <- (data$end - data$start) >= min.contig.size
      data <- data[mask,]
    }  
  }
  
  plt.df <- tidyr::separate(data, col = info, sep = '_', into = c('contig', 'order','cluster.ID'))
  plt.df$seqnames <- factor(plt.df$seqnames, levels=chroms)
  plt.df$order <- as.numeric(plt.df$order)
  ## Keep only standard chromosomes
  plt.df <- plt.df[plt.df$seqnames %in% chroms,]
  
  #plt.df$cluster.ID <- factor(plt.df$cluster.ID, levels = mixedsort(unique(plt.df$cluster.ID)))
  
  ## Prepare ideogram plot
  seq.len <- GenomeInfoDb::seqlengths(bsgenome)[chroms]
  ideo.df <- data.frame(seqnames=names(seq.len), length=seq.len)
  ideo.df$seqnames <- factor(ideo.df$seqnames, levels=chroms)
  ## Set chromosome cluster colors
  n.colors <- length(unique(plt.df$cluster.ID))
  qual.col.pals <- RColorBrewer::brewer.pal.info[RColorBrewer::brewer.pal.info$category == 'qual',]
  col.vector <- unlist(mapply(RColorBrewer::brewer.pal, qual.col.pals$maxcolors, rownames(qual.col.pals)))
  col.vector <- sample(col.vector, n.colors)
  ## Set chromosome order colors
  plt.df$ord.color <- ""
  for (chr in unique(plt.df$seqnames)) {
    chr.idx <- which(plt.df$seqnames == chr)
    colors <- gray.colors(max(plt.df$order[chr.idx]))
    plt.df$ord.color[chr.idx] <- colors[plt.df$order[chr.idx]]
  }
  
  ## Plot ideogram
  if (report == 'clustering') {
    plt <- ggplot2::ggplot() + geom_rect(data = ideo.df, aes(xmin=0, xmax=length, ymin=0, ymax=1), fill="white", color="black") +
      facet_grid(seqnames ~ ., switch = 'y') +
      geom_rect(data=plt.df, aes(xmin=start, xmax=end, ymin=0, ymax=1, fill=cluster.ID)) +
      scale_x_continuous(expand = c(0,0)) +
      scale_fill_manual(values = col.vector) + 
      theme_void() +
      theme(axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank()) +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
      theme(strip.text.y = element_text(angle = 180))
  } else if (report == 'ordering') {
    plt <- ggplot2::ggplot() + geom_rect(data = ideo.df, aes(xmin=0, xmax=length, ymin=0, ymax=1), fill=NA, color="black") +
      facet_grid(seqnames ~ ., switch = 'y') +
      geom_rect(data=plt.df, aes(xmin=start, xmax=end, ymin=0, ymax=1), fill=plt.df$ord.color) +
      scale_x_continuous(expand = c(0,0)) +
      theme_void() +
      theme(axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank()) +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
      theme(strip.text.y = element_text(angle = 180))
  } else if (report == 'orienting') {
    plt <- ggplot2::ggplot() + geom_rect(data = ideo.df, aes(xmin=0, xmax=length, ymin=0, ymax=1), fill=NA, color="black") +
      facet_grid(seqnames ~ ., switch = 'y') +
      geom_rect(data=plt.df, aes(xmin=start, xmax=end, ymin=0, ymax=1, fill=dir)) +
      scale_fill_manual(values = c('chocolate1', 'cadetblue4')) +
      scale_x_continuous(expand = c(0,0)) +
      theme_void() +
      theme(axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank()) +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
      theme(strip.text.y = element_text(angle = 180))
  } else {
    message("Please choose to report either 'clustering' or 'ordering'!!!")
  }  
  
  ## Plot blacklisted regions in white if defined
  if (!is.null(blacklist)) {
    blacklist.df <- as.data.frame(blacklist)
    blacklist.df <- blacklist.df[blacklist.df$seqnames %in% chroms,]
    plt <- plt + geom_rect(data=blacklist.df , aes(xmin=start, xmax=end, ymin=0, ymax=1), fill='white')
  }
  
  ## Reverse x-axis
  if (reverse.x) {
    plt <- plt + scale_x_reverse()
  }
  
  ## Return final plot
  return(plt)
}

#' ...
#'
#' @param dist.matrix ...
#' @param col.low ...
#' @param col.high ...
#' @return A \code{ggplot} object.
#' @author David Porubsky
#' @export
#'
plotDistanceMatrix <- function(dist.matrix, col.low="chartreuse4", col.high="cadetblue1") {
  dist.matrix.long <- reshape2::melt(dist.matrix)
  plt <- ggplot2::ggplot(dist.matrix.long, aes(x = Var2, y = Var1)) + 
    geom_raster(aes(fill = value)) + 
    scale_fill_gradient(low = col.low, high = col.high) +
    coord_fixed() +
    xlab("") +
    ylab("") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
  ## Return a final plot
  return(plt)
}


#' Plot assembly statistics
#'
#' @param infile A file that contains assembled contigs in specific format.
#' @param format Use 'bam' for contigs aligned to the reference or 'fai' for fasta index file or 'GRanges' for \code{\link{GRanges-class}} object.
#' @param title Add title to the plot.
#' @return A \code{ggplot} object.
#' @importFrom Rsamtools scanBamHeader
#' @importFrom scales comma
#' @author David Porubsky
#' @export
#'
plotAssemblyStat <- function(infile=NULL, format='bam', title=NULL) {
  if (format == 'bam') {
    ## Get contigs/scaffolds names and sizes from BAM
    file.header <- Rsamtools::scanBamHeader(infile)[[1]]
    chrom.lengths <- file.header$targets
    plt.df <- data.frame(ctg.len = sort(chrom.lengths))
  } else if (format == 'fai') {
    ## Get contigs/scaffolds names and sizes from fasta index
    fai.tab <- read.table(infile)
    plt.df <- data.frame(ctg.len = sort(fai.tab$V2))
  } else if (format == 'GRanges') {
    plt.df <- data.frame(ctg.len = sort(width(infile)))
  } else {
    message("Unsupported format, please use 'bam', 'fai' or 'GRanges' !!!")
  }
  
  ## Produce summary plot
  len.sorted <- rev(sort(as.numeric(plt.df$ctg.len)))
  N50 <- len.sorted[cumsum(len.sorted) >= sum(len.sorted)*0.5][1]
  N90 <- len.sorted[cumsum(len.sorted) >= sum(len.sorted)*0.9][1]
  total.size <- sum(len.sorted)/1000000000
  total.size <- round(total.size, digits = 2)
  total.size <- paste0('Total size = ', total.size, 'Gb')
  total.contigs <- paste0('Total contigs = ', nrow(plt.df))
  
  plt.df$x <- 1:nrow(plt.df)
  plt <- ggplot2::ggplot() + geom_point(data = plt.df, aes(x=x, y=ctg.len)) +
    geom_hline(yintercept = 1000000, linetype='dashed', color='red') +
    geom_hline(yintercept = N50, color='chartreuse4') +
    geom_hline(yintercept = N90, color='darkgoldenrod3') +
    geom_text(aes(x=0, y=Inf, label=total.size), color='black', vjust=2, hjust=0.1) +
    geom_text(aes(x=0, y=Inf, label=total.contigs), color='black', vjust=4, hjust=0.1) +
    geom_text(aes(x=0, y=N50, label=paste0('N50 = ', N50, 'bp')), color='black', vjust=-0.5, hjust=0.1) +
    geom_text(aes(x=0, y=N90, label=paste0('N90 = ', N90, 'bp')), color='black', vjust=-0.5, hjust=0.1) +
    scale_y_continuous(trans = 'log10', labels = scales::comma) +
    xlab("Size ordered contigs") +
    ylab("Contig length (log10)") +
    theme_bw()
  ## Add title if defined
  if (!is.null(title) & is.character(title)) {
    plt <- plt + ggtitle(title)
  }
  return(plt)
}


#' Plot distribution of cluster assignment probabilities.
#'
#' @param em.prob A \code{matrix} of probability assignments per contig and per cluster.
#' @inheritParams counts2ranges
#' @return A \code{ggplot} object.
#' @author David Porubsky
#' @export
#'
plotEMprobs <- function(em.prob=NULL, prob.th=0) {
  ptm <- startTimedMessage("Plotting probability distribution")
  max.probs <- apply(em.prob, 1, function(x) x[which.max(x)])
  max.probs.df <- data.frame(values=max.probs)
  suppressWarnings(
    plt <- ggplot2::ggplot(max.probs.df, aes(max.probs.df$values)) +
      geom_histogram(bins = 50) +
      scale_y_continuous(trans = 'log10') +
      xlab("Distribution of cluster assignment probabilities") +
      ylab("Counts (log10)") +
      theme_bw()
  )  
  if (prob.th > 0) {
    plt <- plt + geom_vline(xintercept = prob.th, color="red")
  }
  stopTimedMessage(ptm)
  return(plt)
}


#' Plot distribution of cluster sizes.
#'
#' @param clustered.gr A \code{\link{GRanges-class}} object with a contig region and their cluster assignment in the 'ID' metacolumn.
#' @return A \code{ggplot} object.
#' @importFrom dplyr %>%
#' @author David Porubsky
#' @export
#'
plotClusteredContigSizes <- function(clustered.gr=NULL) {
  ptm <- startTimedMessage("Plotting cluster sizes")
  ## Prepare data for plotting
  clustered.df <- as.data.frame(clustered.gr)
  plt.df <- clustered.df %>% dplyr::group_by(ID, dir) %>% dplyr::summarise(length=sum(width)) %>%
    dplyr::mutate(total.len = sum(length)) %>% dplyr::arrange(desc(total.len))
  plt.df$ID <- factor(plt.df$ID, levels = unique(plt.df$ID))
  
  ## Get chromosome breaks and labels
  max.len <- signif(max(c(plt.df$total.len)), digits = 2)
  breaks <- seq(from = 0, to = max.len, length.out = 6)
  labels <- breaks / 1000000
  labels <- paste0(labels, 'Mb')
  
  ## Make plot
  plt <- ggplot2::ggplot(data=plt.df, aes(x=ID, y=length, fill=dir)) +
    geom_col() +
    scale_fill_manual(values = c('cadetblue4','darkgoldenrod3'), name="Direction") +
    scale_y_continuous(breaks = breaks, labels = labels, name="Cluster size (Mb)", expand = c(0,0)) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
    xlab("")
  stopTimedMessage(ptm)
  return(plt)
}  
