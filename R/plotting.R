#' Plot data quality measures.
#'
#' Takes imported data table and plots some relevant data quality measures.
#'
#' @param summary.tab Imported data table.
#' @author David Porubsky
#' @export

plotQualMeasure <- function(summary.tab) {

  mapp.dist.tab <- summary.tab$mapp.stat.counts
  mapp.gaps.tab <- summary.tab$mapp.gaps.stat
  SScov.stat.m <- summary.tab$SScov.stat          
  ord <- order(match(rownames(SScov.stat.m), mapp.dist.tab$PBreadNames))
  SScov.stat.m <- SScov.stat.m[ord,]
  
  #plotting distribution of SSreads mapped to PB reads
  read.map.dist.mean <- round(mean(mapp.dist.tab$SSread.perPB))
  read.dist.plt <- ggplot(mapp.dist.tab) + geom_line(aes(x=c(1:nrow(mapp.dist.tab)),y=SSread.perPB)) + theme_bw() + geom_hline(yintercept = read.map.dist.mean, color="red") + geom_text(aes(nrow(mapp.dist.tab), read.map.dist.mean, label = read.map.dist.mean, vjust = -1), color="red") + xlab("Sorted PB reads") + ylab("# of mapped SS reads")
  read.dist.plt.log <- ggplot(mapp.dist.tab) + geom_line(aes(x=c(1:nrow(mapp.dist.tab)),y=SSread.perPB)) + theme_bw() + geom_hline(yintercept = read.map.dist.mean, color="red") + geom_text(aes(nrow(mapp.dist.tab), read.map.dist.mean, label = read.map.dist.mean, vjust = -1), color="red") + xlab("Sorted PB reads") + ylab("# of mapped SS reads (log10)") + scale_y_continuous(trans = "log10")
  
  #plotting distribution of SSlibs represented per PB read
  SSlib.perPB.mean <- round(mean(mapp.dist.tab$SSlib.perPB))
  SSlib.perPB.plt <- ggplot(mapp.dist.tab) + geom_line(aes(x=c(1:nrow(mapp.dist.tab)),y=SSlib.perPB)) + theme_bw() + geom_hline(yintercept = SSlib.perPB.mean, color="red") + geom_text(aes(nrow(mapp.dist.tab), SSlib.perPB.mean, label = SSlib.perPB.mean, vjust = -1), color="red") + xlab("PB reads sorted by the number\nof SSreads per PBread") + ylab("# of SS libs per PB read")
  SSlib.perPB.plt.log <- ggplot(mapp.dist.tab) + geom_line(aes(x=c(1:nrow(mapp.dist.tab)),y=SSlib.perPB)) + theme_bw() + geom_hline(yintercept = SSlib.perPB.mean, color="red") + geom_text(aes(nrow(mapp.dist.tab), SSlib.perPB.mean, label = SSlib.perPB.mean, vjust = -1), color="red") + xlab("PB reads sorted by the number\nof SSreads per PBread") + ylab("# of SS libs per PB read (log10)") + scale_y_continuous(trans = "log10")
  
  #plotting sum of gaps per PB read normalized by SSread counts mapped to a given PB read
  gaps.perPB.top1perc <- round(quantile(mapp.dist.tab$gaps.perPB.norm,prob=1-1/100)) #99th quantile
  gaps.perPB.plt <- ggplot(mapp.dist.tab) + geom_line(aes(x=c(1:nrow(mapp.dist.tab)),y=gaps.perPB.norm)) + theme_bw() + geom_hline(yintercept = gaps.perPB.top1perc, color="red") + geom_text(aes(nrow(mapp.dist.tab), gaps.perPB.top1perc, label = gaps.perPB.top1perc, vjust = -1), color="red") + xlab("PB reads sorted by the number\nof SSreads per PBread") + ylab("Normalized sum of gaps per PB read")
  
  #plotting counts of SS reads per lib per PB read
  mapp.SScov.stat.mean <- mean(SScov.stat.m[SScov.stat.m>0])
  SScov.stat.m[SScov.stat.m==0] <- NA
  SScov.stat.df <- data.frame(value=rowMeans(SScov.stat.m, na.rm = T))
  SScov.stat.plt <- ggplot(SScov.stat.df) + geom_line(aes(x=c(1:nrow(SScov.stat.df)),y=value)) + theme_bw() + geom_hline(yintercept = mapp.SScov.stat.mean, color="red") + geom_text(aes(nrow(SScov.stat.df), mapp.SScov.stat.mean, label = mapp.SScov.stat.mean, vjust = -1), color="red") + xlab("PB reads sorted by the mean number\nof SSreads per SSlib per PBread") + ylab("Mean counts of SSreads per SSlib per PB read")
  
  #plotting match with gaps per SS read flagged TRUE or FALSE based on mapping agreement
  gaps.perSS.mean <- round(mean(mapp.gaps.tab$matchWithgaps))
  accur.counts <- table(mapp.gaps.tab$mapp.accur)
  falses <- round((accur.counts[1]/sum(accur.counts))*100, digits = 2)
  trues <- round((accur.counts[2]/sum(accur.counts))*100, digits = 2)
  false.lab <- paste('MISS ', falses, '%', sep = "")
  true.lab <- paste('MATCH ', trues, '%', sep = "")
  gaps.perSS.plt <- ggplot(mapp.gaps.tab) + geom_linerange(aes(x=c(1:nrow(mapp.gaps.tab)),ymin=0, ymax=matchWithgaps, color=mapp.accur)) + theme_bw() + geom_hline(yintercept = gaps.perSS.mean, color="red") + geom_text(aes(nrow(mapp.gaps.tab), gaps.perSS.mean, label = gaps.perSS.mean, vjust = -1), color="red") + xlab("Sorted SS reads by the length\nof SS read alignment (bp)") + ylab("(bp) SS read alignment with gaps (log10)") + scale_y_continuous(trans = "log10") + scale_color_manual(labels = c(false.lab, true.lab), values = c("red", "green"))
 
  #plot PB read length distribution
  PBreadLenDist <- ggplot(summary.tab$PBreadLenDist) + geom_linerange(aes(x=midpoints, ymin=0, ymax=freq), size=3)
  
  suppressWarnings( plt <- plot_grid(read.dist.plt, read.dist.plt.log, SSlib.perPB.plt, SSlib.perPB.plt.log, SScov.stat.plt, gaps.perPB.plt, PBreadLenDist, ncol = 2) )
  return(plt)
}


#' Plot heatmap of responsibilities of each PB reads for each cluster as a probability value.
#'
#' @param pVal.df A \code{data.frame} of probability values.
#' @param colOrder A \code{vector} of indices representing of cluster order on heatmap.
#' @param num.clusters Number of cluster present in probability table.
#' @author David Porubsky
#' @export

plotHeatmap <- function(pVal.df=NULL, colOrder=NULL, num.clusters=NULL) {

  #order clusters based on most likely chromosome partners                               
  if (!is.null(colOrder)) {                               
    pVal.df[,1:length(colOrder)] <- pVal.df[,colOrder]
    colnames(pVal.df)[1:num.clusters] <- colOrder
  }
  
  if (!is.null(num.clusters)) {

    chr.ids <- names(sort(table(pVal.df$PBchrom), decreasing = T))
    chr.ids <- gsub('^chr', '', chr.ids)
    #chr.ids <- sort(as.numeric(chr.ids))
    chr.colors <- rep(c("gray48","gray72"), ceiling(length(chr.ids)/2))
    chr.colors <- chr.colors[1:length(chr.ids)]
    names(chr.colors) <- chr.ids
    
    pVal.df$PBchrom <- gsub('^chr', '', pVal.df$PBchrom)
    pVal.df$PBchrom <- factor(pVal.df$PBchrom, levels=chr.ids)
    pVal.df <- pVal.df[order(pVal.df$PBchrom),]
    
    #set unexpected directionality flags to 1
    pVal.df$PBflag[pVal.df$PBflag != 0 & pVal.df$PBflag != 16] <- 1
    
    ha1 = rowAnnotation(df = data.frame(chr = pVal.df$PBchrom), col = list(chr=chr.colors))
    ha2 = rowAnnotation(df = data.frame(PB.dir = pVal.df$PBflag), col = list(PB.dir = c('0'="chocolate1", '16'="chartreuse3", '1'="white")))
    hm <- Heatmap(pVal.df[,c(1:num.clusters)], name = "Probs", cluster_columns = F, cluster_rows = F, show_row_names = FALSE)
    hm + ha1 + ha2
  } else {
    message("num.clusters not specified!!!\n")
  }  
}  


#' Plot accuracy of clustring in respect to known location of each PacBio read.
#'
#' @inheritParams plotHeatmap
#' @param thresh A \code{vector} of probability thresholds to obtain accuracy of clustering. 
#' @author David Porubsky
#' @export

plotClustAccuracy <- function(pVal.df=NULL, num.clusters=NULL, thresh=c(0.5,0.6,0.7,0.8,0.9,0.99)) {
  pVals <- pVal.df[,c(1:num.clusters)]
  chr.rows <- pVal.df$PBchrom
  
  acc.l <- list()
  for (th in thresh) {
    max.pVal <- apply(pVals, 1, max)
    mask <- max.pVal >= th
    ord <- apply(pVals[mask,], 1, which.max)
    chr.clusts <- split(chr.rows[mask], ord)
    chr.clusts <- lapply(chr.clusts, unlist) #Check why is this problem???
    clust.acc <- getClusterAcc(chr.clusts)
    frac.corr <- sum(clust.acc$stat$trues)/(sum(clust.acc$stat$trues) + sum(clust.acc$stat$miss))
    eval.reads <- table(mask)
    eval.reads <- eval.reads[2]/sum(eval.reads)
    acc.l[[as.character(th)]] <- c(frac.corr, eval.reads)
  } 
  
  #get hard clust accuracy
  chr.clusts <- split(chr.rows, pVal.df$hardClust)
  chr.clusts <- lapply(chr.clusts, unlist) #Check why is this problem???
  clust.acc <- getClusterAcc(chr.clusts)
  frac.corr <- sum(clust.acc$stat$trues)/(sum(clust.acc$stat$trues) + sum(clust.acc$stat$miss))
  
  m <- do.call(rbind, acc.l)
  df <- data.frame(values=m[,1], eval=m[,2], thresh = factor(thresh, levels=rev(thresh)))
  HC <- data.frame(values=frac.corr, eval=1, thresh='HardClust')
  
  plt <- ggplot(df) + geom_point(aes(x=values, y=eval, color=thresh), size=5) + geom_linerange(aes(ymin=-Inf, x=values, ymax=eval, color=thresh)) + scale_y_continuous(limits = c(0,1)) + scale_color_manual(values = brewer.pal(n=7, name="Set1"), name="Prob threshold") + ylab("(%) evaluated PB reads") + xlab("(%) correctly assigned PB reads")
  plt <- plt + geom_point(data=HC, aes(x=values, y=eval, color=thresh), size=5, inherit.aes = F) + geom_linerange(data=HC, aes(ymin=-Inf, x=values, ymax=eval, color=thresh), inherit.aes = F)
  return(plt)
}


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
#' @param report Plot either 'clustering' or 'ordering' of the contigs. Default: 'clustering'
#' @return A \code{ggplot} object.
#' @importFrom RColorBrewer brewer.pal.info brewer.pal
#' @importFrom tidyr separate
#' @author David Porubsky
#' @export
#' 
plotClusteredContigs <- function(bedfile, min.mapq=10, bsgenome=NULL, blacklist=NULL, report='clustering') {
  ## Use standard chromosomes only
  chroms <- paste0('chr', c(1:22, 'X','Y'))
  ## Read-in mapped congtigs to the human reference genome
  data <- read.table(bedfile, stringsAsFactors = FALSE)
  colnames(data) <- c('seqnames', 'start', 'end', 'info', 'mapq', 'dir')
  ## Filter contigs by mapping quality
  if (min.mapq > 0) {
    data <- data[data$mapq >= min.mapq,]
  }
  plt.df <- tidyr::separate(data, col = info, sep = '_', into = c('contig', 'order','cluster.ID'))
  plt.df$seqnames <- factor(plt.df$seqnames, levels=chroms)
  plt.df$order <- as.numeric(plt.df$order)
  ## Keep only standard chromosomes
  plt.df <- plt.df[plt.df$seqnames %in% chroms,]
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
  } else {
    message("Please choose to report either 'clustering' or 'ordering'!!!")
  }  
  
  ## Plot blacklisted regions in white if defined
  if (!is.null(blacklist)) {
    blacklist.df <- as.data.frame(blacklist)
    blacklist.df <- blacklist.df[blacklist.df$seqnames %in% chroms,]
    plt <- plt + geom_rect(data=blacklist.df , aes(xmin=start, xmax=end, ymin=0, ymax=1), fill='white')
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
