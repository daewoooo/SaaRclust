#' Obtain summary measures of minimap output.
#'
#' Takes imported data table and export relevant data quality measures.
#'
#' @param data.tab Imported data table.
#' @author David Porubsky
#' @export

getQualMeasure <- function(data.tab) {
  #get number of SS reads per PB read
  SSread.perPB <- sort(table(data.tab$PBreadNames), decreasing = T)
  summary.df <- data.frame(PBreadNames = names(SSread.perPB), SSread.perPB = as.vector(SSread.perPB), stringsAsFactors = F)
  
  #get number of SS libs per PB read
  SSlib.perPB <- split(as.factor(data.tab$SSlibNames), data.tab$PBreadNames)
  SSlib.perPB.counts <- lapply(SSlib.perPB, function(x) table(x))
  SSlib.perPB.m <- do.call(rbind,SSlib.perPB.counts)
  
  SSlib.perPB.counts <- sapply(SSlib.perPB, function(x) length(unique(x)))
  summary.df$SSlib.perPB <- SSlib.perPB.counts[summary.df$PBreadNames]
  
  #get logical vector comapring mapping accuracy
  mapp.accur.comp <- data.tab$SSchrom == data.tab$PBchrom
  mapp.gaps <- data.tab$MatchedBasesWithGaps 
  mapp.gaps.df <- data.frame(matchWithgaps=mapp.gaps, mapp.accur=mapp.accur.comp)
  mapp.gaps.df <- mapp.gaps.df[order(mapp.gaps.df$matchWithgaps, decreasing = T),]
  mapp.accur <- table(mapp.accur.comp)
  
  #get number of gaps in SS read alignments
  match.gaps <- split(data.tab$MatchedBasesWithGaps, data.tab$PBreadNames)
  match.gaps.sums <- sapply(match.gaps, sum)
  summary.df$gaps.perPB <- match.gaps.sums[summary.df$PBreadNames]
  summary.df$gaps.perPB.norm <- summary.df$gaps.perPB / as.vector(SSread.perPB)
  
  #get PB read distribution hist
  hist.data <- hist(data.tab$PBreadLen, breaks = 100)
  hist.df <- data.frame(midpoints= hist.data$mids, freq= hist.data$counts)
  
  return(list(mapp.stat.counts = summary.df, mapp.gaps.stat=mapp.gaps.df, mapp.accur=mapp.accur, SScov.stat=SSlib.perPB.m, PBreadLenDist=hist.df))
}

#' Check clustering accuracy
#'
#' Takes assigned identity of each PB read and compares it to known location.
#'
#' @param clusters List of clusters with predicted chromosome identity for each PB read
#' @author David Porubsky
#' @export

getClusterAcc <- function(clusters) {
  stat <- list()
  miss.reads <- list()
  for (i in 1:length(clusters)) { 
    cluster <- clusters[[i]]
    max.chr <- names(which.max(table(cluster)))
    trues <- length(cluster[cluster == max.chr])
    miss <- length(cluster[cluster != max.chr])
    stat.df <- data.frame(trues=trues, miss=miss, id=i)
    
    miss.PBreads <- names(cluster[cluster != max.chr])
    miss.chr <- unname(cluster[cluster != max.chr])
    if (length(miss.PBreads)>0) {
      tab <- data.frame(miss.PBreads=miss.PBreads, miss.chr=miss.chr, id=i)
      miss.reads[[length(miss.reads)+1]] <- tab
    }
    stat[[length(stat)+1]] <- stat.df  
  }
  stat <- do.call(rbind, stat)
  miss.reads <- do.call(rbind, miss.reads)
  return(list(stat=stat, miss.reads=miss.reads))
}

                               
#' Computes the performance of the hard clustering algorithm
#' @param hard.clust hard clustring result
#' @param pb.chr The true chromosomes of PB reads
#' @param pb.flag The flags of PB reads
#' @param tab.filt table of filtered read counts
#' @param female a binary argument showing the sex of the individual
#' @author Maryam Ghareghani
#' @export

hardClustAccuracy <- function(hard.clust, pb.chr, pb.flag, tab.filt, female=TRUE)
{
  # filter and keep PB reads that have only defined (true) chromosome names with flag 0 or 16
  if(female)
  {
    filt.chrom <- which(grepl('^chr[0-9X][0-9]?$',tab.filt$PBchrom))
  } else {
    filt.chrom <- which(grepl('^chr[0-9XY][0-9]?$',tab.filt$PBchrom))
  }
  filt.flag <- union(which(tab.filt$PBflag == 0), which(tab.filt$PBflag == 16))
  filt <- intersect(filt.chrom, filt.flag)
  tab <- tab.filt[filt,]
  
  # define the classe (true clusters)
  classes <- paste0(pb.chr, "_", pb.flag)
  # some extra filtering ---> take only those PB reads that have flag 0 or 16
  names(classes) <- names(pb.chr)
  # filter out the PB reads that have more than one chr or flag
  oneChrPBreads <- which(sapply(pb.chr, length)==1)
  oneFlagPBreads <- which(sapply(pb.flag, length)==1)
  classes <- classes[intersect(oneChrPBreads, oneFlagPBreads)]
  
  classes <- classes[which(names(classes) %in% tab$PBreadNames)]
  clusters <- hard.clust[which(names(hard.clust) %in% names(classes))]
  ord <- order(names(classes))
  classes <- classes[ord]
  clusters <- clusters[ord]
  
  t <- table(classes, clusters)
  clust.assigned <- apply(t, 2, function(x) which(x == max(x))[1])
  clust.purity <- sum(sapply(1:length(clust.assigned), function(i) t[clust.assigned[i],i]))/length(clusters)
  
  missed <- setdiff(1:nrow(t), clust.assigned)
  names(missed) <- rownames(t)[missed]
  clust.purity.perClust <- sapply(1:length(clust.assigned), function(i) t[clust.assigned[i],i]/sum(t[,i]))
  
  df <- data.frame(clustIDs=clust.assigned, clustNames=rownames(t)[clust.assigned], accuracy=clust.purity.perClust)
  df <- df[order(clust.purity.perClust),]
  
  
  list(acc=clust.purity, missed.clusters = missed)
}

#' Rescale theta values for WC cell type
#'
#' Set WC probs based on WW and CC probs distribution (if ratio of WW and CC == 1 then WC region with high prob)
#'
#' @param theta.l List of estmated theta values for every single cell.
#' @author David Porubsky
#' @export

rescaleTheta <- function(theta.l) {

  new.theta <- list()
  for (i in 1:length(theta.l)) {
    theta.cell <- theta.l[[i]]
    ratios <- theta.cell[,1] / theta.cell[,2]
    mask <- ratios>0.8 & ratios<1.2
    theta.cell[mask,1] <- 0.05
    theta.cell[mask,2] <- 0.05
    theta.cell[mask,3] <- 0.9
    new.theta[[i]] <- theta.cell
  }
  return(new.theta)
}  



#' Simulate random theta estimates for cell type
#'
#' Function to simulate random theta values in order to initialize EM algorithm.
#'
#' @param num.cells Number of cells used in the analysis.
#' @param num.clusters Expected number of clusters. (for 22 autosomes == 44 clusters)
#' @author David Porubsky
#' @export

randomTheta <- function(num.cells=100, num.clusters=44) {
  theta.l <- list()
  for (j in 1:num.cells) {
    cell.theta <- list()
    for (k in 1:num.clusters) {
      random.type <- runif(3)
      max.type <- which.max(random.type)
      if (max.type == 1) {
        theta <- c(0.9, 0.05, 0.05)
      } else if (max.type == 2) {
        theta <- c(0.05, 0.9, 0.05)
      } else {
        theta <- c(0.05, 0.05, 0.9)
      }
      cell.theta[[k]] <- theta
    }
    theta.l[[j]] <- do.call(rbind, cell.theta)
  }
  return(theta.l)
}


#' Rescale theta values for WC cell type.
#'
#' Adjust estimated theta values based on expected number of WC chromosomes given random segregation of parental homologues.
#'
#' @param theta.param A \code{list} of estimated cell types for each cluster and each cell.
#' @param theta.Expected A \code{vector} of expected sums of probabilities to observe WW,CC or WC chromosomes across all cells.
#' @author David Porubsky
#' @export

thetaRescale <- function(theta.param, theta.Expected) {
  theta.rescaled <- list()
  for (i in 1:length(theta.param)) {
    theta.cell <- theta.param[[i]]
    theta.ObservedSums <- colSums(theta.cell)
    corr.factor <- theta.Expected / theta.ObservedSums
    theta.param.rescaled <- t(t(theta.cell) * corr.factor)
    theta.param.rescaled.norm <- theta.param.rescaled / rowSums(theta.param.rescaled)
    theta.rescaled[[i]] <- theta.param.rescaled.norm
  }
  return(theta.rescaled)
}

#' Alternative sum function to reduce numerical error caused by finite precision of floating point numbers.
#'
#' @param x A \code{vector} of values tp sum up.
#' @export

kahansum <- function(x) {
  ks <- 0
  c <- 0
  for(i in 1:length(x)) {
    y <- x[i] - c
    kt <- ks + y
    c = (kt - ks) - y
    ks = kt
  }
  ks
}





exportGenomicLocations <- function(soft.clust, prob.th=0.6) {
  
  getMaxPvals <- function(pval.vector) {
    srt.pval <- sort(pval.vector, decreasing = T)
    breakpoint <- which.min(diff(srt.pval))
    max.pvals <- srt.pval[1:breakpoint]
    max.pvals.clust <- which(pval.vector %in% max.pvals)
    return(max.pvals.clust)
  }  
  
  max.pVal <- apply(soft.clust, 1, max)
  mask <- max.pVal >= prob.th
  ord <- apply(soft.clust[mask,], 1, which.max)
  
  read.IDs <- rep(list(0), nrow(soft.clust)) 
  read.IDs[mask] <- ord 
  
  best.idx <- apply(soft.clust[!mask,], 1, getMaxPvals)
  read.IDs[!mask] <- best.idx 
  
  names(read.IDs) <- rownames(soft.clust)

  return(read.IDs)
  #best.idx <- t( apply(soft.clust[!mask,], 1, function(x) which(x %in% sort(x, decreasing = T)[1:5])) )
  #best.pval <- t( apply(soft.clust[!mask,], 1, function(x) x[which(x %in% sort(x, decreasing = T)[1:5])]) )
}





getClusterIdentity <- function(soft.clust, chr.rows) {
  max.Clust <- apply(soft.clust, 1, which.max)
  clustByChrom <- split(max.Clust, chr.rows)
  clustIdPerChrom <- lapply(clustByChrom, function(x) as.numeric(names(sort(table(x), decreasing=T)[1:2])))
  return(clustIdPerChrom)
}

