#' Hard clustering using k-means
#'
#' This function expects output from custom minimap test dataset that contains original locations of mapped reads in the genome.
#'
#' @param tab.l A \code{list} of PB alignmetns separated per cell.
#' @inheritParams SaaRclust
#' @return A \code{list} of estimated theta values for every cluster and cell.
#' @author David Porubsky
#' @export

hardClust <- function(tab.l=NULL, num.clusters=NULL, alpha=0.1) {

  message("Hard clustering")
  ptm <- startTimedMessage("    Counting directional reads ...") 
  ratios.l <- list()
  counts.l <- list()
  for (j in 1:length(tab.l)) {
    lib.name <- names(tab.l[j])
    #message("\tWorking on ",lib.name)
    lib.aligns <- tab.l[[j]]
    
    #count directional reads per PB read (SLOWEST)
    #aligns.per.read <- split(lib.aligns$strand, lib.aligns$PBreadNames)
    #counts <- t(sapply(aligns.per.read, function(x) table(x)))

    #count directional reads per PB read (FASTEST)  
    counts <- data.table(lib.aligns)[,table(strand),by=PBreadNames] #even faster option TEST
    cov.PBreads <- counts$PBreadNames
    uncov.PBreads <- levels(cov.PBreads)[!levels(cov.PBreads) %in% cov.PBreads]
    counts <- rbind(matrix(counts$V1, ncol=2, byrow = T), matrix(rep(0, 2*length(uncov.PBreads)), ncol=2) )
    rownames(counts) <- c(as.character(unique(cov.PBreads)), uncov.PBreads)
    counts <- counts[order(match(rownames(counts),levels(lib.aligns$PBreadNames))),]
    
    #count directional reads per PB read (MEDIUM)
    #counts <- aggregate(strand ~ PBreadNames, data=lib.aligns, FUN=table)
    #cov.PBreads <- counts$PBreadNames
    #uncov.PBreads <- levels(cov.PBreads)[!levels(cov.PBreads) %in% cov.PBreads]
    #counts <- rbind(counts[,2], matrix(rep(0, 2*length(uncov.PBreads)), ncol=2) )
    #rownames(counts) <- c(as.character(cov.PBreads), uncov.PBreads)
    #counts <- counts[order(match(rownames(counts),levels(lib.aligns$PBreadNames))),]
  
    ratios <- (counts[,2]-counts[,1])/(counts[,2]+counts[,1]) #calculate ratio of WW reads
    ratios[is.nan(ratios)] <- 0
    ratios.l[[j]] <- ratios
    counts.l[[j]] <- counts
  }
  stopTimedMessage(ptm)

  ptm <- startTimedMessage("    Kmeans clustering for ",num.clusters," clusters ...") 
  ratios.m <- do.call(cbind, ratios.l)
  ratios.m[ratios.m<0] <- -1
  ratios.m[ratios.m>0] <- 1
  km <- suppressWarnings( kmeans(ratios.m, centers = num.clusters, nstart = 100) )
  ord <- km$cluster
  #ratios.m.ord <- ratios.m[order(ord),]
  stopTimedMessage(ptm)

  ptm <- startTimedMessage("    Estimate theta values ...") 
  theta.estim <- list()
  for (j in 1:length(counts.l)) {
    minus.c <- split(counts.l[[j]][,1], ord)
    plus.c <- split(counts.l[[j]][,2], ord)
    
    #minus.counts <- sapply(minus.c, sum)
    #plus.counts <- sapply(plus.c, sum)
    #probs <- countProb2(minusCounts = minus.counts, plusCounts = plus.counts)
  
    clust.prob <- mapply(function(X,Y) { countProb(X,Y) }, X=minus.c, Y=plus.c)
    clust.prob.norm <- lapply(clust.prob, function(x) colSums(log(x)))
    estimates <- sapply(clust.prob.norm, which.max)
    
    #Assign cell type probs based on the majority type in each cluster
    probs <- list()
    for (i in 1:length(clust.prob)) {
      estim <- estimates[i]
      if (estim == 1) {
        theta <- c(1-alpha, alpha/2, alpha/2)
      } else if (estim == 2) {
        theta <- c(alpha/2, 1-alpha, alpha/2)
      } else {
        theta <- c(alpha/2, alpha/2, 1-alpha)
      }
      probs[[i]] <- theta
    }
    probs <- do.call(rbind, probs)
    theta.estim[[j]] <- probs
  }
  stopTimedMessage(ptm)
  
  #return(list(theta.estim=theta.estim, clust.id=ord, raw.counts=counts.l))
  return(list(theta.estim=theta.estim, clust.id=ord))
}  
