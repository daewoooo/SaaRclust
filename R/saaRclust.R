#' EM function
#'
#' This function performs basic steps of EM algorithm.
#'
#' @param tab.l A \code{list} of PB alignmetns separated per cell.
#' @param theta.param A \code{list} of estimated cell types for each cluster and each cell.
#' @param pi.param A \code{vector} of estimated sizes of each cluster based on initial hard clustering.
#' @param num.iter Set number of iteration to EM algorithm.
#' @param raw.counts Counts of +/- reads for evrry PB read.
#' @return A \code{list} of various exported results [pVal, pVal.logL, log.l, theta.param, pi.param].
#' @author David Porubsky, Maryam Ghareghani
#' @export

saarclust <- function(tab.l, theta.l=NULL, pi.param=NULL, num.iter=100, raw.counts=NULL) {

  message("Soft clustering")
  message("Running EM") 
  #counts.l <- list()
  BN.probs.l <- list()
  log.like.l <- list()

  for (i in 1:num.iter) {
    ptm <- startTimedMessage("    Iteration num: ",i," ...")  
    #message("Iteration num: ", i)

    clusters.per.cell <- list()
    clust.gammas.colsums.l <- list()
    clust.gammas.rowsums.l <- list()
    clust.gammas.norm.l <- list()
    #loop to obtain initial prob measures
    for (j in 1:length(tab.l)) {
      lib.name <- names(tab.l[j])
      #message("\tWorking on ",lib.name)
      lib.aligns <- tab.l[[j]]
      aligns.per.read <- split(lib.aligns$strand, lib.aligns$PBreadNames)
  
      params <- theta.l[[j]]
  
      #count number of W and C readc per PB read (only in first iteration)
      #if (i == 1) {
      #  counts <- t(sapply(aligns.per.read, function(x) table(x))) #TODO explore perhaps more efficient way to count +/- reads
      #} else {
        #counts <- counts.l[[j]]
      #}
      counts <- raw.counts[[j]]  
  
      #calc BN probs
      if (i == 1) {
        BN.probs <- countProb(minusCounts = counts[,1], plusCounts = counts[,2])
      } else {
        BN.probs <- BN.probs.l[[j]]
      }  
  
      #calc probs given theta
      clusters <- list()
      for (i1 in 1:nrow(params)) {
        clusters[[i1]] <- t(t(BN.probs)*params[i1,])
      }
      #Store rowsums of unnormalized cell type probs for every PB read in every cluster
      clusters.per.cell[[j]] <- lapply(clusters, rowSums)
  
      #calculate gamma function
      clust.gammas <- gammaFunction(clust.prob=clusters, pi = pi.param)
  
      #update theta (normalize gammas to 1)
      clust.gammas.norm <- clust.gammas$gammas.colsums/rowSums(clust.gammas$gammas.colsums)
  
      #store raw SS reads counts per PB read per SS library
      #counts.l[[j]] <- counts
      clust.gammas.colsums.l[[j]] <- clust.gammas$gammas.colsums
      clust.gammas.rowsums.l[[j]] <- clust.gammas$gammas.rowsums
      clust.gammas.norm.l[[j]] <- clust.gammas.norm
      BN.probs.l[[j]] <- BN.probs
    }

    #update pi
    #take sum over all cell types
    pi.update <- rowSums(sapply(clust.gammas.colsums.l, rowSums))
    #pi normalize to 1
    pi.norm.update <- pi.update/sum(pi.update)
    #update previous pi
    pi.param <- pi.norm.update

    #update previous theta
    theta.l <- clust.gammas.norm.l

    #calc likelihood function for any number of clusters
    clust.prod <- list()
    soft.probs <- list()
    for (i in 1:length(pi.param)) {
      clust.cell <- lapply(clusters.per.cell, `[[`, i)
      clust.rowsums.prod <- Reduce("*", clust.cell)
      clust.soft.prob <- clust.rowsums.prod * pi.param[i]
      clust.prod[[i]] <- clust.rowsums.prod
      soft.probs[[i]] <- clust.soft.prob
    }
    cluts.tab.update <- do.call(cbind, clust.prod)
    log.like <- sum(log(rowSums(cluts.tab.update * pi.param))*-1)
    log.like.l[[length(log.like.l)+1]] <- log.like
    
  stopTimedMessage(ptm)
  }
  
  #perform soft clustering
  soft.probs.tab <- do.call(cbind, soft.probs)
  soft.probs.tab.norm <- soft.probs.tab/rowSums(soft.probs.tab)
  cluts.tab <- do.call(cbind, clust.prod)
  cluts.tab.logL <- cluts.tab/rowSums(cluts.tab)
  #cluts.tab <- Reduce("*", clust.gammas.rowsums.l)
  #cluts.tab.norm <- cluts.tab/rowSums(cluts.tab)
  
  return(list(soft.pVal=soft.probs.tab.norm, pVal.logL=cluts.tab.logL, log.l=unlist(log.like.l), theta.param=theta.l, pi.param=pi.param))
  message("DONE!!!")
}  


