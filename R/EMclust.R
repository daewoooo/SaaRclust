#' EM function
#'
#' This function performs basic steps of EM algorithm.
#'
#' @param counts.l A \code{list} of PB alignmetns separated per cell.
#' @param theta.param A \code{list} of estimated cell types for each cluster and each cell.
#' @param pi.param A \code{vector} of estimated sizes of each cluster based on initial hard clustering.
#' @param num.iter Set number of iteration to EM algorithm.
#' @inheritParams countProb
#' @return A \code{list} of various exported results [pVal, pVal.logL, log.l, theta.param, pi.param].
#' @author David Porubsky, Maryam Ghareghani
#' @export

#saarclust <- function(tab.l, theta.l=NULL, pi.param=NULL, num.iter=100, raw.counts=NULL) {
EMclust <- function(counts.l, theta.param=NULL, pi.param=NULL, num.iter=100, alpha=0.1, logL.th=1) {

  if (num.iter>1) {
    message("Running EM") 
  } else {
    message("Soft clustering")
  }  
  
  #counts.l <- list()
  BN.probs.l <- list()
  log.like.l <- list()

  for (i in 1:num.iter) {
    ptm <- startTimedMessage("    Iteration num: ",i," ")  
    #message("Iteration num: ", i)

    clusters.per.cell <- list()
    clust.gammas.colsums.l <- list()
    clust.gammas.rowsums.l <- list()
    clust.gammas.norm.l <- list()
    #loop to obtain initial prob measures
    for (j in 1:length(counts.l)) {
      #lib.name <- names(tab.l[j])
      #message("\tWorking on ",lib.name)
      #lib.aligns <- tab.l[[j]]
      #aligns.per.read <- split(lib.aligns$strand, lib.aligns$PBreadNames)
  
      params <- theta.param[[j]]
  
      #get counts of W and C reads per PB read
      counts <- counts.l[[j]]
  
      #calc BN probs
      if (i == 1) {
        BN.probs <- countProb(minusCounts = counts[,1], plusCounts = counts[,2], alpha=alpha)
      } else {
        BN.probs <- BN.probs.l[[j]]
      }  
  
      #calc probs given theta
      clusters <- list()
      for (i1 in 1:nrow(params)) {
        clusters[[i1]] <- t(t(BN.probs)*params[i1,])
      }
      #Store rowsums of unnormalized BN probs (cell type) for every PB read in every cluster
      clusters.per.cell[[j]] <- lapply(clusters, rowSums)
  
      #calculate gamma function
      cellNum <- length(tab.l)
      pi.param.scaled <- ( pi.param^(1/cellNum) ) #scale pi.param given the number of strandseq cells
      clusters.scaled <- Map("*", clusters, pi.param.scaled) #multiply BN probs(multiplied by theta) with scaled pi.param
      clust.gammas <- gammaFunction(clust.prob = clusters.scaled, pi = pi.param, cellNum = cellNum)
  
      #update theta
      #clust.gammas.norm <- clust.gammas$gammas.colsums/rowSums(clust.gammas$gammas.colsums)
      clust.gammas.norm <- clust.gammas/rowSums(clust.gammas) #normalize gammas to 1
      
      #counts.l[[j]] <- counts  #store raw SS reads counts per PB read per SS library
      #clust.gammas.colsums.l[[j]] <- clust.gammas$gammas.colsums
      clust.gammas.colsums.l[[j]] <- clust.gammas
      #clust.gammas.rowsums.l[[j]] <- clust.gammas$gammas.rowsums
      clust.gammas.norm.l[[j]] <- clust.gammas.norm
      BN.probs.l[[j]] <- BN.probs
    }

    #keep old pi and theta param
    pi.param.old <- pi.param
    theta.param.old <- theta.param
    
    #update pi param
    #take sum over all cell types
    pi.update <- rowSums(sapply(clust.gammas.colsums.l, rowSums))
    #pi normalize to 1
    pi.norm.update <- pi.update/sum(pi.update)
    #update
    pi.param <- pi.norm.update

    #update theta param
    theta.param <- clust.gammas.norm.l
    
    #Adjust theta according to contraints
    #theta.Expected <- nrow(theta.l[[1]]) * c(0.25,0.25,0.5)
    #theta.l <- thetaRescale(theta.l = theta.l, theta.Expected = theta.Expected)

    #preprocess data in order to perform soft clustering and calculate likelihood function
    clust.prod <- list()
    soft.probs <- list()
    for (k in 1:length(pi.param)) {
      clust.cell <- lapply(clusters.per.cell, `[[`, k)
      clust.rowsums.prod <- Reduce("*", clust.cell)
      clust.soft.prob <- clust.rowsums.prod * pi.param[k]
      clust.prod[[k]] <- clust.rowsums.prod
      soft.probs[[k]] <- clust.soft.prob
    }
    cluts.tab.update <- do.call(cbind, clust.prod)
    
    #calc likelihood function
    log.like <- sum(log(rowSums(cluts.tab.update * pi.param)))*(-1)
    cluts.tab.update <- cluts.tab.update * pi.param
    #log.like <- kahansum(log(apply(cluts.tab.update, 1, kahansum)))*(-1)
    
    #Check the difference between previous and last results of likelihood function
    if (length(log.like.l) > 0) {
      log.l.diff <- log.like.l[[length(log.like.l)]] - log.like
      message("[",log.l.diff,"]", appendLF = FALSE)
      if (!is.null(logL.th)) {
        if (log.l.diff <= logL.th) {
          stopTimedMessage(ptm)
          message("CONVERGENCE!!!")
          break
        }  
      }
    }  
    
    if (! is.infinite(log.like)){
      log.like.l[[length(log.like.l)+1]] <- log.like
    } else {message("Warning: infinite objective function")}
    
  stopTimedMessage(ptm)
  }
  
  
  #perform soft clustering
  #soft.probs.tab.norm <- NULL
  #if (num.iter==1) {
  soft.probs.tab <- do.call(cbind, soft.probs)
  soft.probs.tab.norm <- soft.probs.tab/rowSums(soft.probs.tab)
  rownames(soft.probs.tab.norm) <- rownames(counts.l[[1]])
    #cluts.tab <- do.call(cbind, clust.prod)
    #cluts.tab.logL <- cluts.tab/rowSums(cluts.tab)
    #cluts.tab <- Reduce("*", clust.gammas.rowsums.l)
    #cluts.tab.norm <- cluts.tab/rowSums(cluts.tab)
  #}
  
  message("DONE!!!")  
  return(list(soft.pVal=soft.probs.tab.norm, log.l=unlist(log.like.l), theta.param=theta.param, pi.param=pi.param))
}  

