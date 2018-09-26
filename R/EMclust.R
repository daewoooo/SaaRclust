#' EM (Expectation maximization) function
#'
#' This function performs basic steps of EM algorithm.
#'
#' @param counts.l A \code{list} of PB alignmetns separated per cell.
#' @param theta.param A \code{list} of estimated cell types for each single cell. (rows=clusters, cols=strand states)
#' @param pi.param A \code{vector} of estimated sizes of each cluster based on initial hard clustering.
#' @param num.iter Set number of iteration to EM algorithm.
#' @param logL.th Set the difference between objective function from the current and the previous interation for EM algorithm to converge.
#' @inheritParams countProb
#' @return A \code{list} of various exported results [pVal, pVal.logL, log.l, theta.param, pi.param].
#' @importFrom matrixStats logSumExp rowLogSumExps
#' @author David Porubsky, Maryam Ghareghani
#' @export

#saarclust <- function(tab.l, theta.l=NULL, pi.param=NULL, num.iter=100, raw.counts=NULL) {
EMclust <- function(counts.l, theta.param=NULL, pi.param=NULL, num.iter=100, alpha=0.1, logL.th=1, log.scale=FALSE) {

  if (num.iter>1) {
    message("Running EM") 
  } else {
    message("Rescaling theta parameter")
  }
  
  ## Log scale theta and pi parameters if log.scale=TRUE
  if (log.scale) {
    theta.param <- lapply(theta.param, log)
    pi.param <- log(pi.param)
  }
  
  #counts.l <- list()
  BN.probs.l <- list()
  log.like.l <- list()

  ## Run user defined number of EM iterations ##
  for (i in 1:num.iter) {
    ptm <- startTimedMessage("    Iteration num: ",i," ")  

    clusters.per.cell <- list()
    clust.gammas.colsums.l <- list() ## list of matrices per single cell: matrix rows=clusters, cols=strand states
    #clust.gammas.rowsums.l <- list()
    theta.update.l <- list() ## list of matrices per single cell: matrix rows=clusters, cols=strand states
    num.removed.reads <- 0
    ## Loop over all Strand-seq cells
    for (j in 1:length(counts.l)) {
      #lib.name <- names(tab.l[j])
      #message("\tWorking on ",lib.name)
      #lib.aligns <- tab.l[[j]]
      #aligns.per.read <- split(lib.aligns$strand, lib.aligns$PBreadNames)
      
      ## Get theta estimates for a given cell
      #TODO: rename params to cell.params
      params <- theta.param[[j]]
    
      ## Get counts of W and C reads per read/genomic segment for a given cell
      counts <- counts.l[[j]]
  
      ## Calculate BN probabilities
      if (i == 1) {
        BN.probs <- countProb(minusCounts = counts[,1], plusCounts = counts[,2], alpha=alpha, log=log.scale)
        BN.probs.l[[j]] <- BN.probs
      } else {
        BN.probs <- BN.probs.l[[j]]
      }  
  
      ## Remove BN probs and corresponding PB reads which probability is assigned NaN value or when all probabilities are zero [TODO!!! resolve this problem]
      if (log.scale) {
        mask <- which( apply(BN.probs, 1, function(x) any(is.nan(x) | logSumExp(x)==-Inf)) )
      } else {
        mask <- which( apply(BN.probs, 1, function(x) any(is.nan(x) | sum(x)==0)) )
      }
        
      if (length(mask)>0) { #remove PB reads with extremely high StrandS read counts
        message(length(mask), " NAN numbers after binomial probability calculation!!!")
        num.removed.reads <- num.removed.reads + length(mask)
        BN.probs <- BN.probs[-mask,]
        counts.l[[j]] <- counts[-mask,]
      }
        
      ## Multiplication of BN.probs by theta parameter for a given cell ##
      clusters <- list()
      ## clusters: represent list of matrices per cluster. Matrix: rows=reads/genomic segments, cols=strand states
      #TODO: rename clusters to ...
      for (i1 in 1:nrow(params)) {
        ## In the log scale change multiplication changes to summation
        if (log.scale) {
          clusters[[i1]] <- t(t(BN.probs)+params[i1,])
        } else {
          clusters[[i1]] <- t(t(BN.probs)*params[i1,])
        }
      }
      
      ## Take sums over strand states (cols) for every read/genomic segment (rows) of each cluster ##
      ## clusters.per.cell[[[j]]: a list of vectors per cluster. Each vector represent rowSums of each cluster (matrix)
      if (log.scale) {
        #clusters.per.cell[[j]] <- lapply(clusters, function(x) apply(x, 1, logSumExp))
        clusters.per.cell[[j]] <- lapply(clusters, rowLogSumExps)
      } else {
        clusters.per.cell[[j]] <- lapply(clusters, rowSums)
      }  
      
      ## Scale pi.param given the number of Strand-seq cells ##
      cellNum <- length(counts.l)
      if (log.scale) {
        pi.param.scaled <- ( pi.param*(1/cellNum) )
      } else {
        pi.param.scaled <- ( pi.param^(1/cellNum) )
      }  
      
      ## Multiply BN probs (multiplied previously by theta) with scaled pi.param ##
      ## clusters.scaled: list of matrices per cluster for a given cell
      ## Map multiplies each element of the list with correspoding element of the vector
      if (log.scale) {
        clusters.scaled <- Map("+", clusters, pi.param.scaled)
      } else {
        clusters.scaled <- Map("*", clusters, pi.param.scaled)
      }
      
      ## Calculate gamma function ##
      clust.gammas.l <- gammaFunction(clust.prob=clusters.scaled, pi.scaled=pi.param.scaled, cellNum=cellNum, log.scale=log.scale)
  
      ## Take sum over all reads/genomic segments ##
      if (log.scale) { 
        #clust.gammas.colsums <- lapply(clust.gammas.l, function(x) apply(x, 2, logSumExp))
        clust.gammas.colsums <- lapply( clust.gammas.l, function(x) rowLogSumExps(t(x)) )
      } else {  
        clust.gammas.colsums <- lapply(clust.gammas.l, colSums) ## clust.gammas.colsums: matrix with rows=clusters, cols=strand states
      }
      clust.gammas.colsums <- do.call(rbind, clust.gammas.colsums)
      
      ## Updating theta parameters
      if (log.scale) {
        #theta.update <- clust.gammas.colsums - apply(clust.gammas.colsums, 1, logSumExp)
        theta.update <- clust.gammas.colsums - rowLogSumExps(clust.gammas.colsums)
      } else {
        theta.update <- clust.gammas.colsums / rowSums(clust.gammas.colsums)
      }  
      
      #counts.l[[j]] <- counts  #store raw SS reads counts per PB read per SS library
      #clust.gammas.colsums.l[[j]] <- clust.gammas$gammas.colsums
      clust.gammas.colsums.l[[j]] <- clust.gammas.colsums
      #clust.gammas.rowsums.l[[j]] <- clust.gammas$gammas.rowsums
      theta.update.l[[j]] <- theta.update
    } #end of the loop over all Strand-seq cells
    
    ## Report number of removed reads
    if (num.removed.reads > 0) {
      message("Warning: Removing ", num.removed.reads, " PacBio reads with all zero or NaN binom probability!!!", appendLF = FALSE)
    }  
    
    ## Keep old pi and theta parameters
    #pi.param.old <- pi.param
    #theta.param.old <- theta.param
    
    ## Update pi parameter ##
    ## Take sum over all cell types
    if (log.scale) {
      #strand.states.sums <- sapply(clust.gammas.colsums.l, function(x) apply(x, 1, logSumExp)) #sums over all strand states per single cell. Results in matrix with rows=clusters, cols=cells
      #pi.update <- apply(strand.states.sums, 1, logSumExp)
      pi.update <- rowLogSumExps(sapply(clust.gammas.colsums.l, rowLogSumExps)) #sums over all strand states per single cell. Results in matrix with rows=clusters, cols=cells
      
      ## Normalize pi parameter to 1 and update pi
      pi.norm.update <- pi.update - logSumExp(pi.update)
    } else {
      ## first rowSums = sums over all strand states per single cell. Results in matrix with rows=clusters, cols=cells
      ## second rowSums = sums over single cells
      ## pi.update: vector of lentgh of number of clusters. Nominator of pi updating formula
      pi.update <- rowSums(sapply(clust.gammas.colsums.l, rowSums))
      
      ## Normalize pi parameter to 1 and update pi
      pi.norm.update <- pi.update / sum(pi.update)
    }  
    pi.param <- pi.norm.update
    
    ## Assign upated theta param as a current theta parameter
    theta.param <- theta.update.l
    
    ## Adjust theta according to contraints
    #theta.Expected <- nrow(theta.l[[1]]) * c(0.25,0.25,0.5)
    #theta.l <- thetaRescale(theta.l = theta.l, theta.Expected = theta.Expected)

    ## Preprocess data in order to perform soft clustering and calculate likelihood function ##
    #clust.prod <- list()
    soft.probs <- list()
    for (k in 1:length(pi.param)) {
      clust.cell <- lapply(clusters.per.cell, `[[`, k)
      if (log.scale) {
        ## accessing cluster (k) in clusters.per.cell object 
        clust.rowsums.prod <- Reduce("+", clust.cell)
        clust.soft.prob <- clust.rowsums.prod + pi.param[k] #clust.soft.prob: vector of cluster assignment likelihoods per reads/genomic segments
      } else {
        ## accessing cluster (k) in clusters.per.cell object 
        clust.rowsums.prod <- Reduce("*", clust.cell)
        clust.soft.prob <- clust.rowsums.prod * pi.param[k] #clust.soft.prob: vector of cluster assignment likelihoods per reads/genomic segments
      }
      #clust.prod[[k]] <- clust.rowsums.prod
      soft.probs[[k]] <- clust.soft.prob
    }
    cluts.tab.update <- do.call(cbind, soft.probs) #rows=reads/genomic segments, cols=clusters
    
    ## Calculate the likelihood function ##
    if (log.scale) {
      cluts.tab.sums <- apply(cluts.tab.update, 1, logSumExp)
      log.like <- sum(cluts.tab.sums)*(-1)
    } else {
      cluts.tab.sums <- rowSums(cluts.tab.update)
      log.like <- sum(log(cluts.tab.sums))*(-1)
    }  
    #cluts.tab.sums[cluts.tab.sums == 0] <- min(cluts.tab.sums[cluts.tab.sums != 0]) #set zero sums to the lowest sums!!! (to approximate likelihood function)
    #log.like <- sum(log(cluts.tab.sums))*(-1)
    #cluts.tab.update <- cluts.tab.update * pi.param
    #log.like <- kahansum(log(apply(cluts.tab.update, 1, kahansum)))*(-1)
    
    ## Check the difference between previous and last results of likelihood function
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
    
    if (!is.infinite(log.like)) {
      log.like.l[[length(log.like.l)+1]] <- log.like
    } else {
      message("Message: LogL function numerical problem", appendLF = FALSE)
    }
    
  stopTimedMessage(ptm)
  }
  
  
  ## Perform soft clustering
  #soft.probs.tab.norm <- NULL
  #if (num.iter==1) {
  ## Scaling soft probabilities to 1
  if (log.scale) {
    #soft.probs.tab.norm <- cluts.tab.update - apply(cluts.tab.update, 1, logSumExp)
    soft.probs.tab.norm <- cluts.tab.update - rowLogSumExps(cluts.tab.update)
    soft.probs.tab.norm <- exp(soft.probs.tab.norm) #get non-log scale probabilities
    #convert theta.param and pi.param back to non-log scale probabilities
    theta.param <- lapply(theta.param, exp)
    pi.param <- exp(pi.param)
  } else {
    soft.probs.tab.norm <- cluts.tab.update / rowSums(cluts.tab.update)
  }
  rownames(soft.probs.tab.norm) <- rownames(counts.l[[1]])
    #cluts.tab <- do.call(cbind, clust.prod)
    #cluts.tab.logL <- cluts.tab/rowSums(cluts.tab)
    #cluts.tab <- Reduce("*", clust.gammas.rowsums.l)
    #cluts.tab.norm <- cluts.tab/rowSums(cluts.tab)
  #}
  
  message("DONE!!!")  
  return(list(soft.pVal=soft.probs.tab.norm, log.l=unlist(log.like.l), theta.param=theta.param, pi.param=pi.param))
}  

