#' Calculate probability of PacBio reads having WW, CC or WC state.
#'
#' Import counts of directional Strand-seq reads aligned to PacBio reads and calculates binomial probabilites of PacBio reads being WW, CC or WC state.
#'
#' @param minusCounts Minus (Watson) read counts aligned to PacBio reads.
#' @param plusCounts Plus (Crick) read counts aligned to PacBio reads.
#' @param alpha Estimated level of background Strand-seq reads.
#' @return A \code{list} of binomial probabilities for a given counts of plus and minus reads
#' @author David Porubsky, Maryam Ghareghani
#' @export

countProb <- function(minusCounts, plusCounts, alpha=0.1) {
  
  #TODO: it's equivalent to dbinom function 
  
  minusCounts <- as.numeric(minusCounts)
  plusCounts <- as.numeric(plusCounts)
  
  #calculata that given PB read is WW
  prob.w <- (1-alpha)^minusCounts
  prob.c <- alpha^plusCounts 
  prob.ww <- prob.w * prob.c
  
  #calculata that given PB read is CC
  prob.w <- alpha^minusCounts 
  prob.c <- (1-alpha)^plusCounts
  prob.cc <- prob.w * prob.c
  
  #calculata that given PB read is WC
  prob.w <- 0.5^minusCounts
  prob.c <- 0.5^plusCounts
  prob.wc <- prob.w * prob.c
  
  prob.mix <- prob.wc
  
  prob.m <- choose(n=minusCounts+plusCounts, k=plusCounts) * cbind(prob.ww, prob.cc, prob.mix)
  
  return(prob.m)
}

#' Calculate gamma function of EM algorithm
#'
#' Takes binomial probabilites for every cluster ...
#'
#' @param clust.prob A \code{list} of binomial probabilities (multiplied by theta & pi) of each PB read in every possible cluster
#' @param pi Paramter estimate of the size of each cluster
#' @return A \code{list} of ...
#' @author David Porubsky, Maryam Ghareghani
#' @export

gammaFunction <- function(clust.prob=NULL, pi=NULL, cellNum=NULL) {
  
  clust.sums <- rowSums(Reduce("+", clust.prob))
  tab.colsums <- list()
  tab.rowsums <- list()
  #loop over all clusters
  for (i in 1:length(clust.prob)) {
    clust.p <- clust.prob[[i]]
    clust.gamma <- ( clust.p*( pi[i]^(1/cellNum) ) ) / clust.sums
    clust.gamma.colsums <- colSums(clust.gamma) # take the sum over all PB reads  
    #clust.gamma.rowsums <- rowSums(clust.gamma)
    tab.colsums[[i]] <- clust.gamma.colsums
    #tab.rowsums[[i]] <- clust.gamma.rowsums
  }
  gammas.colsums <- do.call(rbind, tab.colsums)
  #gammas.rowsums <- do.call(cbind, tab.rowsums)
  
  return(gammas.colsums)
  #return(list(gammas.colsums=gammas.colsums, gammas.rowsums=gammas.rowsums))
} 
