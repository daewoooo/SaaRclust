#' Order contigs based on shared strand states
#' 
#' This function takes \code{data.frame} of strand states per contig [rows] and per cell [columns]
#' and tries to order contigs based on shared strand states among all cells.
#'
#' @param contig.states A \code{data.frame} of strand states per contig and per cell.
#' @param dist.matrix A symetric \code{matrix} of similarities for each contig pair.
#' @param method Set one of the method to solve TSP, default: 'cheapest_insertion'.
#' @param trials Number of random trials to solve TSP.
#' @param filt.cols If set to \code{TRUE}, will remove columns with the same strand-state across all contigs.
#' @importFrom TSP TSP solve_TSP
#' @importFrom cluster daisy
#' @author David Porubsky
#' @export
#' 
orderContigsTSP <- function(contig.states=NULL, dist.matrix=NULL, method='nearest_insertion', trials=1000, filt.cols=TRUE) {
  ptm <- startTimedMessage("    Running TSP using '", method, "' method and ", trials, " trials")
  if (is.data.frame(contig.states) & is.null(dist.matrix)) {
    ## Remove columns that have the same strand state across all contigs ('uninformative cells')
    if (filt.cols) {
      mask <- apply(contig.states, 2, function(x) length(unique(x)) > 1)
      if (length(mask[mask == TRUE]) > 1) {
        contig.states.filt <- contig.states[,mask]
      } else {
        message("    Parameter 'filt.cols' would leave only one cell, skipping ...")
        contig.states.filt <- contig.states
      }
      ## Transform columns of input contig.states into a factor
      contig.states.f <- data.frame(apply(contig.states.filt, 2, as.factor))
    } else {
      ## Transform columns of input contig.states into a factor
      contig.states.f <- data.frame(apply(contig.states, 2, as.factor))
    }  
    ## Calculate distance of each pair of contigs
    dists <- as.matrix(cluster::daisy(contig.states.f))
  } else if (is.null(contig.states) & is.matrix(dist.matrix)) {
    dists <- dist.matrix
  } else {
    message("Please submit contig.states or dist.matrix. See package documentation.")
  }
  
  ## Add a dummy node with equal (max) distance to all other nodes
  dists <- cbind(dists, rep(1, nrow(dists)))
  dists <- rbind(dists, rep(1, ncol(dists)))
  rownames(dists)[nrow(dists)] <- 'dummy'
  colnames(dists)[nrow(dists)] <- 'dummy'
  
  ## Find 'best' traversal through all contigs using TSP package
  contig.tsp <- TSP::TSP(dists)
  contigs.order <- TSP::solve_TSP(contig.tsp, method=method, control=list(rep=trials))
  ## Remove the dummy node
  contigs.order <- labels(contigs.order)[-c(which(labels(contigs.order)=='dummy'))]
  ## Order contig table for export
  if (is.data.frame(contig.states) & is.null(dist.matrix)) {
    ordered.contigs.table = contig.states[contigs.order,]
    order.vector = match(contigs.order, rownames(contig.states))
  } else if (is.null(contig.states) & is.matrix(dist.matrix)) {
    ordered.contigs.table = dist.matrix[contigs.order,]
    order.vector = match(contigs.order, rownames(dist.matrix))
  }
  stopTimedMessage(ptm)
  ## Return ordered contigs
  return(list(ordered.contigs = contigs.order, 
              order.vector = order.vector,
              ordered.contigs.table = ordered.contigs.table))
}


#' Order contigs based on shared strand states
#' 
#' This function uses ContiBAIT implementation of MC algorithm to solve TSP
#' Copyright (c) 2015, Kieran O'Neill, Mark Hills, Mike Gottlieb
#' All rights reserved.
#'
#' @param contig.states A \code{data.frame} of strand states per contig and per cell.
#' @param randomAttempts Number of random trials to solve contig ordering.
#' @useDynLib SaaRclust
#' @import Rcpp TSP
#' @author Kieran O'Neill, Mark Hills, Mike Gottlieb (Modified by David Porubsky)
#' @export
#' 
orderContigsGreedy <- function(contig.states, randomAttempts=1000) {
  ptm <- startTimedMessage("    Running ContiBAIT ordering using ", randomAttempts, " random starts")
  
  best.order <- .Call('orderContigsGreedy', as.matrix(contig.states))
  best.table <- contig.states
  
  for (i in seq_len(randomAttempts)) {
    temp.table <- as.matrix(contig.states[sample(nrow(contig.states)),])
    temp.order <- .Call('orderContigsGreedy', temp.table)
    
    if (temp.order$score < best.order$score) {
      #if (verbose) {message('     -> Found better ordering!')}  
      best.order <- temp.order
      best.table <- temp.table
    }
  }
  stopTimedMessage(ptm)
  ## Return ordered contigs
  ordered.contigs <- row.names(best.table)[best.order$order]
  return(list(ordered.contigs=ordered.contigs,
              order.vector = match(ordered.contigs, rownames(contig.states)),
              ordered.contigs.table = contig.states[ordered.contigs,], 
              score = best.order$score))
}
