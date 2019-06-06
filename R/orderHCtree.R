#' Order branches of hierarchicaly clustered tree
#' 
#' This function takes as input \code{matrix} of strand states per contig [rows] and per cell [columns]
#' and attments to order branches of hierarchically clustered tree around shared nodes.
#'
#' @param sts.matrix A \code{matrix} of strand states per contig and per cell.
#' @param coinh.matrix A \code{matrix} of distances for each contig pair.
#' @return A \code{list} of hierarchically ordered contigs.
#' @author David Porubsky
#' @export
#' 
orderHCtree <- function(sts.matrix=NULL, coinh.matrix=NULL) {
  if (is.null(sts.matrix) | is.null(coinh.matrix)) {
    message("Please provide parameters of this function. See function documentation.")  
  }
  ## Remove columns that have the same strand state across all contigs
  mask <- apply(sts.matrix, 2, function(x) length(unique(x)) > 1)
  sts.matrix <- sts.matrix[,mask]
  ## Get hierarchical ordering of a given cluster
  contig.dist <- dist(sts.matrix)
  hc.clust <- hclust(contig.dist)
  hc.order <- hc.clust$order
  ## Get clustered tree cuts
  hc.tree.cuts <- cutree( hc.clust, 2:(length(hc.order) - 1) )
  hc.tree.cuts <- hc.tree.cuts[hc.order,]
  ## Loop over every tree cut and attempt to provide better order
  for (i in ncol(hc.tree.cuts):1) {
    cuts <- hc.tree.cuts[,i]
    ## Get indices of clustered subtrees
    sub.tree.idx <- which(duplicated(cuts) | duplicated(cuts, fromLast=TRUE))
    ## Split indices by clustred subtree
    sub.tree.idx.l <- split(sub.tree.idx, cuts[sub.tree.idx])
    ## Keep only set of contigs which are immediate neighbours
    mask <- lapply(sub.tree.idx.l, function(x) sum(diff(x)) < length(x))
    sub.tree.idx.l <- sub.tree.idx.l[which(mask == TRUE)]
    ## Skip iteration if there is no set of contigs to sort
    if (length(sub.tree.idx.l) == 0) {next}
    ## Loop through subtrees
    for (j in 1:length(sub.tree.idx.l)) {
      #print(j)
      sub.tree <- sub.tree.idx.l[[j]]
      ## Get indices of neighboring contigs
      neighbor1 <- sub.tree[1] - 1
      neighbor2 <- sub.tree[length(sub.tree)] + 1
      if (neighbor1 < 1 | neighbor1 > nrow(hc.tree.cuts)) {neighbor1 <- neighbor2}
      if (neighbor2 < 1 | neighbor2 > nrow(hc.tree.cuts)) {neighbor2 <- neighbor1}
      ## Get scores between neighbours
      sub.tree.names <- rownames(hc.tree.cuts)[sub.tree]
      neighbor.names <- rownames(hc.tree.cuts)[c(neighbor1, neighbor2)]
      ## Calculate distance to neighboring contig on the left
      d.n1.c1 <- coinh.matrix[neighbor.names[1], sub.tree.names[1]]
      d.n1.c2 <- coinh.matrix[neighbor.names[1], sub.tree.names[2]]
      dist.n1 <- d.n1.c1 + d.n1.c2
      ## Calculate distance to neighboring contig on the right
      d.n2.c1 <- coinh.matrix[neighbor.names[2], sub.tree.names[1]]
      d.n2.c2 <- coinh.matrix[neighbor.names[2], sub.tree.names[2]]
      dist.n2 <- d.n2.c1 + d.n2.c2
      ## Reverse order of a set of contigs if distance in flanking contigs to its immediate
      ## neigbouring contigs is smaller in reverse order
      if (d.n1.c1 > d.n1.c2 | d.n2.c2 > d.n2.c1) {
        rownames(hc.tree.cuts)[sub.tree] <- rev(rownames(hc.tree.cuts)[sub.tree])
        hc.order[sub.tree] <- rev(hc.order[sub.tree])
      }
    }
  }
  ## Export reordered an original hierarchical tree
  return(list(new.order=hc.order, new.contig.order=rownames(hc.tree.cuts)))
}  
