#' Calculate probability of PacBio reads having WW, CC or WC state.
#'
#' Import counts of directional Strand-seq reads aligned to PacBio reads and calculates binomial probabilites of PacBio reads being WW, CC or WC state.
#'
#' @param minusCounts Minus (Watson) read counts aligned to PacBio reads.
#' @param plusCounts Plus (Crick) read counts aligned to PacBio reads.
#' @param alpha Estimated level of background Strand-seq reads.
#' @return A \code{list} of binomial probabilities for a given counts of plus and minus reads
#' @author David Porubsky, Maryam Ghareghani

countProbOld <- function(minusCounts, plusCounts, alpha=0.1) {

  minusCounts <- as.numeric(minusCounts)
  plusCounts <- as.numeric(plusCounts)

  #calculate that given PB read is WW
  prob.w <- (1-alpha)^minusCounts
  prob.c <- alpha^plusCounts 
  prob.ww <- prob.w * prob.c

  #calculate that given PB read is CC
  prob.w <- alpha^minusCounts 
  prob.c <- (1-alpha)^plusCounts
  prob.cc <- prob.w * prob.c

  #calculate that given PB read is WC
  prob.w <- 0.5^minusCounts
  prob.c <- 0.5^plusCounts
  prob.wc <- prob.w * prob.c

  prob.mix <- prob.wc

  prob.m <- choose(n=minusCounts+plusCounts, k=plusCounts) * cbind(prob.ww, prob.cc, prob.mix)

  return(prob.m)
}


#' Import old output from the minimap
#'
#' This function expects output from custom minimap test dataset that contains original locations of mapped reads in the genome.
#'
#' @param infile A path to the minimap file to load.
#' @return A \code{data.frame}.
#' @author David Porubsky

importOldTestData <- function(infile=NULL, removeDuplicates = TRUE) {  #TODO modify this function for input where genomic location of PB reads is unknown
  
  ptm <- startTimedMessage("Reading the data ...") 
  #data <- read.table(infile, header=F) #TODO test data.table package for faster data import
  
  #filetype = summary( file(infile) )$class #If it's gzipped, filetype will be gzfile
  if (summary( file(infile) )$class == 'gzfile') {
    data <- data.table::fread(paste0('gunzip -cq ',infile), header=F, verbose = F, showProgress = F)
  } else {
    data <- data.table::fread(infile, header=F, verbose = F, showProgress = F)
  }
  
  SSreadIDs <- as.character(data$V1)
  PBreadIDs <- as.character(data$V6)
  #SSreadIDs <- as.character(data$V6)
  #PBreadIDs <- as.character(data$V1)
  
  SSreadNames.fields <- data.table::tstrsplit(SSreadIDs, "_")
  PBreadNames.fields <- data.table::tstrsplit(PBreadIDs, "_")
  
  SSreadNames <-  SSreadNames.fields[[1]] 
  SSlibNames <- paste(SSreadNames.fields[[2]], SSreadNames.fields[[3]], SSreadNames.fields[[4]], sep="_")
  #SSlibNames <- paste(SSreadNames.fields[[2]], SSreadNames.fields[[3]], sep="_")
  
  SSflag <- SSreadNames.fields[[5]] 
  SSchrom <- SSreadNames.fields[[6]] 
  SSpos <- SSreadNames.fields[[7]] 
  #SSflag <- SSreadNames.fields[[4]] 
  #SSchrom <- SSreadNames.fields[[5]] 
  #SSpos <- SSreadNames.fields[[6]] 
  
  PBreadNames <-  paste(PBreadNames.fields[[1]], PBreadNames.fields[[2]], PBreadNames.fields[[3]], PBreadNames.fields[[4]], PBreadNames.fields[[5]], PBreadNames.fields[[6]], PBreadNames.fields[[7]], sep="_")       
  PBflag <- PBreadNames.fields[[8]] 
  PBchrom <- PBreadNames.fields[[9]] 
  PBpos <- PBreadNames.fields[[10]] 
  
  tab <- data.frame(SSreadNames=SSreadNames, 
                    SSlibNames=SSlibNames, 
                    SSflag=SSflag, 
                    SSchrom=SSchrom, 
                    SSpos=SSpos, 
                    SSreadLen=data$V2,
                    strand=factor(data$V5),
                    PBreadNames=PBreadNames,
                    PBflag=PBflag,
                    PBchrom=PBchrom,
                    PBpos=PBpos,
                    PBreadLen=data$V7,
                    MatchedBases=data$V10,
                    MatchedBasesWithGaps=data$V11,
                    stringsAsFactors = F
  )
  
  if (removeDuplicates) {
    bit.flag <- bitwAnd(1024, as.numeric(tab$SSflag))
    mask <- bit.flag == 0 	
    tab <- tab[mask,]
  }  	
  
  stopTimedMessage(ptm)
  return(tab)
}


#' Get the feature vector based on the WmiunsC ratios
#'
#' This function expects output from custom minimap test dataset that contains original locations of mapped reads in the genome.
#'
#' @param counts.l A \code{list} of directional read counts per PB read per library.
#' @inheritParams SaaRclust
#' @return A \code{list} of estimated theta values for every cluster and cell.
#' @author Maryam Ghareghani

WminusCratiosFeatures <- function(counts.l=NULL) {
  ratios.l <- list()
  for (j in 1:length(counts.l)) {
    #lib.name <- names(tab.l[j])
    #message("\tWorking on ",lib.name)
    counts <- counts.l[[j]]
    
    ratios <- (counts[,2]-counts[,1])/(counts[,2]+counts[,1]) #calculate ratio of WW reads
    ratios[is.nan(ratios)] <- 0
    ratios.l[[j]] <- ratios
  }
  
  ratios.m <- do.call(cbind, ratios.l)
  #ratios.m[ratios.m<0] <- -1
  #ratios.m[ratios.m>0] <- 1
  
  ratios.m
}

#' Get the feature vector based on the max{W,C} ratios
#'
#' This function expects output from custom minimap test dataset that contains original locations of mapped reads in the genome.
#'
#' @param counts.l A \code{list} of directional read counts per PB read per library.
#' @inheritParams SaaRclust
#' @return A \code{list} of estimated theta values for every cluster and cell.
#' @author Maryam Ghareghani

maxWandCratiosFeatures <- function(counts.l=NULL) {
  ratios.l <- list()
  for (j in 1:length(counts.l)) {
    #lib.name <- names(tab.l[j])
    #message("\tWorking on ",lib.name)
    counts <- counts.l[[j]]
    
    ratios <- counts[,2]/(counts[,2]+counts[,1])  #calculate ratio of WW reads
    ratios <- sapply(ratios, function(x) max(x, 1-x))
    ratios[is.nan(ratios)] <- 0
    ratios.l[[j]] <- ratios
  }
  
  ratios.m <- do.call(cbind, ratios.l)
  #ratios.m[ratios.m<0] <- -1
  #ratios.m[ratios.m>0] <- 1
  
  ratios.m
}

#' Get the feature vector based on both W and C ratios
#'
#' This function expects output from custom minimap test dataset that contains original locations of mapped reads in the genome.
#'
#' @param counts.l A \code{list} of directional read counts per PB read per library.
#' @inheritParams SaaRclust
#' @return A \code{list} of estimated theta values for every cluster and cell.
#' @author Maryam Ghareghani

WandCratiosFeatures <- function(counts.l=NULL) {
  ratios.l <- list()
  for (j in 1:length(counts.l)) {
    #lib.name <- names(tab.l[j])
    #message("\tWorking on ",lib.name)
    counts <- counts.l[[j]]
    
    ratios <- counts/(counts[,2]+counts[,1]) #calculate ratio of WW reads
    ratios[is.nan(ratios)] <- 0
    ratios.l[[j]] <- ratios
  }
  
  ratios.m <- do.call(cbind, ratios.l)
  #ratios.m[ratios.m<0] <- -1
  #ratios.m[ratios.m>0] <- 1
  
  ratios.m
}


#' Rescale theta values for WC cell type
#'
#' Set WC probs based on WW and CC probs distribution (if ratio of WW and CC == 1 then WC region with high prob)
#'
#' @param theta.l A \code{list} of estmated theta values for every single cell.
#' @author David Porubsky

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


#' Merge splitted clusters
#'
#' This function merges initialy divided clusters that belongs to the same chromosome.
#'
#' @param cluster2merge ...
#' @param soft.pVal ...
#' @author David Porubsky

mergeSplitedClusters <- function(cluster2merge=NULL, soft.pVal=NULL) {
  
  merged.pVals <- list()
  for (i in 1:length(cluster2merge)) {
    merge.idx <- cluster2merge[[i]]
    merge.pVal <- apply(soft.pVal[,merge.idx], 1, max)
    merged.pVals[[i]] <- merge.pVal
  }  
  
  return(do.call(cbind,merged.pVals))
} 