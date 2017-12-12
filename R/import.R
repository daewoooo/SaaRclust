#' Import output from the minimap
#'
#' This function expects output from custom minimap test dataset that contains original locations of mapped reads in the genome.
#'
#' @param infile A path to the minimap file to load.
#' @return A \code{data.frame}.
#' @author David Porubsky
#' @export

importTestData <- function(infile=NULL) {  #TODO modify this function for input where genomic location of PB reads is unknown

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
  #SSlibNames <- paste(SSreadNames.fields[[2]], SSreadNames.fields[[3]], SSreadNames.fields[[4]], sep="_")
  SSlibNames <- paste(SSreadNames.fields[[2]], SSreadNames.fields[[3]], sep="_")

  #SSflag <- SSreadNames.fields[[5]] 
  #SSchrom <- SSreadNames.fields[[6]] 
  #SSpos <- SSreadNames.fields[[7]] 
  SSflag <- SSreadNames.fields[[4]] 
  SSchrom <- SSreadNames.fields[[5]] 
  SSpos <- SSreadNames.fields[[6]] 

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
                  MatchedBases=data$V10,
                  MatchedBasesWithGaps=data$V11,
                  stringsAsFactors = F
  )
  stopTimedMessage(ptm)
  return(tab)
}


#' Filter input data from the minimap
#'
#' This function filters loaded data from the minimap.
#'
#' @param inputData A \code{data.frame} loaded by \code{\link[SaaRclust]{importTestData}}
#' @param quantileSSreads A quantile range for number of SSreads mapped to PB read. (default: 0.4-0.6)
#' @return A filtered \code{data.frame}.
#' @author David Porubsky
#' @export

filterInput <- function(inputData=NULL, quantileSSreads=c(0.1,0.9)) {
  
  ptm <- startTimedMessage("Filtering the data ...") 
  #get number of SS reads per PB read
  SSread.perPB <- sort(table(inputData$PBreadNames), decreasing = T) #this won't be needed when output will be already sorted by PBreads

  #filter reads based on the mean number of SS reads aligned per each PB read
  quantile.range <- quantile(SSread.perPB, probs = quantileSSreads)
  maskNames <- names(SSread.perPB)[SSread.perPB >= quantile.range[1] & SSread.perPB <= quantile.range[2]]
  inputData.filt <- inputData[inputData$PBreadNames %in% maskNames,]
  
  #filter reads based on the number of SS libs per PB read
  SSlib.perPB <- split(as.factor(inputData.filt$SSlibNames), inputData.filt$PBreadNames)
  SSlib.perPB.counts <- sapply(SSlib.perPB, function(x) length(unique(x)))
  maskNames <- names(SSlib.perPB.counts)[SSlib.perPB.counts >= median(SSlib.perPB.counts)]
  inputData.filt <- inputData[inputData$PBreadNames %in% maskNames,]
  
  #remove SS reads mapped with huge gaps (eg. summed gaps 100bp)
  gaps.perSS.mean <- round(mean(inputData.filt$MatchedBasesWithGaps))
  inputData.filt <- inputData.filt[inputData.filt$MatchedBasesWithGaps <= gaps.perSS.mean,]
  
  stopTimedMessage(ptm)
  return(inputData.filt)
}


