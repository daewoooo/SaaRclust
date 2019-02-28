#' Import pre-processed output from the minimap
#'
#' This function expects output from custom minimap test dataset that contains original locations of mapped reads in the genome.
#'
#' @param infile A path to the minimap file to load.
#' @return A \code{data.frame}.
#' @importFrom data.table fread
#' @author David Porubsky
#' @export 

importData <- function(infile=NULL) {  #TODO modify this function for input where genomic location of PB reads is unknown
  
  ptm <- startTimedMessage("Reading the data") 
  #data <- read.table(infile, header=F) #TODO test data.table package for faster data import
  
  #filetype = summary( file(infile) )$class #If it's gzipped, filetype will be gzfile
  if (summary( file(infile) )$class == 'gzfile') {
    data <- data.table::fread(paste0('zcat ',infile), header=T, verbose = F, showProgress = F)
    # select columns:
    data <- data[, .(SSreadNames, SSlibNames, SSflag, SSchrom, SSpos, strand, PBreadNames, PBflag, PBchrom, PBpos, PBreadLen, TargetCoordStart, TargetCoordend, MatchedBasesWithGaps)]
  } else {
    data <- data.table::fread(infile, header=T, verbose = F, showProgress = F)
  }
  
  #make sure strand info is represented as factor variable
  data$strand <- factor(data$strand)
  
  stopTimedMessage(ptm)
  return(data)
}


#' Import output from the minimap
#'
#' This function expects output from custom minimap test dataset that contains original locations of mapped reads in the genome.
#'
#' @param infile A path to the minimap file to load.
#' @param removeDuplicates Set to \code{TRUE} if you want to remove duplicate reads based on the flag.
#' @return A \code{data.frame}.
#' @importFrom data.table fread tstrsplit
#' @author David Porubsky

importTestData <- function(infile=NULL, removeDuplicates = TRUE) {  #TODO modify this function for input where genomic location of PB reads is unknown

  ptm <- startTimedMessage("Reading the data") 
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


#' Filter input data from the minimap
#'
#' This function filters loaded data from the minimap.
#'
#' @param inputData A \code{data.frame} loaded by \code{\link[SaaRclust]{importTestData}}
#' @param quantileSSreads A quantile range for number of SSreads mapped to PB read. (default: 0.4-0.9)
#' @param minSSlibs A range for the minimal and maximal number of StrandS libs being represented per PB read.
#' @return A filtered \code{data.frame}.
#' @importFrom dplyr group_by summarise
#' @author David Porubsky
#' @export


filterInput <- function(inputData=NULL, quantileSSreads=c(0,0.9), minSSlibs=c(20,25)) {
  
  ptm <- startTimedMessage("Filtering the data")
  
  #remove SS reads mapped with huge gaps (eg. summed gaps 100bp)
  gaps.perSS.mean <- round(mean(inputData$MatchedBasesWithGaps))
  inputData.filt <- inputData[inputData$MatchedBasesWithGaps <= gaps.perSS.mean,]
  
  #remove SS reads mapped to multiple PB reads
  #rle.SSreads <- rle(inputData.filt$SSreadNames)
  #mask <- rle.SSreads$lengths == 1
  #single.SSreads <- unique(inputData.filt$SSreadNames)[mask]
  ##multi.SSreads <- unique(inputData.filt$SSreadNames)[!mask]
  ##single <- inputData.filt[inputData.filt$SSreadNames %in% single.SSreads,]
  ##multi <- inputData.filt[inputData.filt$SSreadNames %in% multi.SSreads,]
  #inputData.filt <- inputData.filt[inputData.filt$SSreadNames %in% single.SSreads,]
  
  #get number of SS reads per PB read
  if (!is.null(quantileSSreads)) {
    SSread.perPB <- sort(table(inputData.filt$PBreadNames), decreasing = T) #this won't be needed when output will be already sorted by PBreads
  
    #filter reads based on the mean number of SS reads aligned per each PB read
    quantile.range <- quantile(SSread.perPB, probs = quantileSSreads)
    filt <- SSread.perPB >= quantile.range[1] & SSread.perPB <= quantile.range[2]
    upperQ.reads <- names(SSread.perPB)[!filt]
    maskNames <- names(SSread.perPB)[filt]
    inputData.filt <- inputData.filt[inputData.filt$PBreadNames %in% maskNames,]
  }
  
  #filter reads based on the number of SS libs per PB read [FAST]
  inputData.filt %>% dplyr::group_by(PBreadNames) %>% dplyr::summarise(counts = length(unique(SSlibNames))) -> SSlib.perPB.counts
  maskNames <- SSlib.perPB.counts$PBreadNames[SSlib.perPB.counts$counts >= minSSlibs[1] & SSlib.perPB.counts$counts <= minSSlibs[2]]
  inputData.filt <- inputData.filt[inputData.filt$PBreadNames %in% maskNames,]
  
  #filter reads based on the number of SS libs per PB read [SLOW]
  #SSlib.perPB <- split(as.factor(inputData.filt$SSlibNames), inputData.filt$PBreadNames)
  #SSlib.perPB.counts <- sapply(SSlib.perPB, function(x) length(unique(x)))
  #maskNames <- names(SSlib.perPB.counts)[SSlib.perPB.counts >= minSSlibs]
  #inputData.filt <- inputData.filt[inputData.filt$PBreadNames %in% maskNames,]
  
  stopTimedMessage(ptm)
  return(list(tab.filt=inputData.filt, upperQ.reads=upperQ.reads))
}


#' Import representative PacBio alignments.
#'
#' This function selects the best representative PacBio alignements from multiple chunks
#' in order to get the best estimates of theta values using Kmeans.
#'
#' @param inputfolder  A folder name where minimap files are stored.
#' @param numAlignments  A required number of representative PacBio alignments.
#' @inheritParams filterInput
#' @return A \code{data.frame}.
#' @author David Porubsky
#' @export 

getRepresentativeAlignments <- function(inputfolder=NULL, numAlignments=30000, quantileSSreads=c(0.2,0.8), minSSlibs=c(20,25)) {
  
  ptm <- startTimedMessage("Getting representative alignments\n") 
  file.list <- list.files(path = inputfolder, pattern = "maf.gz$", full.names = TRUE)
  
  bestAligns <- list()
  countAligns <- 0
  for (file in file.list) {
    filename <- basename(file)
    ptm <- proc.time()
    
    suppressMessages( suppressWarnings( tab.in <- importData(infile=file) ) )
    #suppressMessages( suppressWarnings( tab.in <- importTestData(infile=file, removeDuplicates=TRUE) ) )
    suppressMessages( tab.filt.l <- filterInput(inputData=tab.in, quantileSSreads=quantileSSreads, minSSlibs=minSSlibs) )
    tab.filt <- tab.filt.l$tab.filt
    
    numAligns <- length(unique(tab.filt$PBreadNames))
    
    message("    Obtained ",numAligns," representative alignments", appendLF = F); 
    stopTimedMessage(ptm)
    
    countAligns <- countAligns + numAligns
    bestAligns[[filename]] <- tab.filt
    #interupt the loop in case required amount of representative alignement was reached
    if (!is.null(numAlignments) & countAligns >= numAlignments) {
      break
    }
  }
  
  bestAligns.tab <- do.call(rbind, bestAligns)
  rownames(bestAligns.tab) <- NULL
  #bestAligns.tab <- bestAligns.tab[sample(nrow(bestAligns.tab)),] #shuffle rows in tab
  sample <- unique(bestAligns.tab$PBreadNames)[1:numAlignments] #get only required amount of representative alignments
  bestAligns.tab <- bestAligns.tab[bestAligns.tab$PBreadNames %in% sample,]
  
  return(bestAligns.tab)
}


#' Import BAM(s) and count reads
#'
#' Import aligned reads from a multiple BAM files and counts directional reads in specified genomic locations.
#' Results are stored in a \code{list} of matrices with each element of a \code{list} representing counts for single BAM file.
#'
#' @param bamfolder A folder where BAM files to be processed are stored.
#' @param chromosomes If only a subset of the chromosomes should be binned, specify them here.
#' @param pairedEndReads Set to \code{TRUE} if you have paired-end reads in your file.
#' @param bin.length A length of a bin to count reads in.
#' @return A \code{list} of matrices (columns: minus (W) and plus (C) counts; rows: genomic regions).
#' @importFrom data.table data.table
#' @author David Porubsky
#' @export

importBams <- function(bamfolder=bamfolder, chromosomes=NULL, bin.length=1000000) {
  ## List bams present in a directory
  bamfiles <- list.files(bamfolder, pattern = '.bam$', full.names = T)
  
  counts.l <- list()
  for (j in 1:length(bamfiles)) {
    bam <- bamfiles[j]
    bam.name <- basename(bam)
    message("Processing ", bam.name)
    
    ## Load reads from BAM into GRanges object
    #fragments <- bam2GRanges(file = bam, chromosomes=chromosomes, pairedEndReads = T, min.mapq = 10, keep.duplicate.reads = F, what = 'mapq')
    suppressWarnings( fragments <- readBamFileAsGRanges(file=bam, chromosomes=chromosomes, pairedEndReads=TRUE, min.mapq=10, remove.duplicate.reads=TRUE, pair2frgm=FALSE, filtAlt=TRUE) )
    ## Sort fragments by seqlevels
    fragments <- GenomicRanges::sort(fragments, ignore.strand=TRUE)
    
    if (bin.length) {
      ## Get genome bins
      bins.gr <- unlist( GenomicRanges::tileGenome(seqlengths = seqlengths(fragments), tilewidth = bin.length) )
      hits <- IRanges::findOverlaps(fragments, bins.gr, select = "first") #TODO: make sure the same read can't end up in two neighbouring bins!!!
      
      #runLength(seqnames(fragments))[which(runValue(seqnames(fragments)) == 'Super-Scaffold_365')]
      
      ## Add missing values if there is no reads in any given bin
      bin.num <- 1:length(bins.gr)
      missing.bin <- which(!bin.num %in% hits)
      
      ## Append bin location to the contig name
      if (length(missing.bin) > 0) {
        fragments$ID <- rep(as.character(bins.gr)[-missing.bin], table(hits))
      } else {
        fragments$ID <- rep(as.character(bins.gr), table(hits))
      }  
      
    } else {
      fragments$ID <- seqnames(fragments)
    }
    
    ## Transform GRanges object into a data.frame
    fragments.df <- as(fragments[,'ID'], 'data.frame')
    fragments.df$ID <- factor(fragments.df$ID, levels=as.character(bins.gr))
    fragments.df$strand <- factor(as.character(fragments.df$strand))
    
    #table(fragments.df[fragments.df$seqnames == 'Super-Scaffold_365',]$strand)
    #tmp <- fragments.df[fragments.df$seqnames == 'Super-Scaffold_365',]
    
    ## Count reads
    counts <- data.table::data.table(fragments.df)[,table(strand),by=ID]
    cov.reads <- unique(counts$ID)
    uncov.reads <- levels(cov.reads)[!levels(cov.reads) %in% cov.reads]
    counts <- rbind(matrix(counts$V1, ncol=2, byrow = T), matrix(rep(0, 2*length(uncov.reads)), ncol=2) )
    rownames(counts) <- c(as.character(unique(cov.reads)), uncov.reads)
    counts <- counts[order(match(rownames(counts),levels(fragments.df$ID))),]
    #counts.l[[j]] <- counts
    counts.l[[bam.name]] <- counts  #TEST
  }
  return(counts.l)
}
