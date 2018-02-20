############################################################################################################################################
#This function calculates clustering accuracy over different probability cutoffs.
#Depth of coverage is reported as a cumulative numeber of bases sequenced (sum of PB read lenghts/genome size)

inputfolder <- "/media/daewoooo/WORK/Clustering_project/WholeGenomeAnalysis/SaaRclust_results_DataQualFilt/"
thresholds <- c(0, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.99)
ClustersAccuracyPerChrPerDir(inputfolder=inputfolder, thresholds=thresholds, minLib=10) -> accplt.obj
destination <- file.path(inputfolder, "accPlot_minLib0.RData") 
save(file = destination, accplt.obj)

ClustersAccuracyPerChrPerDir <- function(inputfolder=NULL, thresholds=NULL, minLib=NULL) {
  Clusters2process <- list.files(file.path(inputfolder, 'Clusters'), pattern = "clusters.RData", full.names = TRUE)
  Quals2process <- list.files(file.path(inputfolder, 'RawData'), pattern = "dataQuals.RData", full.names = TRUE)
  PBlen2process <- list.files("/media/daewoooo/WORK/Clustering_project/WholeGenomeAnalysis/PBreadLen/", pattern = "gz", full.names = TRUE)
  
  allClusters <- list()
  for (i in 1:length(Clusters2process)) {  
    data.file <- get(load(Clusters2process[i]))
    data.qual <- get(load(Quals2process[i]))
    PB.read.len <- data.table::fread(paste0('gunzip -cq ', PBlen2process[i]), header=T, verbose = F, showProgress = F)
    fileID <- basename(Clusters2process[i])
    message("Processing file: ",fileID)
    
    #Sort data quals according to PB order in clusters
    SSlib.perPB <- data.qual$SSlib.perPB
    SSlib.perPB <- SSlib.perPB[match(rownames(data.file$soft.pVal), SSlib.perPB$PBreadNames),]
    pb.minLib <- SSlib.perPB$counts
    
    #Select required PB read lengths
    PB.read.len <- PB.read.len[match(rownames(data.file$soft.pVal), PB.read.len$PBreadNames),]
    pb.readLen <- PB.read.len$PBreadLen
    
    ##check accuracy only for autosomes and sex chrmosomes
    mask <- which(grepl('^chr[0-9X][0-9]?$', data.file$PBchrom))
    
    #get clusters IDs corresponding to a given chromosome
    chr.rows <- data.file$PBchrom[mask]
    chr.flag <- data.file$PBflag[mask]
    prob.tab <- data.file$soft.pVal[mask,]
    pb.minLib <- pb.minLib[mask]
    pb.readLen <- pb.readLen[mask]
    
    #filter out duplicates
    mask <- which(chr.flag == 16 | chr.flag == 0) 
    chr.rows <- chr.rows[mask]
    chr.flag <- chr.flag[mask]
    prob.tab <- prob.tab[mask,]
    pb.minLib <- pb.minLib[mask]
    pb.readLen <- pb.readLen[mask]
    
    #Find WC cluster in all cells
    theta.sums <- Reduce("+", data.file$theta.param)
    remove.clust <- which.max(theta.sums[,3])
    #Remove probabilities for always WC cluster
    prob.tab <- prob.tab[,-remove.clust]
    
    #Remove PB reads represneted by SSlib less than minLib
    filt <- pb.minLib >= minLib
    prob.tab <- prob.tab[filt,]
    chr.rows <- chr.rows[filt]
    chr.flag <- chr.flag[filt]
    pb.readLen <- pb.readLen[filt]
    
    #get clusters IDs corresponding to a given chromosome
    Clust.IDs <- getClusterIdentityPerChrPerDir(soft.clust=prob.tab, chr.rows=chr.rows, chr.flag=chr.flag)
    
    clust.acc.l <- list()
    for (prob.th in thresholds) {
      message("    Set threshold: ", prob.th)
      max.prob <- apply(prob.tab, 1, max)
      mask <- max.prob >= prob.th
      pb.readLen.sub <- pb.readLen[mask]
      
      Clust.locations <- apply(prob.tab[mask,], 1, which.max) 
      
      #calculate clustering accuracy in comparison to expected values
      clust.acc <- Clust.locations == Clust.IDs[mask]
      acc.th <- table(clust.acc)
      
      clust.acc.l[[1+length(clust.acc.l)]] <- c(prob.th=prob.th, acc.th.match=unname(acc.th[2]), acc.th.sum=sum(acc.th), allReads=length(chr.rows), seq.bases=sum(as.numeric(pb.readLen.sub)))  
    }
    allClusters[[fileID]] <- as.data.frame( do.call(rbind, clust.acc.l) )
  } 
  #sum all counts over all data frames (per position)
  clust.acc.df <- Reduce("+", allClusters)
  
  #calcualte accuracy percentages
  clust.acc.df$prob.th <- thresholds
  clust.acc.df$th.acc <- clust.acc.df$acc.th.match / clust.acc.df$acc.th.sum
  clust.acc.df$th.clustReads <- clust.acc.df$acc.th.sum / clust.acc.df$allReads
  
  #get genome size
  library("biovizBase")
  hg38Ideogram <- getIdeogram("hg38", cytoband = FALSE)
  hg38Ideogram <- keepSeqlevels(hg38Ideogram, paste0('chr', c(1:22,'X')))
  genome.size <- sum(as.numeric(seqlengths(hg38Ideogram)))
  clust.acc.df$depth <- ceiling(clust.acc.df$seq.bases/genome.size)
  
  acc.plt <- ggplot(clust.acc.df) + geom_point(aes(x=th.acc, y=th.clustReads), color="deepskyblue4", size=10) + geom_linerange(aes(ymin=-Inf, x=th.acc, ymax=th.clustReads),color="deepskyblue4") + scale_y_continuous(limits = c(0.8, 1)) + scale_x_continuous(limits = c(0,1)) + scale_y_continuous(limits = c(0,1)) + ylab("(%) evaluated PB reads") + xlab("(%) correctly assigned PB reads") + geom_text(aes(x=th.acc, y=th.clustReads), label=c('all', thresholds[-1]), color="white") + theme_bw()  + theme_bw() + geom_text(aes(x=th.acc, y=th.clustReads+0.05), label=paste0(clust.acc.df$depth, "x"), color="black")
  message("DONE!!!")
  return(list(acc.plot=acc.plt, plot.table=clust.acc.df))
} 


############################################################################################################################################
#This function calculates clustering accuracy over different threshold for minimal SSlibs represented per PBread

inputfolder <- "/media/daewoooo/WORK/Clustering_project/WholeGenomeAnalysis/SaaRclust_results_DataQualFilt/"
minLibs <- c(0,5,10,15,20,30)
ClustersAccuracyPerChrPerDir_minLib(inputfolder=inputfolder, minLibs=minLibs) -> accPlot_minLibFilt.obj
destination <- file.path(inputfolder, " accPlot_minLibFilt.RData") 
save(file = destination, accPlot_minLibFilt.obj)

ClustersAccuracyPerChrPerDir_minLib <- function(inputfolder=NULL, minLibs=NULL) {
  Clusters2process <- list.files(file.path(inputfolder, 'Clusters'), pattern = "clusters.RData", full.names = TRUE)
  Quals2process <- list.files(file.path(inputfolder, 'RawData'), pattern = "dataQuals.RData", full.names = TRUE)
  PBlen2process <- list.files("/media/daewoooo/WORK/Clustering_project/WholeGenomeAnalysis/PBreadLen/", pattern = "gz", full.names = TRUE)
  
  allClusters <- list()
  for (i in 1:length(Clusters2process)) {  
    data.file <- get(load(Clusters2process[i]))
    data.qual <- get(load(Quals2process[i]))
    PB.read.len <- data.table::fread(paste0('gunzip -cq ', PBlen2process[i]), header=T, verbose = F, showProgress = F)
    fileID <- basename(Clusters2process[i])
    message("Processing file: ",fileID)
    
    #Sort data quals according to PB order in clusters
    SSlib.perPB <- data.qual$SSlib.perPB
    SSlib.perPB <- SSlib.perPB[match(rownames(data.file$soft.pVal), SSlib.perPB$PBreadNames),]
    pb.minLib <- SSlib.perPB$counts
    
    #Select required PB read lengths
    PB.read.len <- PB.read.len[match(rownames(data.file$soft.pVal), PB.read.len$PBreadNames),]
    pb.readLen <- PB.read.len$PBreadLen
    
    ##check accuracy only for autosomes and sex chrmosomes
    mask <- which(grepl('^chr[0-9X][0-9]?$', data.file$PBchrom))
    
    #get clusters IDs corresponding to a given chromosome
    chr.rows <- data.file$PBchrom[mask]
    chr.flag <- data.file$PBflag[mask]
    prob.tab <- data.file$soft.pVal[mask,]
    pb.minLib <- pb.minLib[mask]
    pb.readLen <- pb.readLen[mask]
    
    #filter out duplicates
    mask <- which(chr.flag == 16 | chr.flag == 0) 
    chr.rows <- chr.rows[mask]
    chr.flag <- chr.flag[mask]
    prob.tab <- prob.tab[mask,]
    pb.minLib <- pb.minLib[mask]
    pb.readLen <- pb.readLen[mask]
    
    #Find WC cluster in all cells
    theta.sums <- Reduce("+", data.file$theta.param)
    remove.clust <- which.max(theta.sums[,3])
    #Remove probabilities for always WC cluster
    prob.tab <- prob.tab[,-remove.clust]
    
    #get total amount of PB reads
    total.PBreads <- length(rownames(prob.tab)) 
    
    #Set various filer limits for minomum number of SS libs per PB read
    clust.acc.minLib.l <- list()
    for (minLib in minLibs) {
      message("   Set minLib: ", minLib)
      #Remove PB reads represneted by SSlib less than minLib
      filt <- pb.minLib >= minLib
      prob.tab.sub <- prob.tab[filt,]
      chr.rows.sub <- chr.rows[filt]
      chr.flag.sub <- chr.flag[filt]
      pb.readLen.sub <- pb.readLen[filt]
      
      #get clusters IDs corresponding to a given chromosome
      Clust.IDs <- getClusterIdentityPerChrPerDir(soft.clust=prob.tab.sub, chr.rows=chr.rows.sub, chr.flag=chr.flag.sub)
      
      clust.acc.l <- list()

      max.prob <- apply(prob.tab.sub, 1, max)
      Clust.locations <- apply(prob.tab.sub, 1, which.max) 
        
      #calculate clustering accuracy in comparison to expected values
      clust.acc <- Clust.locations == Clust.IDs
      acc.th <- table(clust.acc)

      clust.acc.minLib.l[[1+length(clust.acc.minLib.l)]] <- c(acc.th.match=unname(acc.th[2]), acc.th.sum=sum(acc.th), allReads=total.PBreads, seq.bases=sum(as.numeric(pb.readLen.sub)), minLib=minLib)
    } 
    allClusters[[fileID]] <- as.data.frame( do.call(rbind, clust.acc.minLib.l) )
  }
  
  #sum all counts over all data frames (per position)
  clust.acc.df <- Reduce("+", allClusters)
  
  #calcualte accuracy percentages
  clust.acc.df$minLib <- minLibs
  clust.acc.df$th.acc <- clust.acc.df$acc.th.match / clust.acc.df$acc.th.sum
  clust.acc.df$th.clustReads <- clust.acc.df$acc.th.sum / clust.acc.df$allReads
  
  #get genome size
  library("biovizBase")
  hg38Ideogram <- getIdeogram("hg38", cytoband = FALSE)
  hg38Ideogram <- keepSeqlevels(hg38Ideogram, paste0('chr', c(1:22,'X')))
  genome.size <- sum(as.numeric(seqlengths(hg38Ideogram)))
  clust.acc.df$depth <- ceiling(clust.acc.df$seq.bases/genome.size)
  
  #acc.plt <- ggplot(clust.acc.df) + geom_point(aes(x=th.acc, y=th.clustReads), color="deepskyblue4", size=10) + geom_linerange(aes(ymin=-Inf, x=th.acc, ymax=th.clustReads),color="deepskyblue4") + scale_x_continuous(limits = c(0,1)) + scale_y_continuous(limits = c(0,1)) + ylab("(%) evaluated PB reads") + xlab("(%) correctly assigned PB reads") + geom_text(aes(x=th.acc, y=th.clustReads), label=c('all', minLibs[-1]), color="white") + theme_bw()
  acc.plt <- ggplot(clust.acc.df) + geom_point(aes(x=th.acc, y=th.clustReads), color="deepskyblue4", size=10) + geom_linerange(aes(ymin=-Inf, x=th.acc, ymax=th.clustReads),color="deepskyblue4") + scale_x_continuous(limits = c(0,1)) + scale_y_continuous(limits = c(0,1)) + ylab("(%) evaluated PB reads") + xlab("(%) correctly assigned PB reads") + geom_text(aes(x=th.acc, y=th.clustReads), label=c('all', minLibs[-1]), color="white") + theme_bw() + geom_text(aes(x=th.acc+0.05, y=th.clustReads), label=paste0(clust.acc.df$depth, "x"), color="black")
  
  message("DONE!!!")
  return(list(acc.plot=acc.plt, plot.table=clust.acc.df))
} 



############################################################################################################################################
#In this part we itetrate through PBread lenghs to report different values: Total depth of coverage as well as PB read length distribution

#Get total depth of coverage before filtering of 10% of PB reads with highest SS read counts
PBlen2process <- list.files("/media/daewoooo/WORK/Clustering_project/WholeGenomeAnalysis/PBreadLen/", pattern = "gz", full.names = TRUE)
cov.bases.all <- list()
reads.all <- list()
read.len.dist <- list()
for (i in 1:length(PBlen2process)) {  
  PB.read.len <- data.table::fread(paste0('gunzip -cq ', PBlen2process[i]), header=T, verbose = F, showProgress = F)
  cov.bases.all[[i]] <- sum( as.numeric( PB.read.len$PBreadLen[!duplicated(PB.read.len$PBreadNames)] ) )
  read.len.dist[[i]] <- as.numeric( PB.read.len$PBreadLen[!duplicated(PB.read.len$PBreadNames)] )
  reads.all[[i]] <- length(PB.read.len$PBreadLen[!duplicated(PB.read.len$PBreadNames)])
}  

#plot PB read length distribution
read.len.dist.all <- unlist(read.len.dist)
read.len.dist.all.df <- as.data.frame(read.len.dist.all)
PBread.lenDist <- ggplot(read.len.dist.all.df, aes(x=read.len.dist.all)) + geom_histogram(binwidth = 100, fill="red") + xlab("PacBio read length (bp)") + ylab("Frequency") + scale_x_continuous(breaks = c(10000, 20000, 40000, 60000, 80000)) + scale_y_continuous(labels=comma) 
save.data <- list(PBread.lenDist=PBread.lenDist, read.len.dist.all.df=read.len.dist.all.df)
save(file = "/media/daewoooo/WORK/Clustering_project/WholeGenomeAnalysis/PBreadLen/PBreadLen.RData", save.data)

total.cov.bases <- sum(unlist(cov.bases.all)) 
#get genome size
library("biovizBase")
hg38Ideogram <- getIdeogram("hg38", cytoband = FALSE)
hg38Ideogram <- keepSeqlevels(hg38Ideogram, paste0('chr', c(1:22,'X')))
genome.size <- sum(as.numeric(seqlengths(hg38Ideogram)))
total.cov.bases/genome.size


############################################################################################################################################
#This function calculates distribution of minimal SSlibs represented per PBread in sets of correctly and incorrectly assigned PB reads.

inputfolder <- "/media/daewoooo/WORK/Clustering_project/WholeGenomeAnalysis/SaaRclust_results_DataQualFilt/"
boxplotDistSSlibsPerPB(inputfolder=inputfolder, thresholds=0) -> boxplt.obj
destination <- file.path(inputfolder, "boxplot.RData") 
save(file = destination, boxplt.obj)

boxplotDistSSlibsPerPB <- function(inputfolder=NULL, thresholds=NULL) {
  Clusters2process <- list.files(file.path(inputfolder, 'Clusters'), pattern = "clusters.RData", full.names = TRUE)
  Quals2process <- list.files(file.path(inputfolder, 'RawData'), pattern = "dataQuals.RData", full.names = TRUE)
  
  allSSlib.perPB <- list()
  for (i in 1:length(Clusters2process)) {  
    data.file <- get(load(Clusters2process[i]))
    data.qual <- get(load(Quals2process[i]))
    fileID <- basename(Clusters2process[i])
    message("Processing file: ",fileID)
    
    #Sort data quals according to PB order in clusters
    SSlib.perPB <- data.qual$SSlib.perPB
    SSlib.perPB <- SSlib.perPB[match(rownames(data.file$soft.pVal), SSlib.perPB$PBreadNames),]
    
    #check accuracy only for autosomes and sex chrmosomes
    mask <- which(grepl('^chr[0-9X][0-9]?$', data.file$PBchrom))
    
    #get clusters IDs corresponding to a given chromosome
    chr.rows <- data.file$PBchrom[mask]
    chr.flag <- data.file$PBflag[mask]
    prob.tab <- data.file$soft.pVal[mask,]
    
    #filter out duplicates
    mask <- which(chr.flag == 16 | chr.flag == 0) 
    chr.rows <- chr.rows[mask]
    chr.flag <- chr.flag[mask]
    prob.tab <- prob.tab[mask,]
    SSlib.perPB <- SSlib.perPB[mask,]
    
    #Find WC cluster in all cells
    theta.sums <- Reduce("+", data.file$theta.param)
    remove.clust <- which.max(theta.sums[,3])
    #Remove probabilities for always WC cluster
    prob.tab <- prob.tab[,-remove.clust]
    
    #get clusters IDs corresponding to a given chromosome
    Clust.IDs <- getClusterIdentityPerChrPerDir(soft.clust=prob.tab, chr.rows=chr.rows, chr.flag=chr.flag)
    
    SSlib.perPB.dist <- list()
    for (prob.th in thresholds) {
      max.prob <- apply(prob.tab, 1, max)
      mask <- max.prob >= prob.th
      Clust.locations <- apply(prob.tab[mask,], 1, which.max)  
      
      #calculate clustering accuracy in comparison to expected values
      clust.acc <- Clust.locations == Clust.IDs[mask]
      acc.th <- table(clust.acc)
      
      #Split counts of SSlib per PB by accuracy vector (correct vs incorrect PB read assignemnts)
      SSlib.perPB.dist[[1+length(SSlib.perPB.dist)]] <- split(SSlib.perPB$counts[mask], clust.acc)
    }
    unlist( sapply(SSlib.perPB.dist, function(x) x['TRUE']) ) -> trues
    unlist( sapply(SSlib.perPB.dist, function(x) x['FALSE']) )-> falses
    allSSlib.perPB[[fileID]] <- list(trues=trues, falses=falses)
  } 
  unlist( sapply(allSSlib.perPB, function(x) x['trues']), use.names = F ) -> trues
  unlist( sapply(allSSlib.perPB, function(x) x['falses']), use.names = F )-> falses
  ID <- rep(c('Correct', 'Incorrect'), c(length(trues), length(falses)))
  counts <- c(trues, falses)
  allSSlib.perPB.df <- data.frame(counts=counts, ID=ID)
  #Plot 
  box.plt <- ggplot(allSSlib.perPB.df, aes(x=ID, y=counts, fill=ID)) + geom_boxplot(outlier.colour="red") + scale_fill_manual(values = c("darkolivegreen3" ,"darkgoldenrod1"), guide="none") + xlab("") + ylab("# of Strand-seq libraries per PB read") + theme_bw()
  
  return(list(acc.plot=box.plt, plot.table=allSSlib.perPB.df))
} 


############################################################################################################################################
#This function calculates probability ranking of true clusters. 
#Rank1 = max probability being the true cluster
#Rank2 = second best probability being the true cluster

inputfolder <- "/media/daewoooo/WORK/Clustering_project/WholeGenomeAnalysis/SaaRclust_results_DataQualFilt/Clusters/"
accuracyRanking(inputfolder = inputfolder ) -> rankingPlt.obj
destination <- file.path(inputfolder, "rankingPlot.RData") 
save(file = destination, rankingPlt.obj)

accuracyRanking <- function(inputfolder=NULL) {
  files2process <- list.files(inputfolder, pattern = "clusters.RData", full.names = TRUE)
  
  allCluster.ranks <- list()
  all.prob.trueClust <- list()
  for (file in files2process) {
    data.file <- get(load(file))
    fileID <- basename(file)
    num.clusters <- length(data.file$pi.param)
    
    #check accuracy only for autosomes and sex chrmosomes
    mask <- which(grepl('^chr[0-9X][0-9]?$', data.file$PBchrom))
    
    #get clusters IDs corresponding to a given chromosome
    chr.rows <- data.file$PBchrom[mask]
    chr.flag <- data.file$PBflag[mask]
    prob.tab <- data.file$soft.pVal[mask,]
    
    #filter out duplicates
    mask <- which(chr.flag == 16 | chr.flag == 0) 
    chr.rows <- chr.rows[mask]
    chr.flag <- chr.flag[mask]
    prob.tab <- prob.tab[mask,]
    
    #Find WC cluster in all cells
    theta.sums <- Reduce("+", data.file$theta.param)
    remove.clust <- which.max(theta.sums[,3])
    #Remove probabilities for always WC cluster
    prob.tab <- prob.tab[,-remove.clust]
    
    Clust.IDs <- getClusterIdentityPerChrPerDir(soft.clust=prob.tab, chr.rows=chr.rows, chr.flag=chr.flag)
    
    #for each row of probability matrix select max probability corresponding to true cluster id
    prob.trueClust <- sapply(1:nrow(prob.tab), function(x) prob.tab[x,Clust.IDs[x]])
    rank.acc <- sapply(1:nrow(prob.tab), function(x) which(sort(prob.tab[x,], decreasing = T) == prob.trueClust[x]))
    rank.acc <- unlist(rank.acc)
    #create empty vector to store data
    ranks.store <- rep(0, num.clusters)
    names(ranks.store) <- 1:num.clusters
    #count ranks per file
    table.ranks <- sort(table(rank.acc), decreasing = T)
    #store table of ranks
    ranks.store[names(table.ranks)] <- as.numeric(table.ranks)
    allCluster.ranks[[fileID]] <- ranks.store
    #store max prob for true cluster
    all.prob.trueClust[[fileID]] <- prob.trueClust
  }
  allCluster.ranks.sums <- Reduce("+", allCluster.ranks)
  
  #plotting
  table.ranks.df <- data.frame(rank.acc=names(allCluster.ranks.sums), Freq=allCluster.ranks.sums)
  table.ranks.df$rank.acc <- factor( table.ranks.df$rank.acc, levels= table.ranks.df$rank.acc)
  ranking.plt <- ggplot(table.ranks.df, aes(x=rank.acc, y=Freq)) + geom_bar(fill="chartreuse4", stat="identity") + xlab("Probability ranking of true cluster") + ylab("Frequency") + theme_bw() + scale_y_continuous(labels=comma)
  
  #get accuracy measure (sum probabilities of true cluster divided by all PB reads)
  all.prob.trueClust.v <- unlist(all.prob.trueClust)
  overal.acc <- sum(all.prob.trueClust.v)/length(all.prob.trueClust.v)
  
  return(list(ranking.plt=ranking.plt, ranking.table=table.ranks.df,  overal.acc=overal.acc))
} 


############################################################################################################################################
#This function calculates various data quality measures such as, distribution of SSlibs per PB raed or distribtuon of SSreads per SSlib per PBread

#load required libraries
library(scales)
library(ggplot2)
library(cowplot)
inputfolder <- "/media/daewoooo/WORK/Clustering_project/WholeGenomeAnalysis/SaaRclust_results_DataQualUnfilt/RawData/"
plotDataQualMeasures(inputfolder) -> dataQual.plt
destination <- file.path(inputfolder, "dataQual.RData") 
save(file = destination, dataQual.plt)

plotDataQualMeasures <- function(inputfolder=NULL) {
  files2process <- list.files(inputfolder, pattern = "dataQuals.RData", full.names = TRUE)
  
  SSreads.perPB <- list()
  SSlib.perPB <- list()
  SSreads.perlib.perPB <- list()
  PBreadLenDist <- list()
  #Load all data from all chunks
  for (file in files2process) {
    data.file <- get(load(file))
    fileID <- basename(file)
    message("Processing file: ",fileID)
    
    SSreads.perPB[[fileID]] <- data.file$SSreads.perPB
    SSlib.perPB[[fileID]] <- data.file$SSlib.perPB
    PBreadLenDist[[fileID]] <- data.file$PBreadLenDist
    #get counts and export only max 100
    SSreads.perlib.perPB[[fileID]] <- sort(table(data.file$SSreads.perlib.perPB),decreasing = T)[1:100]
  }
  SSreads.perPB.all <- do.call(c, SSreads.perPB)
  SSlib.perPB.all <- do.call(rbind, SSlib.perPB)
  SSreads.perlib.perPB.all <- Reduce("+", SSreads.perlib.perPB) #sum up counts across all chunks
  
  #Get summary of PB read length distribtion over all chunks
  max <- which.max(sapply(PBreadLenDist, nrow))
  size.ids <- PBreadLenDist[[max]]$midpoints
  
  addMissingLengths <- function(x) {
    ids2add <- size.ids[!size.ids %in% x[[1]]]
    if (!length(ids2add) == 0) {  
      add.data <- data.frame(midpoints=ids2add, freq=rep(0, length(ids2add)))
      merged.data <- rbind(x, add.data)
      return(merged.data[,2])
    } else {
      return(x[,2])
    }  
  }
  
  PBreadLenDist.modif <- lapply(PBreadLenDist, addMissingLengths)
  PBreadLenDist.modif <- PBreadLenDist.modif[lengths(PBreadLenDist.modif) == length(size.ids)]
  PBreadLenDist.all <- Reduce("+",  PBreadLenDist.modif)
  PPBreadLenDist.df <- data.frame(id=size.ids, counts=PBreadLenDist.all)
  
  #Plot data quality measures
  SSreads.perPB.all.df <- as.data.frame(SSreads.perPB.all)
  quantil0.9 <- quantile(SSreads.perPB.all.df$SSreads.perPB.all, probs = 0.9)
  #SSreads.perPB.all.df <- data.frame(SSreads.perPB=SSreads.perPB.all.df[SSreads.perPB.all.df$SSreads.perPB.all <= quantil0.9,]) 
  SSreads.perPB.plt <- ggplot(SSreads.perPB.all.df, aes(x=SSreads.perPB.all)) + geom_histogram(binwidth = 1000, fill="red") + xlab("# of Strand-seq reads per PB read") + ylab("Frequency (log scale)") + scale_y_log10(labels=comma)
  
  SSlib.perPB.all.df <- data.frame(SSlib.perPB=SSlib.perPB.all$counts)
  SSlib.perPB.plt <- ggplot(SSlib.perPB.all.df, aes(x=SSlib.perPB)) + geom_histogram(binwidth = 1, fill="red") + xlab("# of Strand-seq libraries per PB read") + ylab("Frequency") + scale_y_continuous(labels=comma) 
  
  SSreads.perlib.perPB.df <- as.data.frame(SSreads.perlib.perPB.all)
  SSreads.perlib.perPB.df <- SSreads.perlib.perPB.df[SSreads.perlib.perPB.df$Var1 %in% c(1:20),]
  SSreads.perlib.perPB.plt <- ggplot(SSreads.perlib.perPB.df, aes(x=Var1, y=Freq)) + geom_bar(fill="red", stat="identity") + xlab("# of Strand-seq reads per PB read per library") + ylab("Frequency") + scale_y_continuous(labels=comma) 
  
  PPBreadLenDist.plt <- ggplot(PPBreadLenDist.df, aes(x=id, y=counts)) + geom_bar(fill="red", stat="identity") + xlab("PacBio read length (bp)") + ylab("Frequency") + scale_x_continuous(breaks = c(10000, 20000, 40000, 60000, 80000)) + scale_y_continuous(labels=comma) 
  
  #Merge plots
  main.plt1 <- plot_grid(SSreads.perPB.plt, SSlib.perPB.plt, SSreads.perlib.perPB.plt, PPBreadLenDist.plt, nrow = 1)
  main.plt2 <- plot_grid(SSreads.perPB.plt, SSlib.perPB.plt, SSreads.perlib.perPB.plt, PPBreadLenDist.plt, nrow = 2)
  return(list(main.plt1=main.plt1, main.plt2=main.plt2, quantil0.9=quantil0.9))
}  

#######################################################################################################################################
#Plot distribution of mapped StrandS reads along a single PB read

minimap.file <- "/media/daewoooo/WORK/Clustering_project/WholeGenomeAnalysis/SaaRclust_test/SaaRclust_results/TrashBin/NA12878_WashU_PBreads_chunk9000_upperQreads_upperQreads.gz"
suppressWarnings( tab.in <- importData(infile = minimap.file, removeDuplicates = TRUE) )

#tab.filt  %>% group_by(PBreadNames) %>% summarise(counts = length(unique(SSlibNames))) -> SSlib.perPB.counts
#SSreads.perPB.l <- split(tab.filt$SSlibNames, tab.filt$PBreadNames)
#SSreads.perlib.perPB <- sapply(SSreads.perPB.l, function(x) rle(sort(x))$lengths)

#split minimap table by Pacbio reads
PBhigh.count <- split(tab.in, tab.in$PBreadNames)
counts <- sapply(PBhigh.count, nrow)
#get 100 PB reads with the highest StrandS read counts
readIDs <- names( sort(counts, decreasing = T)[1:100] )
PBhigh.count <- PBhigh.count[readIDs]

for (i in 1:length(PBhigh.count)) {
  PBread.aligns <- PBhigh.count[[i]]
  readID <- names(PBhigh.count[i])
  message("Plotting read ", readID)
  plt <- plotReadAlignments(minimap.tab=PBread.aligns)
  
  #save plot
  folder <- "/home/daewoooo/Desktop/TestPlots"
  numAlignedStrandSreads <- nrow(PBread.aligns)
  filename <- paste0("read",i,"_numAlign",numAlignedStrandSreads,".pdf")
  destination <- file.path(folder, filename)
  ggsave(plot = plt , filename = destination, limitsize = F, width = 10, height = 150)
}


##########################################################################################################################################
#Process reads with very high counts of StrandS reads

upperQreads <- "/media/daewoooo/WORK/Clustering_project/WholeGenomeAnalysis/SaaRclust_test/SaaRclust_results/TrashBin/NA12878_WashU_PBreads_chunk9000_upperQreads_upperQreads.gz"
HC.input <- "/media/daewoooo/WORK/Clustering_project/WholeGenomeAnalysis/SaaRclust_test/SaaRclust_results/Clusters/hardClusteringResults_30000.RData"
processUpperQreads <- function(minimap.file=upperQreads, alpha=0.01, theta.param=NULL, pi.param=NULL, HC.input=HC.input)
  ### Read in minimap output file ###  
  suppressWarnings( tab.in <- importData(infile = minimap.file, removeDuplicates = TRUE) )

  ### Load Hard clustering results and initialize parameters of EM algorithm [temporary solution for snakemake]
  if (is.null(theta.param) | is.null(pi.param)) {
    if (!file.exists(HC.input)) {
      stop("Hard clustering results not available!!!")
    }    
    hard.clust.results <- get(load(HC.input))
    
    #Initialize theta parameter
    theta.param <- hard.clust.results$theta.param
    #Initialize pi parameter
    pi.param <- hard.clust.results$pi.param
  }
  
  #use PB read names as factor
  tab.in$PBreadNames <- factor(tab.in$PBreadNames, levels=unique(tab.in$PBreadNames))
  #split data by SS library
  tab.l <- split(tab.in, tab.in$SSlibNames)
  
  #### Count directional reads ###
  counts.l <- countDirectionalReads(tab.l)
  
  libs.IDs <- unique(tab.in$SSlibNames)
  
  clusts.perPB.perCell <- list()
  #loop over all cells and to calculate binom probs
  for (j in 1:length(counts.l)) {
    #get counts of W and C reads per PB read
    counts <- counts.l[[j]]
    #get theta param for given cell
    params <- theta.param[[j]]
    #calc BN probs
    BN.probs <- countProb(minusCounts = counts[,1], plusCounts = counts[,2], alpha=alpha)
    BN.probs.norm <-  BN.probs/rowSums(BN.probs)
    #calc WC ratio
    ww.ratio <- (counts[,2]-counts[,1])/(counts[,2]+counts[,1]) #calculate ratio of WW reads
    read.states <- ww.ratio
    
    
    #get strand state as the highest prob state
    read.states <- apply(BN.probs.norm, 1, which.max)
    clust.states <- apply(params, 1, which.max)
    
    ww.clust.IDs <- which(clust.states == 1)
    cc.clust.IDs <- which(clust.states == 2)
    
    #get possible clusters for PB reads with ww or cc state
    clusts.perPB <- rep(list(0),length(read.states))
    clusts.perPB[which(read.states == 1)] <- rep(list(ww.clust.IDs), length(which(read.states == 1)))
    clusts.perPB[which(read.states == 2)] <- rep(list(cc.clust.IDs), length(which(read.states == 2)))
    
    clusts.perPB.perCell[[j]] <- clusts.perPB
  }  
    
  clust.summary.perPB <- list()
  for (i in 1:nrow(counts.l[[1]])) {
    PB.clust.IDs <- lapply(clusts.perPB.perCell, `[[`, i)  
    PB.clust.IDs <- PB.clust.IDs[lengths(PB.clust.IDs)>1]
    id <- rownames(counts.l[[1]])[i]
    clust.summary.perPB[[id]] <- sort(table(unlist(PB.clust.IDs)), decreasing = T)
  }  
    
  
#########################################################################################################################################
# Plot clustered data

#load clustered data
clustered.file <- "/media/daewoooo/WORK/Clustering_project/WholeGenomeAnalysis/SaaRclust_test/SaaRclust_results/Clusters/NA12878_WashU_PBreads_chunk9000_clusters.RData"  
soft.clust.obj <- get(load(clustered.file))    
  
#Find WC cluster in all cells
theta.sums <- Reduce("+", soft.clust.obj$theta.param)
remove.clust <- which.max(theta.sums[,3])
theta.param.filt <- lapply(soft.clust.obj$theta.param, function(x) x[-remove.clust,])
clust.order <- findClusterPartners(theta.param=theta.param.filt)
soft.pVal.filt <- soft.clust.obj$soft.pVal[,-remove.clust]
clust.order <- findClusterPartners(theta.param=theta.param.filt)
  
### Plotting data ###
#plot likelihood function
logl.df <- data.frame(log.l=soft.clust.obj$log.l, iter=1:length(soft.clust.obj$log.l))
logl.diff <- round(max(logl.df$log.l) - min(logl.df$log.l), digits = 0)
logL.plt <- ggplot(logl.df) + geom_line(aes(x=iter, y=log.l), color="red") + xlab("EM iterations") + ylab("Likelihood function") + annotate("text",  x=Inf, y = Inf, label = paste("Diff: ",logl.diff), vjust=1, hjust=1)
#plot cluster accuracy
#acc.plt <- plotClustAccuracy(pVal.df = soft.clust.df, num.clusters = num.clusters) Merge with Maryam's function
#plot theta values
theta.plt <- plotThetaEstimates(theta.param=soft.clust.obj$theta.param, title=fileID)
#plot heatmap
soft.clust.df <- as.data.frame(soft.pVal.filt)
soft.clust.df$PBreadNames <- rownames(soft.clust.obj$soft.pVal)
soft.clust.df$PBchrom <- soft.clust.obj$PBchrom
soft.clust.df$PBflag <- soft.clust.obj$PBflag
hm.plt <- plotHeatmap(pVal.df=soft.clust.df, colOrder=clust.order, num.clusters=num.clusters)
  
#Save plots
#destination <- file.path(plots.store, paste0(fileID, "_logL.pdf"))
#ggsave(filename = destination, plot = logL.plt, width = 8, height = 5)
#destination <- file.path(plots.store, paste0(fileID, "_thetaEstim.pdf"))
#ggsave(filename = destination, plot = theta.plt, width = 20, height = 20)
#destination <- file.path(plots.store, paste0(fileID, "_heatmap.pdf"))
#pdf(destination, width = 15, height = 10) 
#hm.plt
#dev.off()
#destination <- file.path(plots.store, paste0(fileID, "_logL.RData"))
#save(file = destination, logL.plt)
#destination <- file.path(plots.store, paste0(fileID, "_thetaEstim.RData"))
#save(file = destination, theta.plt)
#destination <- file.path(plots.store, paste0(fileID, "_heatmap.RData"))
#save(file = destination, hm.plt)

#######################################################################################################################################
#get reads assigned to 'garbage cluster'

#load clustered data
clustered.file <- "/media/daewoooo/WORK/Clustering_project/WholeGenomeAnalysis/SaaRclust_test/SaaRclust_results/Clusters/NA12878_WashU_PBreads_chunk9000_clusters.RData"  
soft.clust.obj <- get(load(clustered.file))  

max.clustIDs <- apply(soft.clust.obj$soft.pVal,1, which.max)
garbage.reads <- rownames(soft.clust.obj$soft.pVal)[max.clustIDs == remove.clust]
minimap.file <- "/media/daewoooo/WORK/Clustering_project/WholeGenomeAnalysis/SaaRclust_test/NA12878_WashU_PBreads_chunk9000.maf.gz"
suppressWarnings( tab.in <- importData(infile = minimap.file, removeDuplicates = TRUE) )
tab.garbage.reads <- tab.in[tab.in$PBreadNames %in% garbage.reads,]
tab.garbage.reads.locations <- unique( tab.garbage.reads[,c('PBreadNames', 'PBchrom', 'PBpos', 'PBreadLen')] )
tab.garbage.reads.locations$PBpos <- as.numeric(tab.garbage.reads.locations$PBpos)
tab.garbage.reads.locations$PBreadLen <- as.numeric(tab.garbage.reads.locations$PBreadLen)
tab.garbage.reads.locations <- tab.garbage.reads.locations[!is.na(tab.garbage.reads.locations$PBpos),]
tab.garbage.reads.locations.gr <- GRanges(seqnames=tab.garbage.reads.locations$PBchrom, range=IRanges(start=tab.garbage.reads.locations$PBpos, width=tab.garbage.reads.locations$PBreadLen), readnames=tab.garbage.reads.locations$PBreadNames)

seg.dup <- read.table("/media/daewoooo/WORK/Clustering_project/WholeGenomeAnalysis/SaaRclust_test/grch38_superdups.max_identity.bed", header=F)
seg.dup.gr <-  GRanges(seqnames=seg.dup$V1, ranges=IRanges(start=seg.dup$V2, end=seg.dup$V3))
