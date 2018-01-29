inputfolder <- "/media/daewoooo/WORK/Clustering_project/WholeGenomeAnalysis/SaaRclust_results_DataQualFilt/"
thresholds <- c(0, 0.1, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)
ClustersAccuracyPerChrPerDir(inputfolder=inputfolder, thresholds=thresholds, minLib=10) -> accplt.obj
destination <- file.path(inputfolder, "accPlot_minLin10.RData") 
save(file = destination, accplt.obj)

ClustersAccuracyPerChrPerDir <- function(inputfolder=NULL, thresholds=NULL, minLib=NULL) {
  Clusters2process <- list.files(file.path(inputfolder, 'Clusters'), pattern = "clusters.RData", full.names = TRUE)
  Quals2process <- list.files(file.path(inputfolder, 'RawData'), pattern = "dataQuals.RData", full.names = TRUE)
  
  allClusters <- list()
  for (i in 1:length(Clusters2process)) {  
    data.file <- get(load(Clusters2process[i]))
    data.qual <- get(load(Quals2process[i]))
    fileID <- basename(Clusters2process[i])
    message("Processing file: ",fileID)
    
    #Sort data quals according to PB order in clusters
    SSlib.perPB <- data.qual$SSlib.perPB
    SSlib.perPB <- SSlib.perPB[match(rownames(data.file$soft.pVal), SSlib.perPB$PBreadNames),]
    pb.minLib <- SSlib.perPB$counts
    
    ##check accuracy only for autosomes and sex chrmosomes
    mask <- which(grepl('^chr[0-9X][0-9]?$', data.file$PBchrom))
    
    #get clusters IDs corresponding to a given chromosome
    chr.rows <- data.file$PBchrom[mask]
    chr.flag <- data.file$PBflag[mask]
    prob.tab <- data.file$soft.pVal[mask,]
    pb.minLib <- pb.minLib[mask]
    
    #filter out duplicates
    mask <- which(chr.flag == 16 | chr.flag == 0) 
    chr.rows <- chr.rows[mask]
    chr.flag <- chr.flag[mask]
    prob.tab <- prob.tab[mask,]
    pb.minLib <- pb.minLib[mask]
    
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
    
    #get clusters IDs corresponding to a given chromosome
    Clust.IDs <- getClusterIdentityPerChrPerDir(soft.clust=prob.tab, chr.rows=chr.rows, chr.flag=chr.flag)
    
    clust.acc.l <- list()
    for (prob.th in thresholds) {
      message("    Set threshold: ", prob.th)
      max.prob <- apply(prob.tab, 1, max)
      mask <- max.prob >= prob.th
      Clust.locations <- apply(prob.tab[mask,], 1, which.max) 
      
      #calculate clustering accuracy in comparison to expected values
      clust.acc <- Clust.locations == Clust.IDs[mask]
      acc.th <- table(clust.acc)
      
      clust.acc.l[[1+length(clust.acc.l)]] <- c(prob.th=prob.th, acc.th.match=unname(acc.th[2]), acc.th.sum=sum(acc.th), allReads=length(chr.rows))  
    }
    allClusters[[fileID]] <- as.data.frame( do.call(rbind, clust.acc.l) )
  } 
  #sum all counts over all data frames (per position)
  clust.acc.df <- Reduce("+", allClusters)
  
  #calcualte accuracy percentages
  clust.acc.df$prob.th <- thresholds
  clust.acc.df$th.acc <- clust.acc.df$acc.th.match / clust.acc.df$acc.th.sum
  clust.acc.df$th.clustReads <- clust.acc.df$acc.th.sum / clust.acc.df$allReads
  
  acc.plt <- ggplot(clust.acc.df) + geom_point(aes(x=th.acc, y=th.clustReads), color="deepskyblue4", size=10) + geom_linerange(aes(ymin=-Inf, x=th.acc, ymax=th.clustReads),color="deepskyblue4") + scale_y_continuous(limits = c(0.8, 1)) + scale_y_continuous(limits = c(0,1)) + ylab("(%) evaluated PB reads") + xlab("(%) correctly assigned PB reads") + geom_text(aes(x=th.acc, y=th.clustReads), label=c('all', thresholds[-1]), color="white") + theme_bw()
  message("DONE!!!")
  return(list(acc.plot=acc.plt, plot.table=clust.acc.df))
} 

inputfolder <- "/media/daewoooo/WORK/Clustering_project/WholeGenomeAnalysis/SaaRclust_results_DataQualFilt/"
thresholds <- c(0, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.99)
minLibs <- c(0,5,10,15,20,30)
ClustersAccuracyPerChrPerDir_minTh_minLib(inputfolder = inputfolder, thresholds=thresholds, minLibs=minLibs) -> accpltminLibminTh.obj

ClustersAccuracyPerChrPerDir_minTh_minLib <- function(inputfolder=NULL, thresholds=NULL, minLibs=NULL) {
  Clusters2process <- list.files(file.path(inputfolder, 'Clusters'), pattern = "clusters.RData", full.names = TRUE)
  Quals2process <- list.files(file.path(inputfolder, 'RawData'), pattern = "dataQuals.RData", full.names = TRUE)
  
  allClusters <- list()
  for (i in 1:length(Clusters2process)) {  
    data.file <- get(load(Clusters2process[i]))
    data.qual <- get(load(Quals2process[i]))
    fileID <- basename(Clusters2process[i])
    message("Processing file: ",fileID)
    
    #Sort data quals according to PB order in clusters
    SSlib.perPB <- data.qual$SSlib.perPB
    SSlib.perPB <- SSlib.perPB[match(rownames(data.file$soft.pVal), SSlib.perPB$PBreadNames),]
    pb.minLib <- SSlib.perPB$counts
    
    ##check accuracy only for autosomes and sex chrmosomes
    mask <- which(grepl('^chr[0-9X][0-9]?$', data.file$PBchrom))
    
    #get clusters IDs corresponding to a given chromosome
    chr.rows <- data.file$PBchrom[mask]
    chr.flag <- data.file$PBflag[mask]
    prob.tab <- data.file$soft.pVal[mask,]
    pb.minLib <- pb.minLib[mask]
    
    #filter out duplicates
    mask <- which(chr.flag == 16 | chr.flag == 0) 
    chr.rows <- chr.rows[mask]
    chr.flag <- chr.flag[mask]
    prob.tab <- prob.tab[mask,]
    pb.minLib <- pb.minLib[mask]
    
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
      message(" Set minLib: ", minLib)
      #Remove PB reads represneted by SSlib less than minLib
      filt <- pb.minLib >= minLib
      prob.tab.sub <- prob.tab[filt,]
      chr.rows.sub <- chr.rows[filt]
      chr.flag.sub <- chr.flag[filt]
      
      #get clusters IDs corresponding to a given chromosome
      Clust.IDs <- getClusterIdentityPerChrPerDir(soft.clust=prob.tab.sub, chr.rows=chr.rows.sub, chr.flag=chr.flag.sub)
      
      clust.acc.l <- list()
      for (prob.th in thresholds) {
        message("    Set threshold: ", prob.th)
        max.prob <- apply(prob.tab.sub, 1, max)
        mask <- max.prob >= prob.th
        Clust.locations <- apply(prob.tab.sub[mask,], 1, which.max) 
        
        #calculate clustering accuracy in comparison to expected values
        clust.acc <- Clust.locations == Clust.IDs[mask]
        acc.th <- table(clust.acc)
        
        clust.acc.l[[1+length(clust.acc.l)]] <- c(prob.th=prob.th, acc.th.match=unname(acc.th[2]), acc.th.sum=sum(acc.th), allReads=total.PBreads, minLib=minLib)  
      }
      clust.acc.minLib.l[[1+length(clust.acc.minLib.l)]] <- as.data.frame( do.call(rbind, clust.acc.l) )
    } 
  allClusters[[fileID]] <- do.call(rbind, clust.acc.minLib.l)
  }
  
  #sum all counts over all data frames (per position)
  clust.acc.df <- Reduce("+", allClusters)

    
  #calcualte accuracy percentages
  clust.acc.df$prob.th <- rep(thresholds, length(minLibs))
  clust.acc.df$minLib <- factor( rep(minLibs, each=length(thresholds)) )
  clust.acc.df$th.acc <- clust.acc.df$acc.th.match / clust.acc.df$acc.th.sum
  clust.acc.df$th.clustReads <- clust.acc.df$acc.th.sum / clust.acc.df$allReads
  
  acc.plt <- ggplot(clust.acc.df) + geom_point(aes(x=th.acc, y=th.clustReads, color=minLib, group=minLib)) + geom_linerange(aes(ymin=-Inf, x=th.acc, ymax=th.clustReads, color=minLib, group=minLib)) + scale_color_manual(values=brewer.pal(9,"Set1"))
  
  message("DONE!!!")
  return(list(acc.plot=acc.plt, plot.table=clust.acc.df))
} 


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

total.cov.bases <- sum(unlist(cov.bases.all)) 
#get genome size
library("biovizBase")
hg38Ideogram <- getIdeogram("hg38", cytoband = FALSE)
hg38Ideogram <- keepSeqlevels(hg38Ideogram, paste0('chr', c(1:22,'X')))
genome.size <- sum(as.numeric(seqlengths(hg38Ideogram)))
total.cov.bases/genome.size



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