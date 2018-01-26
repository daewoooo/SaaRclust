inputfolder <- "/media/daewoooo/WORK/Clustering_project/WholeGenomeAnalysis/SaaRclust_results/Clusters/"
ClustersAccuracyPerChr(inputfolder = inputfolder) -> data.obj
ClustersAccuracyPerChrPerDir(inputfolder = inputfolder) -> data.obj
accuracyRanking(inputfolder = inputfolder) -> data.obj

ClustersAccuracyPerChr <- function(inputfolder=NULL) {
  files2process <- list.files(inputfolder, pattern = "clusters.RData", full.names = TRUE)
  
  ptm <- startTimedMessage("Processing clusters")
  allClusters <- list()
  for (file in files2process) {
    data.file <- get(load(file))
    fileID <- basename(file)
    
    ### Prepara data for plotting ###
    #soft.clust.df <- as.data.frame(data.file$soft.pVal)
    #soft.clust.df$PBreadNames <- rownames(data.file$soft.pVal)
    #soft.clust.df$PBchrom <- data.file$PBchrom
    #soft.clust.df$PBflag <- data.file$PBflag
      
    ### Evaluate soft clustering accuracy ###
    #get best clusters for each PB read given set threshold (reports multiple location when max pVal is lower than prob.th)
    thresholds <- c(0.5, 0.6, 0.7, 0.8, 0.9)
    
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
    
    #get clusters IDs corresponding to a given chromosome
    Clust.IDs <- getClusterIdentity_old(soft.clust= prob.tab, chr.rows=chr.rows)
  
    clust.acc.l <- list()
    for (prob.th in thresholds) {
      Best.clusters <- exportGenomicLocations(soft.clust=prob.tab, prob.th)
      Clust.locations <- Best.clusters$clust.IDs
      mask <- Best.clusters$th.boolean
      names(Clust.locations) <- chr.rows
        
      #replicate known clust IDs per chromosome
      Clust.IDs.expand <- rep(Clust.IDs , as.numeric(table(chr.rows)))
      
      #compare if best clust IDs correspond to known clust IDs for any given PB read
      clust.acc <- sapply(1:length(Clust.locations), function(x) any(Clust.locations[[x]] %in% Clust.IDs.expand[[x]]))
      both.th <- table(clust.acc)
      above.th <- table(clust.acc[mask])
      below.th <- table(clust.acc[!mask])
        
      clust.acc.l[[1+length(clust.acc.l)]] <- c(prob.th=prob.th, both.th.match=unname(both.th[2]), both.th.sum=sum(both.th), above.th.match=unname(above.th[2]), above.th.sum=sum(above.th), below.th.match=unname(below.th[2]), below.th.sum=sum(below.th), allReads=length(chr.rows))  
    }
    allClusters[[fileID]] <- as.data.frame( do.call(rbind, clust.acc.l) )
  } 
  #sum all counts over all data frames (per position)
  clust.acc.df <- Reduce("+", allClusters)
  
  #calcualte accuracy percentages
  clust.acc.df$prob.th <- thresholds
  clust.acc.df$both.th.acc <- clust.acc.df$both.th.match / clust.acc.df$both.th.sum
  clust.acc.df$above.th.acc <- clust.acc.df$above.th.match / clust.acc.df$above.th.sum
  clust.acc.df$below.th.acc <- clust.acc.df$below.th.match / clust.acc.df$below.th.sum
  clust.acc.df$both.th.clustReads <- clust.acc.df$both.th.sum / clust.acc.df$allReads
  clust.acc.df$above.th.clustReads <- clust.acc.df$above.th.sum / clust.acc.df$allReads
  clust.acc.df$below.th.clustReads <- clust.acc.df$below.th.sum / clust.acc.df$allReads
  
  #plotting  
  plt <- ggplot(clust.acc.df) + geom_point(aes(x=above.th.acc, y=above.th.clustReads), color="red", size=10) + geom_linerange(aes(ymin=-Inf, x=above.th.acc, ymax=above.th.clustReads),color="red") + scale_y_continuous(limits = c(0,1)) + ylab("(%) evaluated PB reads") + xlab("(%) correctly assigned PB reads") + geom_text(aes(x=above.th.acc, y=above.th.clustReads), label=thresholds, color="white")
  plt <- plt + geom_point(data=clust.acc.df, aes(x=below.th.acc, y=below.th.clustReads), color="blue", size=10) + geom_linerange(aes(ymin=-Inf, x=below.th.acc, ymax=below.th.clustReads), color="blue") + geom_text(aes(x=below.th.acc, y=below.th.clustReads), label=thresholds, color="white")
  plt <- plt + geom_point(data=clust.acc.df, aes(x=both.th.acc, y=both.th.clustReads), color="green", size=10) + geom_linerange(aes(ymin=-Inf, x=both.th.acc, ymax=both.th.clustReads), color="green") + geom_text(aes(x=both.th.acc, y=both.th.clustReads), label=thresholds, color="white")
  
  stopTimedMessage(ptm)
  return(list(plot=plt, plot.table=clust.acc.df))
}  
  

ClustersAccuracyPerChrPerDir <- function(inputfolder=NULL) {
  files2process <- list.files(inputfolder, pattern = "clusters.RData", full.names = TRUE)
  
  ptm <- startTimedMessage("Processing clusters")
  allClusters <- list()
  for (file in files2process) {
    data.file <- get(load(file))
    fileID <- basename(file)
    
    ### Prepara data for plotting ###
    #soft.clust.df <- as.data.frame(data.file$soft.pVal)
    #soft.clust.df$PBreadNames <- rownames(data.file$soft.pVal)
    #soft.clust.df$PBchrom <- data.file$PBchrom
    #soft.clust.df$PBflag <- data.file$PBflag
    
    ### Evaluate soft clustering accuracy ###
    #get best clusters for each PB read given set threshold (reports multiple location when max pVal is lower than prob.th)
    thresholds <- c(0, 0.5, 0.6, 0.7, 0.8, 0.9)
    
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
    
    #get clusters IDs corresponding to a given chromosome
    #Clust.IDs <- getClusterIdentity_old(soft.clust= prob.tab, chr.rows=chr.rows)
    Clust.IDs <- getClusterIdentity(soft.clust=prob.tab, chr.rows=chr.rows, chr.flag=chr.flag)
    
    clust.acc.l <- list()
    for (prob.th in thresholds) {
      #Best.clusters <- exportGenomicLocations(soft.clust=prob.tab, prob.th)
      #Clust.locations <- Best.clusters$clust.IDs
      #mask <- Best.clusters$th.boolean
      #names(Clust.locations) <- chr.rows
      
      #replicate known clust IDs per chromosome
      #Clust.IDs.expand <- rep(Clust.IDs , as.numeric(table(chr.rows)))
      
      max.prob <- apply(prob.tab, 1, max)
      mask <- max.prob >= prob.th
      Clust.locations <- apply(prob.tab[mask,], 1, which.max)  
      
      #compare if best clust IDs correspond to known clust IDs for any given PB read
      #clust.acc <- sapply(1:length(Clust.locations), function(x) any(Clust.locations[[x]] %in% Clust.IDs.expand[[x]]))
      clust.acc <- Clust.locations == Clust.IDs[mask]
      acc.th <- table(clust.acc)
      #above.th <- table(clust.acc[mask])
      #below.th <- table(clust.acc[!mask])
      
      #clust.acc.l[[1+length(clust.acc.l)]] <- c(prob.th=prob.th, both.th.match=unname(both.th[2]), both.th.sum=sum(both.th), above.th.match=unname(above.th[2]), above.th.sum=sum(above.th), below.th.match=unname(below.th[2]), below.th.sum=sum(below.th), allReads=length(chr.rows))  
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
  
  plt <- ggplot(clust.acc.df) + geom_point(aes(x=th.acc, y=th.clustReads), color="red", size=10) + geom_linerange(aes(ymin=-Inf, x=th.acc, ymax=th.clustReads),color="red") + scale_y_continuous(limits = c(0,1)) + ylab("(%) evaluated PB reads") + xlab("(%) correctly assigned PB reads") + geom_text(aes(x=th.acc, y=th.clustReads), label=c('all', thresholds[-1]), color="white")
  
  #calcualte accuracy percentages
  #clust.acc.df$prob.th <- thresholds
  #clust.acc.df$both.th.acc <- clust.acc.df$both.th.match / clust.acc.df$both.th.sum
  #clust.acc.df$above.th.acc <- clust.acc.df$above.th.match / clust.acc.df$above.th.sum
  #clust.acc.df$below.th.acc <- clust.acc.df$below.th.match / clust.acc.df$below.th.sum
  #clust.acc.df$both.th.clustReads <- clust.acc.df$both.th.sum / clust.acc.df$allReads
  #clust.acc.df$above.th.clustReads <- clust.acc.df$above.th.sum / clust.acc.df$allReads
  #clust.acc.df$below.th.clustReads <- clust.acc.df$below.th.sum / clust.acc.df$allReads
  
  #plotting  
  #plt <- ggplot(clust.acc.df) + geom_point(aes(x=above.th.acc, y=above.th.clustReads), color="red", size=10) + geom_linerange(aes(ymin=-Inf, x=above.th.acc, ymax=above.th.clustReads),color="red") + scale_y_continuous(limits = c(0,1)) + ylab("(%) evaluated PB reads") + xlab("(%) correctly assigned PB reads") + geom_text(aes(x=above.th.acc, y=above.th.clustReads), label=thresholds, color="white")
  #plt <- plt + geom_point(data=clust.acc.df, aes(x=below.th.acc, y=below.th.clustReads), color="blue", size=10) + geom_linerange(aes(ymin=-Inf, x=below.th.acc, ymax=below.th.clustReads), color="blue") + geom_text(aes(x=below.th.acc, y=below.th.clustReads), label=thresholds, color="white")
  #plt <- plt + geom_point(data=clust.acc.df, aes(x=both.th.acc, y=both.th.clustReads), color="green", size=10) + geom_linerange(aes(ymin=-Inf, x=both.th.acc, ymax=both.th.clustReads), color="green") + geom_text(aes(x=both.th.acc, y=both.th.clustReads), label=thresholds, color="white")
  
  stopTimedMessage(ptm)
  return(list(plot=plt, plot.table=clust.acc.df))
}  



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
    
    Clust.IDs <- getClusterIdentity(soft.clust=prob.tab, chr.rows=chr.rows, chr.flag=chr.flag)
    #replicate known clust IDs per chromosome
    #Clust.IDs.expand <- rep(Clust.IDs , as.numeric(table(chr.rows)))
    #Clust.IDs.matrix <- do.call(rbind, Clust.IDs.expand)
    
    #for each row of probability matrix select max probability corresponding to true cluster id
    #max.prob.trueClust <- sapply(1:nrow(prob.tab), function(x) max(prob.tab[x, Clust.IDs.matrix[x,]]))
    prob.trueClust <- sapply(1:nrow(prob.tab), function(x) prob.tab[x,Clust.IDs[x]])
    #having highest probability for the true cluster find the index in sorted probabilities (decreasing)
    #rank.acc <- sapply(1:nrow(prob.tab), function(x) which(sort(prob.tab[x,], decreasing = T) == max.prob.trueClust[x]))
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
  ranking.plt <- ggplot(table.ranks.df, aes(x=rank.acc, y=Freq)) + geom_bar(fill="red", stat="identity") + xlab("Probability ranking of true cluster") + ylab("Frequency")
  
  #get accuracy measure (sum probabilities of true cluster divided by all PB reads)
  all.prob.trueClust.v <- unlist(all.prob.trueClust)
  overal.acc <- sum(all.prob.trueClust.v)/length(all.prob.trueClust.v)
  
  return(list(ranking.plt=ranking.plt, ranking.table=table.ranks.df,  overal.acc=overal.acc))
} 


#heatmap plotting
#prepare data
data.file <- get(load("/media/daewoooo/WORK/Clustering_project/WholeGenomeAnalysis/SaaRclust_results_test/Clusters/NA12878_WashU_PBreads_chunk00_clusters.RData"))
soft.clust.df <- as.data.frame(data.file$soft.pVal)
soft.clust.df$PBreadNames <- rownames(data.file$soft.pVal)
soft.clust.df$PBchrom <- data.file$PBchrom
soft.clust.df$PBflag <- data.file$PBflag

#find WC cluster in all cells
theta.sums <- Reduce("+", data.file$theta.param)
remove.clust <- which.max(theta.sums[,3])

theta.modif <- lapply(data.file$theta.param, function(x) x[-remove.clust,])

#Get pairs of clusters coming from the same chromosome but differs in directionality of PB reads  [EXPERIMENTAL]
clust.order <- findClusterPartners(theta.param=theta.modif)
clust.order <- c(clust.order, length(clust.order)+1)

#swap column which is always WC
a <- soft.clust.df[,remove.clust]
b <- soft.clust.df[,47]
soft.clust.df[,remove.clust] <- b
soft.clust.df[,47] <- a

hm.plt <- plotHeatmap(pVal.df=soft.clust.df, colOrder=clust.order, num.clusters=47)

#plot theta parameter
theta.plt <- plotThetaEstimates(theta.param=data.file$theta.param)

#plot pi parameter
ord <- order(data.file$pi.param, decreasing = T)
Clust.IDs <- getClusterIdentity_old(soft.clust= data.file$soft.pVal, chr.rows=data.file$PBchrom)

ord <- do.call(c, Clust.IDs[paste0('chr', 1:22)])

library("biovizBase")
hg38Ideogram <- getIdeogram("hg38", cytoband = FALSE)
hg38Ideogram <- keepSeqlevels(hg38Ideogram, paste0('chr', 1:22))
hg38Ideogram <- sort(hg38Ideogram)
chr.len <- end(ranges(hg38Ideogram))
chr.len.df <- data.frame(chr.len=rep(chr.len, each=2))

chr.len.df$pi.param <- data.file$pi.param[ord]
chr.len.df$x <- factor(rownames(chr.len.df), levels=rownames(chr.len.df))

ggplot(chr.len.df) + geom_bar(aes(x=x, y=chr.len), fill="black", stat="identity") 

 