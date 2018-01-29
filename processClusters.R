inputfolder <- "/media/daewoooo/WORK/Clustering_project/WholeGenomeAnalysis/SaaRclust_results/Clusters/"
ClustersAccuracyPerChr(inputfolder = inputfolder) -> data.obj
ClustersAccuracyPerChrPerDir(inputfolder = inputfolder) -> data.obj
accuracyRanking(inputfolder = inputfolder) -> data.obj

### Evaluate soft clustering accuracy ###
#get best clusters for each PB read given set threshold (reports multiple location when max pVal is lower than prob.th)
ClustersAccuracyPerChr <- function(inputfolder=NULL) {
  files2process <- list.files(inputfolder, pattern = "clusters.RData", full.names = TRUE)
  
  ptm <- startTimedMessage("Processing clusters")
  allClusters <- list()
  for (file in files2process) {
    data.file <- get(load(file))
    fileID <- basename(file)
      
    thresholds <- c(0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)
    
    #Get accuracy only for autosomes and sex chrmosomes
    mask <- which(grepl('^chr[0-9X][0-9]?$', data.file$PBchrom))
    chr.rows <- data.file$PBchrom[mask]
    chr.flag <- data.file$PBflag[mask]
    prob.tab <- data.file$soft.pVal[mask,]
    #Filter out duplicates
    mask <- which(chr.flag == 16 | chr.flag == 0) 
    chr.rows <- chr.rows[mask]
    chr.flag <- chr.flag[mask]
    prob.tab <- prob.tab[mask,]
    
    #Find WC cluster in all cells
    theta.sums <- Reduce("+", data.file$theta.param)
    remove.clust <- which.max(theta.sums[,3])
    #Remove probabilities for always WC cluster
    prob.tab <- prob.tab[,-remove.clust]
    
    #remove reads with highest probability for WC cluster ???
    #max.clust <- apply(prob.tab, 1, which.max)
    #mask <- max.clust != remove.clust
    #chr.rows <- chr.rows[mask]
    #chr.flag <- chr.flag[mask]
    #prob.tab <- prob.tab[mask,]
    
    #Get clusters IDs corresponding to a given chromosome
    Clust.IDs <- getClusterIdentityPerChr(soft.clust= prob.tab, chr.rows=chr.rows)
  
    clust.acc.l <- list()
    for (prob.th in thresholds) {
      Best.clusters <- exportGenomicLocations(soft.clust=prob.tab, prob.th)
      #Best.clusters <- exportGenomicLocationsAll(soft.clust=prob.tab, prob.th)
      Clust.locations <- Best.clusters$clust.IDs
      #Clust.locations <- Best.clusters
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

#This function compares accuracy for all possible max probs using 'exportGenomicLocationsAll' function
ClustersAccuracyPerChr2 <- function(inputfolder=NULL) {
  files2process <- list.files(inputfolder, pattern = "clusters.RData", full.names = TRUE)
  
  ptm <- startTimedMessage("Processing clusters")
  allClusters <- list()
  for (file in files2process) {
    data.file <- get(load(file))
    fileID <- basename(file)
    
    thresholds <- c(0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)
    
    #Get accuracy only for autosomes and sex chrmosomes
    mask <- which(grepl('^chr[0-9X][0-9]?$', data.file$PBchrom))
    chr.rows <- data.file$PBchrom[mask]
    chr.flag <- data.file$PBflag[mask]
    prob.tab <- data.file$soft.pVal[mask,]
    #Filter out duplicates
    mask <- which(chr.flag == 16 | chr.flag == 0) 
    chr.rows <- chr.rows[mask]
    chr.flag <- chr.flag[mask]
    prob.tab <- prob.tab[mask,]
    
    #Find WC cluster in all cells
    theta.sums <- Reduce("+", data.file$theta.param)
    remove.clust <- which.max(theta.sums[,3])
    #Remove probabilities for always WC cluster
    prob.tab <- prob.tab[,-remove.clust]
    
    #remove reads with highest probability for WC cluster ???
    #max.clust <- apply(prob.tab, 1, which.max)
    #mask <- max.clust != remove.clust
    #chr.rows <- chr.rows[mask]
    #chr.flag <- chr.flag[mask]
    #prob.tab <- prob.tab[mask,]
    
    #Get clusters IDs corresponding to a given chromosome
    Clust.IDs <- getClusterIdentityPerChr(soft.clust= prob.tab, chr.rows=chr.rows)
    
    clust.acc.l <- list()
    for (prob.th in thresholds) {
      Best.clusters <- exportGenomicLocationsAll(soft.clust=prob.tab, prob.th)
      Clust.locations <- Best.clusters
      names(Clust.locations) <- chr.rows
      
      #replicate known clust IDs per chromosome
      Clust.IDs.expand <- rep(Clust.IDs , as.numeric(table(chr.rows)))
      
      #compare if best clust IDs correspond to known clust IDs for any given PB read
      clust.acc <- sapply(1:length(Clust.locations), function(x) any(Clust.locations[[x]] %in% Clust.IDs.expand[[x]]))

      max.prob <- apply(prob.tab, 1, max)
      mask <- max.prob >= prob.th
      above.th <- table(clust.acc[mask])
      
      clust.acc.l[[1+length(clust.acc.l)]] <- c(prob.th=prob.th, above.th.match=unname(above.th[2]), above.th.sum=sum(above.th), allReads=length(chr.rows))  
    }
    allClusters[[fileID]] <- as.data.frame( do.call(rbind, clust.acc.l) )
  } 
  #sum all counts over all data frames (per position)
  clust.acc.df <- Reduce("+", allClusters)
  
  #calcualte accuracy percentages
  clust.acc.df$prob.th <- thresholds
  clust.acc.df$above.th.acc <- clust.acc.df$above.th.match / clust.acc.df$above.th.sum
  clust.acc.df$above.th.clustReads <- clust.acc.df$above.th.sum / clust.acc.df$allReads
  
  #plotting  
  plt <- ggplot(clust.acc.df) + geom_point(aes(x=above.th.acc, y=above.th.clustReads), color="red", size=10) + geom_linerange(aes(ymin=-Inf, x=above.th.acc, ymax=above.th.clustReads),color="red") + scale_y_continuous(limits = c(0,1)) + ylab("(%) evaluated PB reads") + xlab("(%) correctly assigned PB reads") + geom_text(aes(x=above.th.acc, y=above.th.clustReads), label=thresholds, color="white")
 
  stopTimedMessage(ptm)
  return(list(plot=plt, plot.table=clust.acc.df))
}  
  

thresholds <- c(0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)
inputfolder <- "/media/daewoooo/WORK/Clustering_project/WholeGenomeAnalysis/SaaRclust_results_DataQualFilt/Clusters/"
ClustersAccuracyHighestTwoDist(inputfolder=inputfolder, thresholds=thresholds) -> HighestTwoDist_accplt.obj

ClustersAccuracyHighestTwoDist <- function(inputfolder=NULL, thresholds=thresholds) {
  files2process <- list.files(inputfolder, pattern = "clusters.RData", full.names = TRUE)
  
  ptm <- startTimedMessage("Processing clusters")
  allClusters <- list()
  for (file in files2process) {
    data.file <- get(load(file))
    fileID <- basename(file)
    
    #Get accuracy only for autosomes and sex chrmosomes
    mask <- which(grepl('^chr[0-9X][0-9]?$', data.file$PBchrom))
    chr.rows <- data.file$PBchrom[mask]
    chr.flag <- data.file$PBflag[mask]
    prob.tab <- data.file$soft.pVal[mask,]
    #Filter out duplicates
    mask <- which(chr.flag == 16 | chr.flag == 0) 
    chr.rows <- chr.rows[mask]
    chr.flag <- chr.flag[mask]
    prob.tab <- prob.tab[mask,]
    
    #Find WC cluster in all cells
    theta.sums <- Reduce("+", data.file$theta.param)
    remove.clust <- which.max(theta.sums[,3])
    #Remove probabilities for always WC cluster
    prob.tab <- prob.tab[,-remove.clust]
    
    #Get clusters IDs corresponding to a given chromosome
    #Clust.IDs <- getClusterIdentityPerChr(soft.clust= prob.tab, chr.rows=chr.rows)
    #replicate known clust IDs per chromosome
    #Clust.IDs.expand <- rep(Clust.IDs , as.numeric(table(chr.rows)))
    Clust.IDs.expand <- getClusterIdentityPerChrPerDir(soft.clust=prob.tab, chr.rows=chr.rows, chr.flag=chr.flag)
    
    clust.acc.l <- list()
    for (prob.th in thresholds) {
      message("Processing threshold ...", prob.th)
      #get distance between two highest probs
      max2Dist <- apply(prob.tab, 1, function(x) abs(diff(sort(x, decreasing = T)[1:2])))
      
      #filter distance between two highest probs based on set threshold
      mask <- max2Dist >= prob.th
      Clust.IDs.expand.sub <- Clust.IDs.expand[mask]
      prob.tab.sub <- prob.tab[mask,]
      chr.rows.sub <- chr.rows[mask]
      
      #get max prob for filtered data
      Clust.locations <- apply(prob.tab.sub, 1, which.max)
      
      #compare if best clust IDs correspond to known clust IDs for any given PB read
      clust.acc <- sapply(1:length(Clust.locations), function(x) any(Clust.locations[x] %in% Clust.IDs.expand.sub[[x]]))
      
      above.th <- table(clust.acc)
      
      clust.acc.l[[1+length(clust.acc.l)]] <- c(prob.th=prob.th, above.th.match=unname(above.th[2]), above.th.sum=sum(above.th), allReads=length(chr.rows))  
    }
    allClusters[[fileID]] <- as.data.frame( do.call(rbind, clust.acc.l) )
  } 
  #sum all counts over all data frames (per position)
  clust.acc.df <- Reduce("+", allClusters)
  
  #calcualte accuracy percentages
  clust.acc.df$prob.th <- thresholds
  clust.acc.df$above.th.acc <- clust.acc.df$above.th.match / clust.acc.df$above.th.sum
  clust.acc.df$above.th.clustReads <- clust.acc.df$above.th.sum / clust.acc.df$allReads
  
  #plotting  
  plt <- ggplot(clust.acc.df) + geom_point(aes(x=above.th.acc, y=above.th.clustReads), color="deepskyblue4", size=10) + geom_linerange(aes(ymin=-Inf, x=above.th.acc, ymax=above.th.clustReads),color="deepskyblue4") + scale_y_continuous(limits = c(0,1)) + ylab("(%) evaluated PB reads") + xlab("(%) correctly assigned PB reads") + geom_text(aes(x=above.th.acc, y=above.th.clustReads), label=thresholds, color="white")
  
  stopTimedMessage(ptm)
  return(list(plot=plt, plot.table=clust.acc.df))
}  


inputfolder <- "/media/daewoooo/WORK/Clustering_project/WholeGenomeAnalysis/SaaRclust_results_DataQualFilt/"
#Set required probability thresholds
thresholds <- c(0, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)
ClustersAccuracyPerChrPerDir(inputfolder=inputfolder, thresholds=thresholds) -> accplt.obj
destination <- file.path(inputfolder, "accPlot_perChrperDir.RData") 
save(file = destination, accplt.obj)

ClustersAccuracyPerChrPerDir <- function(inputfolder=NULL, thresholds=thresholds) {
  files2process <- list.files(inputfolder, pattern = "clusters.RData", full.names = TRUE)
  
  ptm <- startTimedMessage("Processing clusters")
  allClusters <- list()
  for (file in files2process) {
    data.file <- get(load(file))
    fileID <- basename(file)
    
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
    
    #get clusters IDs corresponding to a given chromosome
    Clust.IDs <- getClusterIdentityPerChrPerDir(soft.clust=prob.tab, chr.rows=chr.rows, chr.flag=chr.flag)
    
    clust.acc.l <- list()
    for (prob.th in thresholds) {
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
  
  acc.plt <- ggplot(clust.acc.df) + geom_point(aes(x=th.acc, y=th.clustReads), color="deepskyblue4", size=10) + geom_linerange(aes(ymin=-Inf, x=th.acc, ymax=th.clustReads),color="deepskyblue4") + scale_y_continuous(limits = c(0,1)) + ylab("(%) evaluated PB reads") + xlab("(%) correctly assigned PB reads") + geom_text(aes(x=th.acc, y=th.clustReads), label=c('all', thresholds[-1]), color="white") + theme_bw()
  
  stopTimedMessage(ptm)
  return(list(acc.plot=acc.plt, plot.table=clust.acc.df))
}  


inputfolder <- "/media/daewoooo/WORK/Clustering_project/WholeGenomeAnalysis/SaaRclust_results_DataQualFilt/"
thresholds <- c(0, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)
ClustersAccuracyPerChrPerDir2(inputfolder=inputfolder, thresholds=thresholds) -> accplt.obj

ClustersAccuracyPerChrPerDir2 <- function(inputfolder=NULL, thresholds=thresholds, minLib=15) {
  Clusters2process <- list.files(file.path(inputfolder, 'Clusters'), pattern = "clusters.RData", full.names = TRUE)
  Quals2process <- list.files(file.path(inputfolder, 'RawData'), pattern = "dataQuals.RData", full.names = TRUE)
  
  ptm <- startTimedMessage("Processing clusters")
  allClusters <- list()
  for (i in 1:length(Clusters2process)) {  
    #data.file <- get(load(file))
    #fileID <- basename(file)
    data.file <- get(load(Clusters2process[i]))
    data.qual <- get(load(Quals2process[i]))
    fileID <- basename(Clusters2process[i])
    
    #Sort data quals according to PB order in clusters
    SSlib.perPB <- data.qual$SSlib.perPB
    #mask.PBnames <- SSlib.perPB$PBreadNames[SSlib.perPB$counts >= minLib]
    SSlib.perPB <- SSlib.perPB[match(rownames(data.file$soft.pVal), SSlib.perPB$PBreadNames),]
    pb.minLib <- SSlib.perPB$counts
    
    #Remove PB reads represneted by SSlib less than minLib
    #mask <- SSlib.perPB$counts >= minLib
    #chr.rows <- data.file$PBchrom[mask]
    #chr.flag <- data.file$PBflag[mask]
    #prob.tab <- data.file$soft.pVal[mask,]
    #SSlib.perPB <- SSlib.perPB[mask,]
    
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
    #Clust.IDs <- getClusterIdentityPerChr(soft.clust=prob.tab, chr.rows=chr.rows)
    Clust.IDs <- getClusterIdentityPerChrPerDir(soft.clust=prob.tab, chr.rows=chr.rows, chr.flag=chr.flag)
    
    clust.acc.l <- list()
    for (prob.th in thresholds) {
      max.prob <- apply(prob.tab, 1, max)
      mask <- max.prob >= prob.th
      Clust.locations <- apply(prob.tab[mask,], 1, which.max) 
      
      #replicate known clust IDs per chromosome
      #Clust.IDs.expand <- rep(Clust.IDs , as.numeric(table(chr.rows[mask])))
      
      #compare if best clust IDs correspond to known clust IDs for any given PB read
      #clust.acc <- sapply(1:length(Clust.locations), function(x) any(Clust.locations[[x]] %in% Clust.IDs.expand[[x]]))
      
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
  
  acc.plt <- ggplot(clust.acc.df) + geom_point(aes(x=th.acc, y=th.clustReads), color="deepskyblue4", size=10) + geom_linerange(aes(ymin=-Inf, x=th.acc, ymax=th.clustReads),color="deepskyblue4") + scale_y_continuous(limits = c(0,1)) + ylab("(%) evaluated PB reads") + xlab("(%) correctly assigned PB reads") + geom_text(aes(x=th.acc, y=th.clustReads), label=c('all', thresholds[-1]), color="white") + theme_bw()
  
  stopTimedMessage(ptm)
  return(list(acc.plot=acc.plt, plot.table=clust.acc.df))
}  


inputfolder <- "/media/daewoooo/WORK/Clustering_project/WholeGenomeAnalysis/SaaRclust_results_DataQualFilt/"
boxplotDistSSlibsPerPB <- function(inputfolder=NULL, thresholds=0) {
  Clusters2process <- list.files(file.path(inputfolder, 'Clusters'), pattern = "clusters.RData", full.names = TRUE)
  Quals2process <- list.files(file.path(inputfolder, 'RawData'), pattern = "dataQuals.RData", full.names = TRUE)
  
  ptm <- startTimedMessage("Processing clusters")
  allSSlib.perPB <- list()
  for (i in 1:length(Clusters2process)) {  
    data.file <- get(load(Clusters2process[i]))
    data.qual <- get(load(Quals2process[i]))
    PB.read.len <- data.table::fread(paste0('gunzip -cq ', PBlen2process[i]), header=T, verbose = F, showProgress = F)
    fileID <- basename(Clusters2process[i])
    
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
  
  stopTimedMessage(ptm)
  return(list(acc.plot=box.plt, plot.table=allSSlib.perPB.df))
} 



inputfolder <- "/media/daewoooo/WORK/Clustering_project/WholeGenomeAnalysis/SaaRclust_results_DataQualFilt/"
boxplotDistPBreadLen <- function(inputfolder=NULL, thresholds=0) {
  Clusters2process <- list.files(file.path(inputfolder, 'Clusters'), pattern = "clusters.RData", full.names = TRUE)
  Quals2process <- list.files(file.path(inputfolder, 'RawData'), pattern = "dataQuals.RData", full.names = TRUE)
  PBlen2process <- list.files("/media/daewoooo/WORK/Clustering_project/WholeGenomeAnalysis/PBreadLen/", pattern = "gz", full.names = TRUE)
  
  ptm <- startTimedMessage("Processing clusters")

  all.PBreadLen <- list()
  for (i in 1:length(Clusters2process)) {  
    data.file <- get(load(Clusters2process[i]))
    data.qual <- get(load(Quals2process[i]))
    PB.read.len <- data.table::fread(paste0('gunzip -cq ', PBlen2process[i]), header=T, verbose = F, showProgress = F)
    fileID <- basename(Clusters2process[i])
    
    #Select required PB read lengths
    PB.read.len <- PB.read.len[match(rownames(data.file$soft.pVal), PB.read.len$PBreadNames),]
    
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
    PB.read.len <- PB.read.len[mask,]
    
    #Find WC cluster in all cells
    theta.sums <- Reduce("+", data.file$theta.param)
    remove.clust <- which.max(theta.sums[,3])
    #Remove probabilities for always WC cluster
    prob.tab <- prob.tab[,-remove.clust]
    
    #get clusters IDs corresponding to a given chromosome
    Clust.IDs <- getClusterIdentityPerChrPerDir(soft.clust=prob.tab, chr.rows=chr.rows, chr.flag=chr.flag)
    
    PB.read.len.dist <- list()
    for (prob.th in thresholds) {
      max.prob <- apply(prob.tab, 1, max)
      mask <- max.prob >= prob.th
      Clust.locations <- apply(prob.tab[mask,], 1, which.max)  
      
      #calculate clustering accuracy in comparison to expected values
      clust.acc <- Clust.locations == Clust.IDs[mask]
      acc.th <- table(clust.acc)
      
      #Split counts of SSlib per PB by accuracy vector (correct vs incorrect PB read assignemnts)
      PB.read.len.dist[[1+length(PB.read.len.dist)]] <- split(PB.read.len$PBreadLen[mask], clust.acc)
    }
    unlist( sapply(PB.read.len.dist, function(x) x['TRUE']) ) -> trues
    unlist( sapply(PB.read.len.dist, function(x) x['FALSE']) )-> falses
    all.PBreadLen[[fileID]] <- list(trues=trues, falses=falses)
  } 
  unlist( sapply(all.PBreadLen, function(x) x['trues']), use.names = F ) -> trues
  unlist( sapply(all.PBreadLen, function(x) x['falses']), use.names = F )-> falses
  ID <- rep(c('Correct', 'Incorrect'), c(length(trues), length(falses)))
  counts <- c(trues, falses)
  all.PBreadLen.df <- data.frame(counts=counts, ID=ID)
  #Plot 
  box.plt <- ggplot(all.PBreadLen.df, aes(x=ID, y=counts, fill=ID)) + geom_boxplot(outlier.colour="red") + scale_fill_manual(values = c("darkolivegreen3" ,"darkgoldenrod1"), guide="none") + xlab("") + ylab("PacBio read length") + theme_bw()
  
  stopTimedMessage(ptm)
  return(list(acc.plot=box.plt, plot.table=all.PBreadLen.df))
} 

inputfolder <- "/media/daewoooo/WORK/Clustering_project/WholeGenomeAnalysis/SaaRclust_results_DataQualFilt/Clusters/"
accuracyRanking(inputfolder = inputfolder ) -> plt.obj
destination <- file.path(inputfolder, "rankingPlot.RData") 
save(file = destination, plt.obj)

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
  ranking.plt <- ggplot(table.ranks.df, aes(x=rank.acc, y=Freq)) + geom_bar(fill="chartreuse4", stat="identity") + xlab("Probability ranking of true cluster") + ylab("Frequency") + theme_bw() + scale_y_continuous(labels=comma)
  
  #get accuracy measure (sum probabilities of true cluster divided by all PB reads)
  all.prob.trueClust.v <- unlist(all.prob.trueClust)
  overal.acc <- sum(all.prob.trueClust.v)/length(all.prob.trueClust.v)
  
  return(list(ranking.plt=ranking.plt, ranking.table=table.ranks.df,  overal.acc=overal.acc))
} 


### Heatmap plotting ###
#prepare data
data.file <- get(load("/media/daewoooo/WORK/Clustering_project/WholeGenomeAnalysis/SaaRclust_results_DataQualFilt/Clusters/NA12878_WashU_PBreads_chunk00_clusters.RData"))
#soft.clust.df <- as.data.frame(data.file$soft.pVal)
soft.clust.df <- as.data.frame(data.file$soft.pVal[,c(10,43,19,2,34,12)])
soft.clust.df$PBreadNames <- rownames(data.file$soft.pVal)
soft.clust.df$PBchrom <- data.file$PBchrom
soft.clust.df$PBflag <- data.file$PBflag

#take only first 5 chromosomes
soft.clust.df <- soft.clust.df[soft.clust.df$PBchrom %in% paste0('chr', 1:3),]

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

### Plot theta parameter ###
theta.plt <- plotThetaEstimates(theta.param=data.file$theta.param)

### Plot pi parameter ###
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

 