inputfolder <- "/media/daewoooo/WORK/Clustering_project/WholeGenomeAnalysis/SaaRclust_results/Clusters/"
processClusters(inputfolder = inputfolder) -> data.obj

processClusters <- function(inputfolder=NULL) {
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
    chr.rows <- data.file$PBchrom
    
    #get clusters IDs corresponding to a given chromosome
    Clust.IDs <- getClusterIdentity(soft.clust=data.file$soft.pVal, chr.rows=chr.rows)
    
    clust.acc.l <- list()
    for (prob.th in thresholds) {
      Best.clusters <- exportGenomicLocations(soft.clust=data.file$soft.pVal, prob.th)
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
  