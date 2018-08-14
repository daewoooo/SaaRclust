

#Export corresponding cluster for true chromosome and directionality

getClusterIdentityPerChrPerDir <- function(soft.clust, chr.rows, chr.flag) {
  max.Clust <- apply(soft.clust, 1, which.max)
  unique.clust.ID <- paste0(chr.rows,"_",chr.flag)
  unique.clust.ID <- factor(unique.clust.ID, unique(unique.clust.ID))
  clustByChromByflag <- split(max.Clust, unique.clust.ID)
  clustIdPerChrom <- lapply(clustByChromByflag, function(x) names(which.max(table(x))))
  clustIDperPB <- rep(as.numeric(unlist(clustIdPerChrom)), table(unique.clust.ID))
  return(clustIDperPB)
}

############################################################################################################################################
#This function calculates clustering accuracy over different probability cutoffs.
#Depth of coverage is reported as a cumulative numeber of bases sequenced (sum of PB read lenghts/genome size)

ClustersAccuracyPerChrPerDir <- function(Clusters2process=NULL, Quals2process=NULL, thresholds=1:10/10, minLib=0) {
	allClusters <- list()
	fileID <- basename(Clusters2process)
	message("Processing file: ",fileID)

	#load required data
	data.file <- get(load(Clusters2process))
	data.qual <- get(load(Quals2process))
	pb.readLen <- data.file$pb.readLen

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
	return(as.data.frame(do.call(rbind, clust.acc.l)))
} 


############################################################################################################################################
#This function calculates distribution of minimal SSlibs represented per PBread in sets of correctly and incorrectly assigned PB reads.

boxplotDistSSlibsPerPB <- function(Clusters2process=NULL, Quals2process=NULL, thresholds=1:10/10) {
  fileID <- basename(Clusters2process)
  message("Processing file: ",fileID)
  
  data.file <- get(load(Clusters2process))
  data.qual <- get(load(Quals2process))
  
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
  return(list(trues=trues, falses=falses))
} 


############################################################################################################################################
#This function calculates probability ranking of true clusters. 
#Rank1 = max probability being the true cluster
#Rank2 = second best probability being the true cluster

accuracyRanking <- function(inputfile=NULL) {
  fileID <- basename(inputfile)
  message("Processing file: ",fileID)
  
  #load required data
  data.file <- get(load(inputfile))
  
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

	return(list(ranks=ranks.store, trueClust=prob.trueClust))
}



VennDiagramStats <- function(inputfolder=NULL, thresholds=c(0.5,0.75)) {
  Clusters2process <- list.files(file.path(inputfolder, 'Clusters'), pattern = "clusters.RData", full.names = TRUE)
  #Quals2process <- list.files(file.path(inputfolder, 'RawData'), pattern = "dataQuals.RData", full.names = TRUE)
  
  venn.stats <- matrix(0L, nrow=3+length(thresholds), ncol=3+length(thresholds))
  rownames(venn.stats) <- c("mapped_normal_chroms", "mapped_unknown_chroms", "unmapped", paste0("threshold",thresholds))
  colnames(venn.stats) <- rownames(venn.stats)
  
  allClusters <- list()
  for (i in 1:length(Clusters2process)) {
    fileID <- basename(Clusters2process[i])
    message("Processing file: ",fileID)
    
    #load required data
    data.file <- get(load(Clusters2process[i]))
    #data.qual <- get(load(Quals2process[i]))
    #pb.readLen <- data.file$pb.readLen
    
    # getting total number of reads and three different sets of PB reads
    num.pb.reads <- length(data.file$PBchrom)
    
    mapped.normal.chroms <- which(grepl('^chr[0-9X][0-9]?$', data.file$PBchrom))
    mapped.unknown.chroms <- setdiff(1:num.pb.reads, mapped.normal.chroms)
    unmapped <- which(data.file$PBflag == 4)
    
    # update intersection matrix
    venn.stats[1,1] = venn.stats[1,1] + length(mapped.normal.chroms)
    venn.stats[2,2] = venn.stats[2,2] + length(mapped.unknown.chroms)
    venn.stats[3,3] = venn.stats[3,3] + length(unmapped)
    
    prob.tab <- data.file$soft.pVal
    
    #Find WC cluster in all cells
    theta.sums <- Reduce("+", data.file$theta.param)
    remove.clust <- which.max(theta.sums[,3])
    #Remove probabilities for always WC cluster
    prob.tab <- prob.tab[,-remove.clust]
    
    filt.pb.reads <- list()
    for (j in 1:length(thresholds)) {
      message("    Set threshold: ", thresholds[j])
      max.prob <- apply(prob.tab, 1, max)
      filt.pb.reads[[j]] <- which(max.prob >= thresholds[j])
      
      # update intersection matrix
      venn.stats[j+3,j+3] <- venn.stats[j+3,j+3]+length(filt.pb.reads[[j]])
    }
    
    # update intersection matrix
    for (j in 1:length(thresholds))
    {
      venn.stats[1,j+3] <- venn.stats[1,j+3] + length(intersect(mapped.normal.chroms, filt.pb.reads[[j]]))
      venn.stats[2,j+3] <- venn.stats[2,j+3] + length(intersect(mapped.unknown.chroms, filt.pb.reads[[j]]))
      venn.stats[3,j+3] <- venn.stats[3,j+3] + length(intersect(unmapped, filt.pb.reads[[j]]))
    }
  }
  
  return(venn.stats)
} 
                       
                       
getProbDiff <- function(inputfolder=NULL) {
  Clusters2process <- list.files(file.path(inputfolder, 'Clusters'), pattern = "clusters.RData", full.names = TRUE)
  #Quals2process <- list.files(file.path(inputfolder, 'RawData'), pattern = "dataQuals.RData", full.names = TRUE)
  
  diff.probs <- NULL
  for (i in 1:length(Clusters2process)) {
    fileID <- basename(Clusters2process[i])
    message("Processing file: ",fileID)
    
    #load required data
    data.file <- get(load(Clusters2process[i]))
    
    prob.tab <- data.file$soft.pVal
    
    #Find WC cluster in all cells
    theta.sums <- Reduce("+", data.file$theta.param)
    remove.clust <- which.max(theta.sums[,3])
    #Remove probabilities for always WC cluster
    prob.tab <- prob.tab[,-remove.clust]
    
    # Normalize prob.tab after removong the garnage cluster
    prob.tab <- prob.tab / rowSums(prob.tab)
    
    diff.probs <- c(diff.probs, sapply(1:nrow(prob.tab), function(i) (max(prob.tab[i,]) - sort(prob.tab[i,], decreasing = T)[2])))
  }
  
  return(diff.probs)
}
