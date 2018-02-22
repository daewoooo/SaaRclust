#' Wrapper function to run saarclust pipeline.
#'
#' @param inputfolder A folder name where minimap files are stored.
#' @param store.bestAlign Store best alignements in RData object.
#' @param HC.only Perform only hard clustering and skip the rest of the pipeline.
#' @param numAlignments ...
#' @param verbose ... 
#' @inheritParams SaaRclust
#' @export
#' @author David Porubsky, Maryam Ghareghani

#inputfolder <- "/media/daewoooo/WORK/Clustering_project/WholeGenomeAnalysis/"

#load the function below into R if you want to run all steps in one command

runSaaRclust <- function(inputfolder=NULL, outputfolder="SaaRclust_results", num.clusters=54, EM.iter=100, alpha=0.01, minLib=10, logL.th=1, theta.constrain=FALSE, store.counts=FALSE, store.bestAlign=TRUE, numAlignments=30000, HC.only=TRUE, verbose=TRUE) {
  
  #=========================#
  ### Create directiories ###
  #=========================#
  
  #Create a master output directory
  outputfolder.destination <- file.path(inputfolder, outputfolder)
  if (!file.exists( outputfolder.destination)) {
    dir.create( outputfolder.destination)
  }
  
  #Directory to store raw read counts and best alignments
  if (store.counts) {
    rawdata.store <- file.path(outputfolder.destination, 'RawData')
    if (!file.exists(rawdata.store)) {
      dir.create(rawdata.store)
    }
  }
  
  if (store.bestAlign) {
    rawdata.store <- file.path(outputfolder.destination, 'RawData')
    if (!file.exists(rawdata.store)) {
      dir.create(rawdata.store)
    }
  }

  #Directory to store processed/clustered data
  Clusters.store <- file.path(outputfolder.destination, 'Clusters')
  if (!file.exists(Clusters.store)) {
    dir.create(Clusters.store)
  }
  
  #Directory to store plots
  plots.store <- file.path(outputfolder.destination, 'Plots')
  if (!file.exists(plots.store)) {
    dir.create(plots.store)
  }

  #Directory to store 'difficult' PacBio reads for later processing [TODO]
  trashbin.store <- file.path(outputfolder.destination, 'TrashBin')
  if (!file.exists(trashbin.store)) {
    dir.create(trashbin.store)
  }
  
  #Consider to have separate pipeline for Hard clustering [TODO]
  #numAlignments <- 50000 #perhaps add this parameter into a main function definition???
  #Load Hard clustering results if they were already created
  destination <- file.path(Clusters.store, paste0("hardClusteringResults_", as.integer(numAlignments),".RData"))
  if (!file.exists(destination)) {
    message("Hard clustering results not available!!!")
    message("Running Hard clustering")
    
    ### Get representative alignments to estimate theta and pi values ###
    destination <- file.path(rawdata.store, paste0("representativeAligns_", as.integer(numAlignments),".RData"))
    #reuse existing data if they were already created and save in a given location
    if (!file.exists(destination)) {
      best.alignments <- getRepresentativeAlignments(inputfolder=inputfolder, numAlignments=numAlignments, quantileSSreads=c(0,0.9), minSSlibs=c(35,Inf))
      if (store.bestAlign) {
        save(file = destination, best.alignments)
      }
    } else {
      best.alignments <- get(load(destination))
    }  
    
    #use PB read names as factor in order to export counts for every PB read (also for zero counts)
    best.alignments$PBreadNames <- factor(best.alignments$PBreadNames, levels=unique(best.alignments$PBreadNames))
    
    #split data by Strand-seq library
    tab.l <- split(best.alignments, best.alignments$SSlibNames)
    
    ### Count directional reads ###
    counts.l <- countDirectionalReads(tab.l)
    
    ### Perform k-means hard clustering method ###
    set.seed(1000) #in order to reproduce hard clustering results
    hardClust.ord <- hardClust(counts.l, num.clusters=num.clusters, nstart = 100)
    
    ### computing the accuracy of the hard clustering before merging lusters ### [OPTIONAL]
    #get PB chrom names from the ordered PB reads
    chr.l <- split(best.alignments$PBchrom, best.alignments$PBreadNames)
    chr.rows <- sapply(chr.l, function(x) x[1])
    #get PB directionality from the ordered PB reads
    pb.flag <- split(best.alignments$PBflag, best.alignments$PBreadNames)
    pb.flag <- sapply(pb.flag, unique)
    
    #get hard clustering accuracy
    #chr.clusts <- split(chr.rows, hardClust.ord)
    #clust.acc <- getClusterAcc(chr.clusts)
    #Create Hard clustering log
    log.destination <- file.path(outputfolder.destination, "hardClust.log")
    beQuiet <- file.create(log.destination)
    
    acc <- hardClustAccuracy(hard.clust = hardClust.ord, pb.chr = chr.rows, pb.flag = pb.flag, tab.filt = best.alignments)
    #print to log file
    write("Hard clustering summary:", file=log.destination, append=TRUE)
    write(paste("Accuracy before merging ", acc$acc), file=log.destination, append=TRUE)
    write(paste("Number of missing clusters =", length(acc$missed.clusters)), file=log.destination, append=TRUE)
    
    #Estimate theta parameter
    theta.estim <- estimateTheta(counts.l, ord=hardClust.ord, alpha=alpha)
    
    #Merge splitted clusters after hard clustering
    hardClust.ord.merged <- mergeClusters(kmeans.clust=hardClust.ord, theta.l=theta.estim, k=47)
    #findSplitedClusters(theta.param = theta.estim) -> to.join
    #hardClust.ord.merged <- hardClust.ord
    #for (i in 1:length(to.join)) {
    #  to.merge <- as.numeric(to.join[[i]])
    #  hardClust.ord.merged[ hardClust.ord.merged %in% to.merge] <- to.merge[1]     
    #}
  
    #Computing the accuracy of the hard clustering after merging
    acc <- hardClustAccuracy(hard.clust = hardClust.ord.merged, pb.chr = chr.rows, pb.flag = pb.flag, tab.filt = best.alignments)
    #print to log file
    write(paste("\nAccuracy after merging ", acc$acc), file=log.destination, append=TRUE)
    write(paste("Number of missing clusters =", length(acc$missed.clusters)), file=log.destination, append=TRUE)
    
    #Re-estimate theta parameter after cluster merging
    theta.estim <- estimateTheta(counts.l, ord=hardClust.ord.merged, alpha=alpha)
    
    #Initialize theta parameter
    theta.param <- theta.estim
    #Estimate pi parameter based on # of PB reads in each cluster
    readsPerCluts <- table(hardClust.ord.merged)
    pi.param <- readsPerCluts/sum(readsPerCluts)
    
    #save hard clustering results into a file
    hard.clust <- list(ord=hardClust.ord.merged, theta.param=theta.param, pi.param=pi.param)
    destination <- file.path(Clusters.store, paste0("hardClusteringResults_", as.integer(numAlignments), ".RData"))
    if (!file.exists(destination)) {
      save(file = destination, hard.clust)
    }  
    
  } else {
    message("Loading Hard clustering results")
    hard.clust <- get(load(destination))
  }
  
  if(!HC.only) {
    #Initialize theta parameter
    theta.param <- hard.clust$theta.param
    #Initialize pi parameter
    pi.param <- hard.clust$pi.param 
    
    #List files to process
    file.list <- list.files(path = inputfolder, pattern = "chunk.+maf", full.names = TRUE)
  
    ### Main loop to process all files using EM algorithm ###
    for (file in file.list) {
      if (verbose) {
        clust.obj <- SaaRclust(minimap.file=file, outputfolder=outputfolder.destination, num.clusters=length(pi.param), EM.iter=EM.iter, alpha=alpha, minLib=minLib, theta.param=theta.param, pi.param=pi.param, logL.th=logL.th, theta.constrain=theta.constrain)
      } else {
        suppressMessages(  clust.obj <- SaaRclust(minimap.file=file, outputfolder=outputfolder.destination, num.clusters=num.clusters, EM.iter=EM.iter, alpha=alpha, theta.param=theta.param, pi.param=pi.param, logL.th=logL.th, theta.constrain=theta.constrain) )
      }
    }
  
    ### Evaluate soft clustering accuracy ###
    #get clusters IDs corresponding to a given chromosome
    Clust.IDs <- getClusterIdentity(soft.clust=clust.obj$EM.data$soft.pVal, chr.rows=clust.obj$Data2plot$PBchrom)
  
    #get best clusters for each PB read given set threshold (reports multiple location when max pVal is lower than prob.th)
    thresholds <- c(0.5, 0.6, 0.7, 0.8, 0.9)
    clust.acc.l <- list()
    chr.rows <- clust.obj$Data2plot$PBchrom
    
    for (prob.th in thresholds) {
      Best.clusters <- exportGenomicLocations(soft.clust=clust.obj$EM.data$soft.pVal, prob.th)
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
      both.th.acc <- unname( both.th[2]/sum(both.th) )
      above.th.acc <- unname( above.th[2]/sum(above.th) )
      below.th.acc <- unname( below.th[2]/sum(below.th) )
      both.th.clustReads <- sum(both.th)/length(chr.rows)
      above.th.clustReads <- sum(above.th)/length(chr.rows)
      below.th.clustReads <- sum(below.th)/length(chr.rows)
      
      clust.acc.l[[1+length(clust.acc.l)]] <- c(prob.th=prob.th, both.th.acc=both.th.acc, above.th.acc=above.th.acc, below.th.acc=below.th.acc, both.th.clustReads=both.th.clustReads, above.th.clustReads=above.th.clustReads, below.th.clustReads=below.th.clustReads)  
    }
    clust.acc.df <- as.data.frame( do.call(rbind, clust.acc.l) )
  
    #TMP code
    #plt <- ggplot(clust.acc.df) + geom_point(aes(x=above.th.acc, y=above.th.clustReads, color=as.character(prob.th)), size=5) + geom_linerange(aes(ymin=-Inf, x=above.th.acc, ymax=above.th.clustReads, color=as.character(prob.th))) + scale_y_continuous(limits = c(0,1)) + scale_color_manual(values = brewer.pal(n=7, name="Set1"), name="Prob threshold") + ylab("(%) evaluated PB reads") + xlab("(%) correctly assigned PB reads")
    #plt <- plt + geom_point(data=clust.acc.df, aes(x=below.th.acc, y=below.th.clustReads, color=as.character(prob.th)), size=5) + geom_linerange(aes(ymin=-Inf, x=below.th.acc, ymax=below.th.clustReads, color=as.character(prob.th))) + scale_y_continuous(limits = c(0,1)) + scale_color_manual(values = brewer.pal(n=7, name="Set1"), name="Prob threshold")
    #plt + geom_point(data=clust.acc.df, aes(x=both.th.acc, y=both.th.clustReads, color=as.character(prob.th)), size=5) + geom_linerange(aes(ymin=-Inf, x=both.th.acc, ymax=both.th.clustReads, color=as.character(prob.th))) + scale_y_continuous(limits = c(0,1)) + scale_color_manual(values = brewer.pal(n=7, name="Set1"), name="Prob threshold")
  
    plt <- ggplot(clust.acc.df) + geom_point(aes(x=above.th.acc, y=above.th.clustReads), color="red", size=10) + geom_linerange(aes(ymin=-Inf, x=above.th.acc, ymax=above.th.clustReads),color="red") + scale_y_continuous(limits = c(0,1)) + ylab("(%) evaluated PB reads") + xlab("(%) correctly assigned PB reads") + geom_text(aes(x=above.th.acc, y=above.th.clustReads), label=thresholds, color="white")
    plt <- plt + geom_point(data=clust.acc.df, aes(x=below.th.acc, y=below.th.clustReads), color="blue", size=10) + geom_linerange(aes(ymin=-Inf, x=below.th.acc, ymax=below.th.clustReads), color="blue") + geom_text(aes(x=below.th.acc, y=below.th.clustReads), label=thresholds, color="white")
    plt + geom_point(data=clust.acc.df, aes(x=both.th.acc, y=both.th.clustReads), color="green", size=10) + geom_linerange(aes(ymin=-Inf, x=both.th.acc, ymax=both.th.clustReads), color="green") + geom_text(aes(x=both.th.acc, y=both.th.clustReads), label=thresholds, color="white")
  
  } else {
    return(hard.clust)
  }

}

