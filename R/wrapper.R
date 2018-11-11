#' Wrapper function to run saarclust pipeline.
#'
#' @param inputfolder A folder name where minimap files are stored.
#' @param store.bestAlign Store best alignements in RData object.
#' @param HC.only Perform only hard clustering and skip the rest of the pipeline.
#' @param numAlignments Required number of best PBvsSS alignmnets to selest for hard clustering.
#' @param verbose Set to \code{TRUE} to print function messages.
#' @param HC.input Filaname where hard clustering results are stored
#' @param cellNum specifies the number of single cells to be used in clustering
#' @inheritParams SaaRclust
#' @inheritParams EMclust
#' @export
#' @author David Porubsky, Maryam Ghareghani

runSaaRclust <- function(inputfolder=NULL, outputfolder="SaaRclust_results", num.clusters=54, EM.iter=100, alpha=0.01, minLib=10, upperQ=0.95, logL.th=1, theta.constrain=FALSE, store.counts=FALSE, store.bestAlign=TRUE, numAlignments=30000, HC.only=TRUE, verbose=TRUE, cellNum=NULL, log.scale=FALSE) {
  
  #=========================#
  ### Create directiories ###
  #=========================#
  
  #Create a master output directory
  outputfolder.destination <- file.path(outputfolder)
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
  
  #Load Hard clustering results if they were already created
  destination <- file.path(Clusters.store, "hardClusteringResults.RData")
  if (!file.exists(destination)) {
    message("Hard clustering results not available!!!")
    message("Running Hard clustering")
    
    ### Get representative alignments to estimate theta and pi values ###
    destination <- file.path(rawdata.store, "representativeAligns.RData")
    #reuse existing data if they were already created and save in a given location
    if (!file.exists(destination)) {
      # setting the min number of SS libs as a cutoff to use PB reads in hard clustering
      cov.cutoff = 35
      if (!is.null(cellNum)) {
        cov.cutoff <- round(cellNum/4)
      }
      best.alignments <- getRepresentativeAlignments(inputfolder=inputfolder, numAlignments=numAlignments, quantileSSreads=c(0,0.9), minSSlibs=c(cov.cutoff,Inf))
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
    #counts.l <- importBams(bamfolder = bamfolder, chromosomes = chromosomes, bin.length = 1000000) TODO: implement as an option
    
    #subsetting single cell libraries
    if (!is.null(cellNum)) {
      counts.l = counts.l[1:cellNum]
    }
    
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
    
    #Create Hard clustering log
    log.destination <- file.path(outputfolder.destination, "hardClust.log")
    beQuiet <- file.create(log.destination)
    
    #get hard clustering accuracy
    acc <- hardClustAccuracy(hard.clust = hardClust.ord, pb.chr = chr.rows, pb.flag = pb.flag, tab.filt = best.alignments)
    #print to log file
    write("Hard clustering summary:", file=log.destination, append=TRUE)
    write(paste("Accuracy before merging ", acc$acc), file=log.destination, append=TRUE)
    write(paste("Number of missing clusters =", length(acc$missed.clusters)), file=log.destination, append=TRUE)
    
    #Estimate theta parameter
    theta.estim <- estimateTheta(counts.l, hard.clust=hardClust.ord, alpha=alpha)
    
    #Merge splitted clusters after hard clustering
    hardClust.ord.merged <- mergeClusters(hard.clust=hardClust.ord, theta.l=theta.estim, k=47)
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
    theta.estim <- estimateTheta(counts.l, hard.clust=hardClust.ord.merged, alpha=alpha)
    
    #Initialize theta parameter
    theta.param <- theta.estim
    #Estimate pi parameter based on # of PB reads in each cluster
    readsPerCluts <- table(hardClust.ord.merged)
    pi.param <- readsPerCluts/sum(readsPerCluts)
    
    #save hard clustering results into a file
    hard.clust <- list(ord=hardClust.ord.merged, theta.param=theta.param, pi.param=pi.param)
    destination <- file.path(Clusters.store, "hardClusteringResults.RData")
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
    file.list <- list.files(path = inputfolder, pattern = "chunk.+maf.gz$", full.names = TRUE)
  
    ### Main loop to process all files using EM algorithm ###
    for (file in file.list) {
      if (verbose) {
        clust.obj <- SaaRclust(minimap.file=file, outputfolder=outputfolder.destination, num.clusters=length(pi.param), EM.iter=EM.iter, alpha=alpha, minLib=minLib, upperQ=upperQ, theta.param=theta.param, pi.param=pi.param, logL.th=logL.th, theta.constrain=theta.constrain, log.scale=log.scale)
      } else {
        suppressMessages(  clust.obj <- SaaRclust(minimap.file=file, outputfolder=outputfolder.destination, num.clusters=num.clusters, EM.iter=EM.iter, alpha=alpha, minLib=minLib, upperQ=upperQ, theta.param=theta.param, pi.param=pi.param, logL.th=logL.th, theta.constrain=theta.constrain, log.scale=log.scale) )
      }
    }
  
  } else {
    return(hard.clust)
  }

}
