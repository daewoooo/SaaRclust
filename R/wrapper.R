#' Wrapper function to run saarclust pipeline.
#'
#' @param inputfolder A folder name where minimap files are stored.
#' @param verbose ... 
#' @inheritParams SaaRclust
#' @export
#' @author David Porubsky

#inputfolder <- "/media/daewoooo/WORK/Clustering_project/WholeGenomeAnalysis/"

#load the function below into R if you want to run all steps in one command

runSaaRclust <- function(inputfolder=NULL, outputfolder="./SaaRclust_results", num.clusters=44, EM.iter=100, alpha=0.1, logL.th=1, theta.constrain=FALSE, store.counts=FALSE, store.bestAlign=TRUE, verbose=TRUE) {
  
  #=========================#
  ### Create directiories ###
  #=========================#
  
  #Create a master output directory
  if (!file.exists(outputfolder)) {
    dir.create(outputfolder)
  }
  
  #Directory to store raw read counts and best alignments
  if (store.counts) {
    rawdata.store <- file.path(outputfolder, 'RawData')
    if (!file.exists(rawdata.store)) {
      dir.create(rawdata.store)
    }
  }
  
  if (store.bestAlign) {
    rawdata.store <- file.path(outputfolder, 'RawData')
    if (!file.exists(rawdata.store)) {
      dir.create(rawdata.store)
    }
  }

  #Directory to store processed/clustered data
  Clusters.store <- file.path(outputfolder, 'Clusters')
  if (!file.exists(Clusters.store)) {
    dir.create(Clusters.store)
  }

  #Directory to store 'difficult' PacBio reads for later processing [TODO]
  #trashbin.store <- file.path(outputfolder, 'TrashBin')
  #if (!file.exists(trashbin.store)) {
  #  dir.create(trashbin.store)
  #}
  
  ### Get representative alignments to estimate theta and pi values ###
  numAlignments <- 30000 #perhaps add this parameter into a main function definition???
  destination <- file.path(rawdata.store, paste0("representativeAligns_",numAlignments,".RData"))
  #reuse existing data if they were already created and save in a given location
  if (!file.exists(destination)) {
    best.alignments <- getRepresentativeAlignments(inputfolder=inputfolder, numAlignments=numAlignments)
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
  
  ### Perform hard clustering method ###
  hardClust.ord <- hardClust(counts.l, num.clusters=num.clusters)
  
  #Estimate theta parameter
  theta.estim <- estimateTheta(counts.l, ord=hardClust.ord, alpha=alpha)
  
  #Merge splitted clusters after hard clustering
  hardClust.ord.merged <- mergeClusters(kmeans.clust=hardClust.ord, theta.l=theta.estim, k = 46)
  
  #Re-estimate theta parameter after cluster merging
  theta.estim <- estimateTheta(counts.l, ord=hardClust.ord.merged, alpha=alpha)
  
  ### Get Hard clustering accuracy ###
  #get PB chrom names from the ordered PB reads
  chr.l <- split(best.alignments$PBchrom, best.alignments$PBreadNames)
  chr.rows <- sapply(chr.l, unique)
  #get PB directionality from the ordered PB reads
  pb.flag <- split(best.alignments$PBflag, best.alignments$PBreadNames)
  pb.flag <- sapply(pb.flag, unique)
  # define the classe (true clusters) -- we may need to add an additional line for removing the PB reads with more than 1 chr or direction
  classes <- paste0(chr.rows, "_", pb.flag)
  names(classes) <- names(chr.rows)
  # accuracy based on chromosome location and directionality
  acc <- maryam_hardClustAccuracy(hard.clust = hard.clust, classes=classes, tab.filt = tab.filt)
  # accuracy based on chrom only
  #acc_chrom <- maryam_hardClustAccuracy(hard.clust = hard.clust.chrom, classes=chr.rows, tab.filt = tab.filt)
  print(acc)
  print(paste("number of missing clusters =", length(acc$missed.clusters)))
  
  #Initialize theta parameter
  theta.param <- theta.estim
  #Estimate pi parameter based on # of PB reads in each cluster
  readsPerCluts <- table(hardClust.ord.merged)
  pi.param <- readsPerCluts/sum(readsPerCluts)
  
  #List files to process
  file.list <- list.files(path = inputfolder, pattern = "chunk.+maf", full.names = TRUE)
  
  ### Main loop to process all files using EM algorithm ###
  for (file in file.list) {
    if (verbose) {
      clust.obj <- SaaRclust(minimap.file=file, outputfolder=outputfolder, num.clusters=num.clusters, EM.iter=EM.iter, alpha=alpha, theta.param=theta.param, pi.param=pi.param, logL.th=logL.th, theta.constrain=theta.constrain)
    } else {
      suppressMessages(  clust.obj <- SaaRclust(minimap.file=file, outputfolder=outputfolder, num.clusters=num.clusters, EM.iter=EM.iter, alpha=alpha, theta.param=theta.param, pi.param=pi.param, logL.th=logL.th, theta.constrain=theta.constrain) )
    }
  }
  
  #load processed data [TODO]
  
}

