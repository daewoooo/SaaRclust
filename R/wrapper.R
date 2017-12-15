#' Wrapper function to run saarclust pipeline.
#'
#' @param inputfolder A folder name where minimap files are stored.
#' @param verbose ... 
#' @inheritParams SaaRclust
#' @export
#' @author David Porubsky


#minimap test files are present in TestData folder
#minimap.file <- "/media/daewoooo/WORK/SS2PacBio_alignment_HG00733/WholeGenomeAnalysis/subsetMinimap.txt"
#minimap.file <- "/media/daewoooo/WORK/SS2PacBio_alignment_HG00733/Test_cluster_chr21&chr22/Minimap_out/SS2Pacbio_minimap_HG00733_k13_w1_L70_f0.01_Chr21andChr22_allSSreads"
#minimap.file <- "/media/daewoooo/WORK/Clustering_project/WholeGenomeAnalysis/NA12878_WashU_PBreads_chunk14.maf.gz"

#load the function below into R if you want to run all steps in one command

runSaaRclust <- function(inputfolder=NULL, outputfolder="./SaaRclust_results", num.clusters=44, EM.iter=100, alpha=0.1, logL.th=1, theta.constrain=FALSE, store.counts=FALSE, verbose=TRUE) {
  
  #=========================#
  ### Create directiories ###
  #=========================#
  
  #Create a master output directory
  if (!file.exists(outputfolder)) {
    dir.create(outputfolder)
  }
  
  #Directory to store raw read counts
  if (store.counts) {
    rawcounts.store <- file.path(outputfolder, 'RawCounts')
    if (!file.exists(rawcounts.store)) {
      dir.create(rawcounts.store)
    }
  }

  #Directory to stare processed/clustered data
  Clusters.store <- file.path(outputfolder, 'Clusters')
  if (!file.exists(Clusters.store)) {
    dir.create(Clusters.store)
  }

  #Directory to store 'difficult' PacBio reads for later processing [TODO]
  #trashbin.store <- file.path(outputfolder, 'TrashBin')
  #if (!file.exists(trashbin.store)) {
  #  dir.create(trashbin.store)
  #}
  
  #List of files to porecess
  file.list <- list.files(path = inputfolder, pattern = "mimimapChunk")
  
  for (file in file.list) {
    if (verbose) {
      clust.obj <- SaaRclust(minimap.file=file, outputfolder=outputfolder, num.clusters=num.clusters, EM.iter=EM.iter, alpha=alpha, logL.th=logL.th, theta.constrain=theta.constrain)
    } else {
      suppressMessages(  clust.obj <- SaaRclust(minimap.file=file, outputfolder=outputfolder, num.clusters=num.clusters, EM.iter=EM.iter, alpha=alpha, logL.th=logL.th, theta.constrain=theta.constrain) )
    }
  }
  
  #load processed data [TODO]
  
}

