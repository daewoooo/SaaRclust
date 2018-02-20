
#' Wrapper function to run saarclust pipeline for a given number of PB reads.
#'
#' @param minimap.file A path to the minimap file to load.
#' @param outputfolder A folder name to export to results.
#' @param num.clusters Expected number of clusters. (for 22 autosomes == 44 clusters)
#' @param minLib Minimal number of StrandS libraries being represent per long PB read
#' @param upperQ ...
#' @param EM.iter Number of iteration to run EM for.
#' @param store.counts Logical if to store raw read counts per PB read
#' @param HC.input Filaname where hard clustering results are stored
#' @inheritParams countProb
#' @export
#' @author David Porubsky

#minimap.file <- "/media/daewoooo/WORK/Clustering_project/WholeGenomeAnalysis/subsetMinimap.txt"
#minimap.file <- "/media/daewoooo/WORK/SS2PacBio_alignment_HG00733/Test_cluster_chr21&chr22/Minimap_out/SS2Pacbio_minimap_HG00733_k13_w1_L70_f0.01_Chr21andChr22_allSSreads"
#minimap.file <- "/media/daewoooo/WORK/Clustering_project/WholeGenomeAnalysis/NA12878_WashU_PBreads_chunk14.maf.gz"

SaaRclust <- function(minimap.file=NULL, outputfolder='SaaRclust_results', num.clusters=47, EM.iter=100, alpha=0.1, minLib=10, upperQ=0.95, theta.param=NULL, pi.param=NULL, logL.th=1, theta.constrain=FALSE, store.counts=FALSE, HC.input=NULL) {

  #get file ID
  fileID <- basename(minimap.file)
  fileID <- strsplit(fileID, "\\.")[[1]][1]
  
  #get directories for export
  rawdata.store <- file.path(outputfolder, 'RawData')
  Clusters.store <- file.path(outputfolder, 'Clusters')
  plots.store <- file.path(outputfolder, 'Plots')
  trashbin.store <- file.path(outputfolder, 'TrashBin')
  
  #Load Hard clustering results and initialize parameters of EM algorithm [temporary solution for snakemake]
  #destination <- file.path(Clusters.store, HC.input)
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
  
  ### Read in minimap output file ###
  suppressWarnings( tab.in <- importData(infile = minimap.file, removeDuplicates = TRUE) )
  #suppressWarnings( tab.in <- importTestData(infile = minimap.file, removeDuplicates = TRUE) ) #SLOW because test data have to be processed differently
  #suppressWarnings( tab.in <- importOldTestData(infile = minimap.file, removeDuplicates = TRUE) ) #use this function to import old test data (HG00733) from HGSVC
  #tab.in <- tab.in[tab.in$SSchrom != 'chrUn' & tab.in$SSchrom != 'chrX',] #applies only for test data
  #tab.in <- tab.in[tab.in$PBchrom %in% paste0('chr', c(18:22)),] #run only sertain chromosomes
  
  ### get some quality measures on imported data ### [OPTIONAL]
  data.qual.measures <- getQualMeasure(tab.in)
  destination <- file.path(rawdata.store, paste0(fileID, "_dataQuals.RData"))
  save(file = destination, data.qual.measures)
  
  ### Filter imported data ###
  tab.filt.l <- filterInput(inputData=tab.in, quantileSSreads = c(0, upperQ), minSSlibs = c(minLib,Inf))
  tab.filt <- tab.filt.l$tab.filt
  
  ### Store upperQ reads in trashBin ###
  upperQ.tab <- tab.in[tab.in$PBreadNames %in% tab.filt.l$upperQ.reads,]
  #upperQ.tab$PBreadNames <- factor(upperQ.tab$PBreadNames, levels=unique(upperQ.tab$PBreadNames))
  #tab.upperQ.l <- split(upperQ.tab, upperQ.tab$SSlibNames)
  #counts.upperQ.l <- countDirectionalReads(tab.upperQ.l)
  ptm <- startTimedMessage("Writing upperQ reads into a file")
  destination <- file.path(trashbin.store, paste0(fileID, "_upperQreads.gz"))
  gzf = gzfile(destination, 'w')
  write.table(x = upperQ.tab, file = gzf, quote = F, row.names = F)
  close(gzf)
  stopTimedMessage(ptm)
  #data.table::fwrite(upperQ.tab, destination)
  #gzip(destination)
  
  #take a smaller chunk of PB reads to process [NOT USED!!!]
  #tab.filt <- tab.filt[sample(nrow(tab.filt)),] #shuffle rows in tab
  #chunk <- unique(tab.filt$PBreadNames)[1:chunk.size]
  #tab.filt <- tab.filt[tab.filt$PBreadNames %in% chunk,]
  
  ### Sorting filtered data table by direction and chromosome ###
  ptm <- startTimedMessage("Sorting data")
  #additional sort by direction
  tab.filt <- tab.filt[order(tab.filt$PBflag),]
  
  #use PB read names as factor
  tab.filt <- tab.filt[order(tab.filt$PBchrom),]
  tab.filt$PBreadNames <- factor(tab.filt$PBreadNames, levels=unique(tab.filt$PBreadNames))
  
  #get PB chrom names from the ordered PB reads
  chr.l <- split(tab.filt$PBchrom, tab.filt$PBreadNames)
  chr.rows <- sapply(chr.l, function(x) x[1])
  
  #get PB directionality from the ordered PB reads
  pb.flag <- split(tab.filt$PBflag, tab.filt$PBreadNames)
  pb.flag <- sapply(pb.flag, unique)
  
  #get PB read position
  #pb.pos <- split(tab.filt$PBpos, tab.filt$PBreadNames)
  #pb.pos <- sapply(pb.pos, unique)
  
  #get PB read length
  pb.readLen <- split(tab.filt$PBreadLen, tab.filt$PBreadNames)
  pb.readLen <- sapply(pb.readLen, unique)
  
  #split data by SS library
  tab.l <- split(tab.filt, tab.filt$SSlibNames)
  stopTimedMessage(ptm)
  
  #### Count directional reads ###
  counts.l <- countDirectionalReads(tab.l)
  
  if (store.counts) {
    destination <- file.path(rawdata.store, paste0(fileID, "_counts.RData"))
    save(file = destination, counts.l)
  }
  
  ### Hard clustering ### [Excluded from this function]
  #hard.clust <- hardClust(counts.l, num.clusters=num.clusters, alpha=alpha)
  
  #get accuracy of hard clustering [OPTIONAL]
  #chr.clusts <- split(chr.rows, hard.clust$clust.id)
  #clust.acc <- getClusterAcc(chr.clusts)
  
  #initialize thetas
  #theta.l <- hard.clust$theta.estim
  #theta.l <- randomTheta(num.cells=100, num.clusters=num.clusters)
  
  #estimate pi based on # of PB reads in each cluster
  #readsPerCluts <- table(hard.clust$clust.id)
  #pi.param <- readsPerCluts/sum(readsPerCluts)
  
  ### EM algorithm ###
  soft.clust.obj <- EMclust(counts.l, theta.param=theta.param, pi.param=pi.param, num.iter=EM.iter, alpha=alpha, logL.th=logL.th)
  #soft.clust.obj2 <- EMclust(tab.l, theta.param=soft.clust.obj$theta.param, pi.param=soft.clust.obj$pi.param, num.iter=EM.iter, alpha=alpha) #to redo clustering
  
  #rescale theta ans run one more iteration to do soft clustering [EXPERIMENTAL]
  if (theta.constrain) {
    theta.expected <- num.clusters * c(0.25,0.25,0.5)
    theta.rescaled <- thetaRescale(theta.param=soft.clust.obj$theta.param, theta.expected=theta.expected)
    soft.clust.obj <- EMclust(counts.l, theta.param=theta.rescaled, pi.param=soft.clust.obj$pi.param, num.iter=1, alpha=alpha, logL.th=logL.th)
  }

  ### Merge clusters ### Part of hard clustering now!!!
  #get splitted clusters with same directionality and coming from the same chromosome [EXPERIMENTAL]
  #split.clust <- findSplitedClusters(theta.param=soft.clust$theta.param)
  #merged.pVals <- mergeSplitedClusters(cluster2merge=split.clust, soft.pVal=soft.clust$soft.pVal)
  #soft.clust$soft.pVal <- soft.clust$soft.pVal[,-unlist(split.clust)]
  #soft.clust$soft.pVal <- cbind(soft.clust$soft.pVal, merged.pVals)
  
  #Get pairs of clusters coming from the same chromosome but differs in directionality of PB reads  [EXPERIMENTAL]
  #clust.order <- findClusterPartners(theta.param=soft.clust.obj$theta.param)

  ### Save final results ###
  #add known chromosome and directionality of PB reads to a final data object
  soft.clust.obj$PBchrom <- as.character(chr.rows)
  soft.clust.obj$PBflag <- as.character(pb.flag)
  soft.clust.obj$pb.readLen <- as.numeric(unlist(pb.readLen))
  #export data 
  destination <- file.path(Clusters.store, paste0(fileID, "_clusters.RData"))
  save(file = destination, soft.clust.obj)

  return(soft.clust.obj)  #add cluster order??? [TODO]
}
