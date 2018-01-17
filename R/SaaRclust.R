
#' Wrapper function to run saarclust pipeline for a given number of PB reads.
#'
#' @param minimap.file A path to the minimap file to load.
#' @param outputfolder A folder name to export to results.
#' @param num.clusters Expected number of clusters. (for 22 autosomes == 44 clusters)
#' @param EM.iter Number of iteration to run EM for.
#' @param store.counts Logical if to store raw read counts per PB read
#' @inheritParams countProb
#' @export
#' @author David Porubsky

#minimap.file <- "/media/daewoooo/WORK/Clustering_project/WholeGenomeAnalysis/subsetMinimap.txt"
#minimap.file <- "/media/daewoooo/WORK/SS2PacBio_alignment_HG00733/Test_cluster_chr21&chr22/Minimap_out/SS2Pacbio_minimap_HG00733_k13_w1_L70_f0.01_Chr21andChr22_allSSreads"
#minimap.file <- "/media/daewoooo/WORK/Clustering_project/WholeGenomeAnalysis/NA12878_WashU_PBreads_chunk14.maf.gz"

SaaRclust <- function(minimap.file=NULL, outputfolder='SaaRclust_analysis', num.clusters=46, EM.iter=100, alpha=0.1, theta.param=theta.param, pi.param=pi.param, logL.th=1, theta.constrain=FALSE, store.counts=FALSE) {

  #get file ID
  fileID <- basename(minimap.file)
  fileID <- strsplit(fileID, "\\.")[[1]][1]
  
  #get directories for export
  rawdata.store <- file.path(outputfolder, 'RawData')
  Clusters.store <- file.path(outputfolder, 'Clusters')
  plots.store <- file.path(outputfolder, 'Plots')
  #trashbin.store <- file.path(outputfolder, 'TrashBin')
  
  ### Read in minimap output file ###
  suppressWarnings( tab.in <- importTestData(infile = minimap.file, removeDuplicates = TRUE) ) #SLOW because test data have to be processed differently
  #suppressWarnings( tab.in <- importOldTestData(infile = minimap.file, removeDuplicates = TRUE) ) #use this function to import old test data (HG00733) from HGSVC
  #tab.in <- tab.in[tab.in$SSchrom != 'chrUn' & tab.in$SSchrom != 'chrX',] #applies only for test data
  #tab.in <- tab.in[tab.in$PBchrom %in% paste0('chr', c(18:22)),] #run only sertain chromosomes
  
  #get some quality measures on imported data [OPTIONAL]
  #tab.in.quals <- getQualMeasure(tab.in) #time consuming!!!
  #qual.plt <- plotQualMeasure(tab.in.quals)
  
  ### Filter imported data ###
  tab.filt <- filterInput(inputData=tab.in, quantileSSreads = c(0, 0.9), minSSlibs = c(10,Inf))
  
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
  
  #split data by SS library
  tab.l <- split(tab.filt, tab.filt$SSlibNames)
  stopTimedMessage(ptm)
  
  #### Count directional reads ###
  counts.l <- countDirectionalReads(tab.l)
  
  if (store.counts) {
    destination <- file.path(rawdata.store, paste0(fileID, ".RData"))
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
    soft.clust.rescale <- SaaRclust(tab.l, theta.param=theta.rescaled, pi.param=pi.param, num.iter=1)
  }

  ### Merge clusters ### Part of hard clustering!!!
  #get splitted clusters with same directionality and coming from the same chromosome [EXPERIMENTAL]
  #split.clust <- findSplitedClusters(theta.param=soft.clust$theta.param)
  #merged.pVals <- mergeSplitedClusters(cluster2merge=split.clust, soft.pVal=soft.clust$soft.pVal)
  #soft.clust$soft.pVal <- soft.clust$soft.pVal[,-unlist(split.clust)]
  #soft.clust$soft.pVal <- cbind(soft.clust$soft.pVal, merged.pVals)
  
  #Get pairs of clusters coming from the same chromosome but differs in directionality of PB reads  [EXPERIMENTAL]
  clust.order <- findClusterPartners(theta.param=soft.clust.obj$theta.param)

  ### Prepara data for plotting ###
  soft.clust.df <- as.data.frame(soft.clust.obj$soft.pVal)
  soft.clust.df$PBreadNames <- levels(tab.filt$PBreadNames)
  soft.clust.df$PBchrom <- as.character(chr.rows)
  soft.clust.df$PBflag <- as.character(pb.flag)
  #soft.clust.df$hardClust <- hard.clust$clust.id

  ### Plotting data ###
  ptm <- startTimedMessage("Exporting plots")
  #plot likelihood function
  logl.df <- data.frame(log.l=soft.clust.obj$log.l, iter=1:length(soft.clust.obj$log.l))
  logl.diff <- round(max(logl.df$log.l) - min(logl.df$log.l), digits = 0)
  logL.plt <- ggplot(logl.df) + geom_line(aes(x=iter, y=log.l), color="red") + xlab("EM iterations") + ylab("Likelihood function") + annotate("text",  x=Inf, y = Inf, label = paste("Diff: ",logl.diff), vjust=1, hjust=1)
  #plot cluster accuracy
  #acc.plt <- plotClustAccuracy(pVal.df = soft.clust.df, num.clusters = num.clusters) Merge with Maryam's function
  #plot theta values
  theta.plt <- plotThetaEstimates(theta.param=soft.clust.obj$theta.param, title=fileID)
  #plot heatmap
  hm.plt <- plotHeatmap(pVal.df=soft.clust.df, colOrder=clust.order, num.clusters=num.clusters)

  #Save plots
  destination <- file.path(plots.store, paste0(fileID, "_logL.pdf"))
  ggsave(filename = destination, plot = logL.plt, width = 8, height = 5)
  destination <- file.path(plots.store, paste0(fileID, "_thetaEstim.pdf"))
  ggsave(filename = destination, plot = theta.plt, width = 20, height = 20)
  destination <- file.path(plots.store, paste0(fileID, "_heatmap.pdf"))
  pdf(destination, width = 15, height = 10) 
  hm.plt
  dev.off()
  stopTimedMessage(ptm)

  ### Save data ###
  destination <- file.path(Clusters.store, paste0(fileID, ".RData"))
  save(file = destination, soft.clust.obj)

  return(list(Data2plot=soft.clust.df, EM.data=soft.clust.obj))  #add cluster order
}
