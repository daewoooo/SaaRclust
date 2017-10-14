
#' Wrapper function to run saarclust pipeline for a given number of PB reads.
#'
#' @param minimap.file A path to the minimap file to load.
#' @param num.clusters Expected number of clusters. (for 22 autosomes == 44 clusters)
#' @param EM.iter Number of iteration to run EM for.
#' @param chunk.size Number of PB reads to precess.
#' @author David Porubsky, Maryam Ghareghani

#NOTE: I suggest to run this code line by line for debugging purposes

#minimap test files are present in TestData folder
#minimap.file <- "/media/daewoooo/WORK/SS2PacBio_alignment_HG00733/WholeGenomeAnalysis/subsetMinimap.txt"
#minimap.file <- "/media/daewoooo/WORK/SS2PacBio_alignment_HG00733/Test_cluster_chr21&chr22/Minimap_out/SS2Pacbio_minimap_HG00733_k13_w1_L70_f0.01_Chr21andChr22_allSSreads"

#load the function below into R if you want to run all steps in one command

runSaaRclust <- function(minimap.file=NULL, num.clusters=44, EM.iter=50, chunk.size=10000) {

  tab.in <- importTestData(infile = minimap.file)
  tab.in <- tab.in[tab.in$SSchrom != 'chrUn' & tab.in$SSchrom != 'chrX',] #applies only for test data
  tab.filt <- filterInput(inputData=tab.in, quantileSSreads=c(0.4,0.9))

  #take a smaller chunk of PB reads to process
  tab.filt <- tab.filt[sample(nrow(tab.filt)),] #shuffle rows in tab
  chunk <- unique(tab.filt$PBreadNames)[1:chunk.size]
  tab.filt <- tab.filt[tab.filt$PBreadNames %in% chunk,]

  #additional sort by direction
  tab.filt <- tab.filt[order(tab.filt$PBflag),]

  #use PB read names as factor
  tab.filt <- tab.filt[order(tab.filt$PBchrom),]
  tab.filt$PBreadNames <- factor(tab.filt$PBreadNames, levels=unique(tab.filt$PBreadNames))

  #get PB chrom names from the ordered PB reads
  chr.l <- split(tab.filt$PBchrom, tab.filt$PBreadNames)
  chr.rows <- sapply(chr.l, unique)

  #get PB directionality from the ordered PB reads
  pb.flag <- split(tab.filt$PBflag, tab.filt$PBreadNames)
  pb.flag <- sapply(pb.flag, unique)

  #get PB read position
  pb.pos <- split(tab.filt$PBpos, tab.filt$PBreadNames)
  pb.pos <- sapply(pb.pos, unique)

  #split data by SS library
  tab.l <- split(tab.filt, tab.filt$SSlibNames)

  ### Hard clustering ###
  hard.clust <- hardClust(tab.l, clusters=num.clusters)
  
  #get accuracy of hard clustering [OPTIONAL]
  #chr.clusts <- split(chr.rows, hard.clust$clust.id)
  #clust.acc <- getClusterAcc(chr.clusts)
  
  #initialize thetas
  theta.l <- hard.clust$theta.estim

  #estimate pi based on # of PB reads in each cluster
  readsPerCluts <- table(hard.clust$clust.id)
  pi.param <- readsPerCluts/sum(readsPerCluts)

  ### EM algorithm ###
  soft.clust <- saarclust(tab.l, theta.l=theta.l, pi.param=pi.param, num.iter=EM.iter, raw.counts=hard.clust$raw.counts)

  #plot resultant clusters
  soft.clust.df <- as.data.frame(soft.clust$soft.pVal)
  soft.clust.df$PBchrom <- chr.rows
  soft.clust.df$PBflag <- pb.flag

  ha1 = rowAnnotation(df = data.frame(chr = unname(chr.rows)))
  ha2 = rowAnnotation(df = data.frame(PB.dir = unname(pb.flag)), col = list(PB.dir = c('0'="chocolate1", '16'="chartreuse3")))
  hm <- Heatmap(soft.clust.df[,c(1:num.clusters)], name = "Probs", cluster_columns = F, cluster_rows = F, show_row_names = FALSE)
  hm + ha1 + ha2
}


  