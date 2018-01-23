inputfolder <- "/media/daewoooo/WORK/Clustering_project/WholeGenomeAnalysis/SaaRclust_results/RawData/"

processClusters <- function(inputfolder=NULL) {
  files2process <- list.files(inputfolder, pattern = "dataQuals.RData", full.names = TRUE)
   
  SSreads.perPB <- list()
  SSlib.perPB <- list()
  PBreadLenDist <- list()
  for (file in files2process) {
    data.file <- get(load(file))
    fileID <- basename(file)
    
    SSreads.perPB[[fileID]] <- data.file$SSreads.perPB
    SSlib.perPB[[fileID]] <- data.file$SSlib.perPB
    PBreadLenDist[[fileID]] <- data.file$PBreadLenDist
  }
  SSreads.perPB.all <- do.call(c, SSreads.perPB)
  SSlib.perPB.all <- do.call(rbind, SSlib.perPB)
  
  #plot data quality
  SSreads.perPB.all.df <- as.data.frame(SSreads.perPB.all)
  quantil0.09 <- quantile(SSreads.perPB.all.df$SSreads.perPB.all, probs = 0.9)
  SSreads.perPB.all.df <- data.frame(SSreads.perPB=SSreads.perPB.all.df[SSreads.perPB.all.df$SSreads.perPB.all <= quantil0.09,])
  SSreads.perPB.plt <- ggplot(SSreads.perPB.all.df, aes(x=SSreads.perPB)) + geom_histogram(binwidth = 1, fill="red")
  
  SSlib.perPB.all.df <- data.frame(SSlib.perPB=SSlib.perPB.all$counts)
  SSlib.perPB.plt <- ggplot(SSlib.perPB.all.df, aes(x=SSlib.perPB)) + geom_histogram(binwidth = 1, fill="red") + geom_vline(xintercept = 10)
  
}  