inputfolder <- "/media/daewoooo/WORK/Clustering_project/WholeGenomeAnalysis/SaaRclust_results/RawData/"
processClusters(inputfolder) -> dataQual.plt

processClusters <- function(inputfolder=NULL) {
  files2process <- list.files(inputfolder, pattern = "dataQuals.RData", full.names = TRUE)
   
  SSreads.perPB <- list()
  SSlib.perPB <- list()
  SSreads.perlib.perPB <- list()
  PBreadLenDist <- list()
  for (file in files2process) {
    data.file <- get(load(file))
    fileID <- basename(file)
    
    SSreads.perPB[[fileID]] <- data.file$SSreads.perPB
    SSlib.perPB[[fileID]] <- data.file$SSlib.perPB
    PBreadLenDist[[fileID]] <- data.file$PBreadLenDist
    #get counts and export only max 100
    SSreads.perlib.perPB[[fileID]] <- sort(table(data.file$SSreads.perlib.perPB),decreasing = T)[1:100]
  }
  SSreads.perPB.all <- do.call(c, SSreads.perPB)
  SSlib.perPB.all <- do.call(rbind, SSlib.perPB)
  SSreads.perlib.perPB.all <- Reduce("+", SSreads.perlib.perPB)
  
  max <- which.max(sapply(PBreadLenDist, nrow))
  size.ids <- PBreadLenDist[[max]]$midpoints
  
  addMissingLengths <- function(x) {
    ids2add <- size.ids[!size.ids %in% x[[1]]]
    if (!length(ids2add) == 0) {  
      add.data <- data.frame(midpoints=ids2add, freq=rep(0, length(ids2add)))
      merged.data <- rbind(x, add.data)
      return(merged.data[,2])
    } else {
      return(x[,2])
    }  
  }
  
  PBreadLenDist.modif <- lapply(PBreadLenDist, addMissingLengths)
  PBreadLenDist.modif <- PBreadLenDist.modif[lengths(PBreadLenDist.modif) == length(size.ids)]
  PBreadLenDist.all <- Reduce("+",  PBreadLenDist.modif)
  PPBreadLenDist.df <- data.frame(id=size.ids, counts=PBreadLenDist.all)
  
  #plot data quality
  SSreads.perPB.all.df <- as.data.frame(SSreads.perPB.all)
  quantil0.09 <- quantile(SSreads.perPB.all.df$SSreads.perPB.all, probs = 0.9)
  SSreads.perPB.all.df <- data.frame(SSreads.perPB=SSreads.perPB.all.df[SSreads.perPB.all.df$SSreads.perPB.all <= quantil0.09,])
  SSreads.perPB.plt <- ggplot(SSreads.perPB.all.df, aes(x=SSreads.perPB)) + geom_histogram(binwidth = 1, fill="red") + xlab("# of StrandS reads per PBread") + ylab("Frequency")
  
  SSlib.perPB.all.df <- data.frame(SSlib.perPB=SSlib.perPB.all$counts)
  SSlib.perPB.plt <- ggplot(SSlib.perPB.all.df, aes(x=SSlib.perPB)) + geom_histogram(binwidth = 1, fill="red") + geom_vline(xintercept = 10) + xlab("# of StrandS libraries per PBread") + ylab("Frequency")
  
  SSreads.perlib.perPB.df <- as.data.frame(SSreads.perlib.perPB.all)
  SSreads.perlib.perPB.df <- SSreads.perlib.perPB.df[SSreads.perlib.perPB.df$Var1 %in% c(1:20),]
  SSreads.perlib.perPB.plt <- ggplot(SSreads.perlib.perPB.df, aes(x=Var1, y=Freq)) + geom_bar(fill="red", stat="identity") + xlab("# of ShortReads per PBread per Library") + ylab("Frequency")

  PPBreadLenDist.plt <- ggplot(PPBreadLenDist.df, aes(x=id, y=counts)) + geom_bar(fill="red", stat="identity") + xlab("PacBio read length (bp)") + ylab("Frequency") + scale_x_continuous(breaks = c(10000, 20000, 40000, 60000, 80000))
  
  main.plt <- plot_grid(SSreads.perPB.plt, SSlib.perPB.plt, SSreads.perlib.perPB.plt, PPBreadLenDist.plt, nrow = 1)
  return(main.plt)
}  