## Load required libraries
library("gsubfn")
library("data.table")
library("RColorBrewer")
library("SaaRclust")


getAccuracy <- function(file){
  message("Working on ",file)
  path <- dirname(file)
  file <- basename(file)

  options <- unlist(strsplit(basename(file), split = "_"))[c(4:7)]

  filePath <- file.path(path, file)
  conn <- file(filePath, open="r")
  lines <-readLines(conn)
  lines.l <- list()
  for (i in 1:length(lines)){
    lines.l[[i]] <- lines[i]
  }
  close(conn)

  minimizers <- strapplyc(lines.l[[4]], "consider: (\\d+)", simplify = T)
  real.time <- strapplyc(lines.l[[7]], "Real time: (\\d+.\\d+)", simplify = T)
  cpu.time <- strapplyc(lines.l[[7]], "CPU: (\\d+.\\d+)", simplify = T)

  k <- strapplyc(options[1], "(\\d+)", simplify = T)
  w <- strapplyc(options[2], "(\\d+)", simplify = T)
  L <- strapplyc(options[3], "(\\d+)", simplify = T)
  f <- strapplyc(options[4], "(\\d+.\\d+)", simplify = T)

  minimap.file <- gsub(file, pattern = '\\.log', replacement = '')
  minimap.file.path <- file.path(path, minimap.file)
  tab.in <- importData(infile = minimap.file.path)

  ## Export data quality measures
  qual.measures <- getQualMeasure(tab.in)

  mean.SSreads.perPB <- round(mean(qual.measures$SSreads.perPB))
  mean.SSlib.perPB <- round(mean(qual.measures$SSlib.perPB$counts))
  mean.SSreads.perlib.perPB <- round(mean(qual.measures$SSreads.perlib.perPB))
  accur.counts <- table(tab.in$SSchrom == tab.in$PBchrom)
  falses <- round((accur.counts[1]/sum(accur.counts))*100, digits = 2)
  trues <- round((accur.counts[2]/sum(accur.counts))*100, digits = 2)

  total.PBreads <- length(unique(tab.in$PBreadNames))
  total.SSreads <- length(unique(tab.in$SSreadNames))

  v <- data.frame(k=as.numeric(k), w=as.numeric(w), L=as.numeric(L), f=as.numeric(f), num.minimizers=as.numeric(minimizers), real.time=as.numeric(real.time), cpu.time=as.numeric(cpu.time), mean.SSreads.perPB=mean.SSreads.perPB, mean.SSlib.perPB=mean.SSlib.perPB, mean.SSreads.perlib.perPB=mean.SSreads.perlib.perPB, total.PBreads=as.numeric(total.PBreads), total.SSreads=as.numeric(total.SSreads), falses=falses, trues=trues)
  
 return(v)
}
