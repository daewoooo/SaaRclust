#!/usr/bin/Rscript

args=commandArgs(TRUE)

source(args[1])
inputfolder <- args[2]
outputfolder <- args[3]

#load ggplot2 and biovizBase library
suppressPackageStartupMessages( library('ggplot2') )
suppressPackageStartupMessages( library('biovizBase') )

thresholds <- c(0, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.99)

message("Venn diagram statistics...")
venn.stats <- VennDiagramStats(inputfolder = inputfolder, thresholds = thresholds)
venn <- file.path(outputfolder, "VennStatistics.txt")
write.table(venn.stats, file = venn, sep = "\t", quote = F)

message("Preparing clustering accuracy lollipop plot ...")
ClustersAccuracyPerChrPerDir(inputfolder=inputfolder, thresholds=thresholds, minLib=5) -> accplt.obj

message("\nPreparing clustering accuracy boxplot ...")
boxplotDistSSlibsPerPB(inputfolder=inputfolder, thresholds=0) -> boxplt.obj

message("\nPreparing accuracy ranking plot ...")
accuracyRanking(inputfolder = inputfolder) -> rankingPlt.obj

message("\nExporting plot data ...")
destination <- file.path(outputfolder, "accPlot.RData") 
save(file = destination, accplt.obj)
destination <- file.path(outputfolder, "accBoxplot.RData") 
save(file = destination, boxplt.obj)
destination <- file.path(outputfolder, "rankingPlot.RData") 
save(file = destination, rankingPlt.obj)
