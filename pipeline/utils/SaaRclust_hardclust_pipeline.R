#!/usr/bin/Rscript

#This Rscript runs kmeans based hard clustering implenented in package SaaRclust.
#author: David Porubsky

args=commandArgs(TRUE)
print(args)

#add user defined path to load needed libraries
.libPaths( c( .libPaths(), args[6]) )

suppressPackageStartupMessages(library(SaaRclust))

output <- runSaaRclust(inputfolder = args[1], outputfolder = args[2], num.clusters=as.numeric(args[3]), alpha = as.numeric(args[4]),  HC.only = TRUE, store.counts = FALSE, store.bestAlign = TRUE, numAlignments = as.numeric(args[5]))

