#!/usr/bin/Rscript

#This Rscript runs EM (Expectation Maximization) based soft clustering algorithm implenented in package SaaRclust.
#author: David Porubsky

args=commandArgs(TRUE)
print(args)

#add user defined path to load needed libraries
.libPaths( c( .libPaths(), args[10]) )

suppressPackageStartupMessages(library(SaaRclust))

output <- SaaRclust(minimap.file=args[1], outputfolder=args[2], num.clusters=args[3], EM.iter=args[4], alpha=as.numeric(args[5]), minLib=as.numeric(args[6]), upperQ=as.numeric(args[7]), logL.th=args[8], theta.constrain=FALSE, HC.input=args[9])


