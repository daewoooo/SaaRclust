#!/usr/bin/env Rscript

#This Rscript installs SaaRclust package from the github repository along with all dependencies.
#author: David Porubsky

Sys.setenv(Renv='PWD')
library(devtools)

source("https://bioconductor.org/biocLite.R")
biocLite("ComplexHeatmap", lib="utils/R-packages")

withr::with_libpaths(new = "utils/R-packages", install_git("git://github.com/daewoooo/SaaRclust.git", branch = "master"), "prefix")
