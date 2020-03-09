<img src="https://github.com/daewoooo/SaaRclust/raw/master/saarclust_logo.png" />

# SaaRclust
R Package to cluster and orient long sequencing reads or contigs into chromosomal scaffolds to facilitate de novo genome assembly. To do so this package utilize single-cell Strand-seq data (Falconer et al., 2012) that are able to preserve structural contiguity of individual homologs.

Collaborators: David Porubsky, Maryam Ghareghani, Peter Ebert and Tobias Marschall

## Installation

### Development version from Github
To install the development version from Github, follow the steps given below. The installation has only been tested on Ubuntu so far, if you need to install on Windows or Mac additional steps might be necessary (e.g. installation of Rtools from https://cran.r-project.org/bin/windows/Rtools/)

1. Install a recent version of R (>=3.5) from https://www.r-project.org/
2. Optional: For ease of use, install Rstudio from https://www.rstudio.com/
3. Open R and install all dependencies. Please ensure that you have writing permissions to install packages. Execute the following lines one by one:

  install.packages("devtools")
	source("http://bioconductor.org/biocLite.R")
	library(devtools)
	install_github("daewoooo/SaaRclust", branch="devel")
	Or alternatively if the above line doesn't work:
	install_git("git://github.com/daewoooo/SaaRclust.git", branch = "devel")
	
#### External Code
This package implements source code of orderContigsGreedy from [ContiBAIT](https://bioconductor.org/packages/contiBAIT) package. 
Copyright (c) 2015, Kieran O'Neill, Mark Hills, Mike Gottlieb

## Report Errors
If you encounter errors of any kind, please report an [issue here](https://github.com/daewoooo/SaaRclust/issues/new).

## NOTE
The 'devel' branch of SaaRclust package is currently under development and contains unpublished work. Any usage for publishing is strictly prohibited without permission.

## References
Ghareghani, M., Porubsky, D., et al. Strand-seq enables reliable separation of long reads by chromosome via expectation maximization. Bioinformatics 34, i115–i123 (2018).

Falconer, Ester, Mark Hills, Ulrike Naumann, Steven S. S. Poon, Elizabeth A. Chavez, Ashley D. Sanders, Yongjun Zhao, Martin Hirst, and Peter M. Lansdorp. 2012. “DNA Template Strand Sequencing of Single-Cells Maps Genomic Rearrangements at High Resolution.” Nature Methods 9 (11): 1107–12.

