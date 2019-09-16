#' Main function of the \pkg{\link{SaaRclust}} package to order and orient contigs of denovo genome assembly.
#'
#' This function is an easy-to-use wrapper to take Strand-seq reads aligned to the denovo genome assembly
#' and cluster each contig into chromosomal scaffold. In addition contigs are oriented and ordered within each
#' cluster/scaffold.
#'
#' @param outputfolder Folder to output the results. If it does not exist it will be created.
#' @param configfile A file specifying the parameters of this function (without \code{bamfolder}, \code{outputfolder} and \code{configfile}). 
#' If a \code{configfile} is specified, it will take priority over the command line parameters.
#' @param min.contig.size A minimal contig size to be processed (default: 100kb).
#' @param store.data.obj A logical indicating whether or not intermediate Rdata objects should be stored.
#' @param reuse.data.obj A logical indicating whether or not existing files in \code{outputfolder} should be reused.
#' @inheritParams importBams
#' @inheritParams hardClust
#' @inheritParams countProb
#' @inheritParams counts2ranges
#' @inheritParams connectDividedClusters
#' @inheritParams orderAndOrientClusters
#' @inheritParams exportPseudoChromosomalScaffolds
#' @return A \code{\link{GRanges-class}} object ???
#' @author David Porubsky
#' @export
#' 
scaffoldDenovoAssembly <- function(bamfolder, outputfolder, configfile=NULL, min.contig.size=100000, bin.size=100000, store.data.obj=TRUE, reuse.data.obj=FALSE, num.clusters=100, alpha=0.1, best.prob=1, prob.th=0, ord.method='TSP', assembly.fasta=NULL, concat.fasta=TRUE, z.limit=3, remove.always.WC=FALSE) {
  ## Get total processing time
  ptm <- proc.time()
  
  ## Set up the directory structure
  datapath <- file.path(outputfolder,'data')
  asmpath <- file.path(outputfolder,'clustered_assembly')
  pltpath <- file.path(outputfolder,'plots')
  
  ## Create directories if they do not exist
  if (!file.exists(outputfolder)) {
    dir.create(outputfolder)
  }
  if (!file.exists(datapath)) { dir.create(datapath) }
  if (!file.exists(asmpath)) { dir.create(asmpath) }
  if (!file.exists(pltpath)) { dir.create(pltpath) }
  
  ## Get user defined config file
  config <- NULL
  if (is.character(configfile)) {
    ## Read config file
    errstring <- tryCatch({
      config <- readConfig(configfile)
      errstring <- ''
    }, error = function(err) {
      errstring <- paste0("Could not read configuration file ",configfile)
    })
    if (errstring!='') {
      stop(errstring)
    }
  }
  
  ## Put all parameters into list and merge with config ##
  params <- list(min.contig.size=min.contig.size, bin.size=bin.size, store.data.obj=store.data.obj, 
                 reuse.data.obj=reuse.data.obj, num.clusters=num.clusters, alpha=alpha, 
                 best.prob=best.prob, prob.th=prob.th, ord.method=ord.method, assembly.fasta=assembly.fasta, 
                 concat.fasta=concat.fasta, z.limit=z.limit, remove.always.WC=remove.always.WC)
  config <- c(config, params[setdiff(names(params), names(config))])
  
  ## Make a copy of the config file ##
  writeConfig(config = config, configfile=file.path(outputfolder, 'SaaRclust.config'))
  
  ## Get contigs/scaffolds names and sizes
  bamFile <- list.files(bamfolder, pattern = ".bam$", full.names = TRUE)[1]
  file.header <- Rsamtools::scanBamHeader(bamFile)[[1]]
  chrom.lengths <- file.header$targets
  ## Keep only contigs/scaffolds >=100Kb
  chrom.lengths <- chrom.lengths[chrom.lengths >= config[['min.contig.size']]]
  chroms.in.data <- names(chrom.lengths)

  ## Get counts per genomic regions ##
  destination <- file.path(datapath, paste0("rawCounts_", config[['bin.size']],"bp_chunks.RData"))
  if (config[['reuse.data.obj']]) {
    if (file.exists(destination)) {
      message("Loading previously generated BAM read counts ...\n", destination)
      counts.l <- get(load(destination))
    } else {
      counts.l <- importBams(bamfolder = bamfolder, chromosomes = chroms.in.data, bin.size = config[['bin.size']])
    }
  } else {
    counts.l <- importBams(bamfolder = bamfolder, chromosomes = chroms.in.data, bin.size = config[['bin.size']])
  }
  ## Store data object
  if (config[['store.data.obj']]) {
    save(counts.l, file = destination)
  }

  ## Perform hard clustering ##
  destination <- file.path(datapath, paste0("hardClust_", config[['num.clusters']], "K_", config[['bin.size']], "bp_chunks.RData"))
  if (config[['reuse.data.obj']]) {
    if (file.exists(destination)) {
      message("Loading previously generated hard clustering results ...\n", destination)
      hardClust.ord <- get(load(destination))
    } else {
      set.seed(1000) ## to reproduce hard clustering results
      hardClust.ord <- hardClust(counts.l, num.clusters=config[['num.clusters']], nstart = 1000)
    }
  } else {
    set.seed(1000) ## to reproduce hard clustering results
    hardClust.ord <- hardClust(counts.l, num.clusters=config[['num.clusters']], nstart = 1000)
  }  
  ## Store data object
  if (config[['store.data.obj']]) {
    save(hardClust.ord, file = destination)
  }
  
  ## Estimate EM parameters
  theta.estim <- estimateTheta(counts.l, hard.clust=hardClust.ord, alpha=config[['alpha']])
  ## Set theta parameter
  theta.param <- theta.estim
  ## Set pi parameter
  readsPerCluts <- table(hardClust.ord)
  pi.param <- readsPerCluts / sum(readsPerCluts)

  ## RUN EM ##
  destination <- file.path(datapath, paste0("softClust_", config[['num.clusters']], "K_", config[['bin.size']], "bp_chunks.RData"))
  if (config[['reuse.data.obj']]) {
    if (file.exists(destination)) {
      message("Loading previously generated soft clustering results ...\n", destination)
      EM.obj <- get(load(destination))
    } else {
      EM.obj <- EMclust(counts.l, theta.param=theta.param, pi.param=pi.param, num.iter=20, alpha=config[['alpha']], logL.th=1, log.scale=TRUE)
    }
  } else {
    EM.obj <- EMclust(counts.l, theta.param=theta.param, pi.param=pi.param, num.iter=20, alpha=config[['alpha']], logL.th=1, log.scale=TRUE)
  }
  ## Store data object
  if (config[['store.data.obj']]) {
    save(EM.obj, file = destination)
  }
  
  ## Get cluster IDs that belong to the same chromosome/scaffold ##
  split.pairs <- connectDividedClusters(theta.param=EM.obj$theta.param, z.limit=config[['z.limit']], remove.always.WC=config[['remove.always.WC']])
  ## Store data object
  destination <- file.path(datapath, paste0("connectedClusters_", config[['bin.size']], "bp_chunks.RData"))
  if (config[['store.data.obj']]) {
    save(split.pairs, file = destination)
  }
  
  ## Assign contigs to clusters based on soft probabilities ##
  clustered.grl <- counts2ranges(counts.l, saarclust.obj=EM.obj, best.prob=config[['best.prob']], prob.th=config[['prob.th']])
  
  ## Order and orient contigs ##
  destination <- file.path(asmpath, paste0("ordered&oriented_", config[['bin.size']], "bp_chunks.tsv"))
  ordered.contigs.gr <- orderAndOrientClusters(clustered.grl=clustered.grl, split.pairs=split.pairs, ord.method=config[['ord.method']], alpha=config[['alpha']], bin.size=config[['bin.size']], filename=destination)
  ## Store data object
  destination <- file.path(asmpath, paste0("ordered&oriented_", config[['bin.size']], "bp_chunks.RData"))
  if (store.data.obj) {
    save(ordered.contigs.gr, file = destination)
  }
  
  ## Export clustered FASTA ##
  if (!is.null(config[['assembly.fasta']]) & is.character(config[['assembly.fasta']])) {
    if (file.exists(config[['assembly.fasta']])) {
      exportPseudoChromosomalScaffolds(clustered.gr=ordered.contigs.gr, assembly.fasta=config[['assembly.fasta']], outputfolder=asmpath, concat.fasta=config[['concat.fasta']])
    } else {
      message("Submitted 'assembly.fasta' does not exists !!!")
    }  
  }
  ## TODO return final clustered object??? 
  
  ## TODO make some plots ???
  
  ## Report total processing time
  time <- proc.time() - ptm
  message("\nTotal analysis time: ", round(time[3],2), "s")
}