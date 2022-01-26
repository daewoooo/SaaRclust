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
#' @param mask.regions Set to \code{TRUE} if regions that appear as WC in majority of cells and low coverage regions should be masked.
#' @param eval.ploidy If set to \code{TRUE} estimated ploidy of each contig will be reported and appended to each contig name.
#' @param allow.contig.cuts If set to \code{TRUE} contigs that were assigned to more than one cluster or direction will be cut.
#' @param em.param.estim A \code{\link{SaaRclust}} object with theta and pi estimates to be used in EM procedure (Soft clustering step).
#' @inheritParams importBams
#' @inheritParams hardClust
#' @inheritParams countProb
#' @inheritParams counts2ranges
#' @inheritParams connectDividedClusters
#' @inheritParams orderAndOrientClusters
#' @inheritParams exportPseudoChromosomalScaffolds
#' @importFrom Rsamtools indexFa FaFile scanBamHeader scanFaIndex
#' @importFrom BiocGenerics table as.list
#' @return NULL
#' @author David Porubsky
#' @export
#' 
#' @examples
#'\dontrun{
#'## Required parameters to run SaaRclust on BAM files stored "bam-data-folder" using default settings.
#'## To export clustred FASTA file, an original FASTA used in BAM alignments has to be submitted as 'assembly.fasta'.
#'scaffoldDenovoAssembly(bamfolder="bam-data-folder", outputfolder="saarclust-output-folder")}
#'
scaffoldDenovoAssembly <- function(bamfolder, outputfolder, configfile=NULL, min.mapq=10, min.contig.size=100000, min.region.to.order=0, chromosomes=NULL, pairedEndReads=TRUE, custom.bins = NULL, bin.size=100000, step.size=NULL, bin.method='fixed', store.data.obj=TRUE, reuse.data.obj=FALSE, num.clusters=100, desired.num.clusters=NULL, max.cluster.length.mbp=0, alpha=0.1, prob.th=0, ord.method='TSP', assembly.fasta=NULL, concat.fasta=TRUE, z.limit=3.29, remove.always.WC=FALSE, mask.regions=FALSE, eval.ploidy=FALSE, allow.contig.cuts=FALSE, em.param.estim=NULL) {
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
  params <- list(min.mapq=min.mapq, min.contig.size=min.contig.size, min.region.to.order=min.region.to.order, chromosomes=chromosomes, pairedEndReads=pairedEndReads, custom.bins=custom.bins, 
                 bin.size=bin.size, store.data.obj=store.data.obj, step.size=step.size, bin.method=bin.method, reuse.data.obj=reuse.data.obj, num.clusters=num.clusters,
                 desired.num.clusters=desired.num.clusters, max.cluster.length.mbp=max.cluster.length.mbp, alpha=alpha, prob.th=prob.th, ord.method=ord.method, 
                 assembly.fasta=assembly.fasta, concat.fasta=concat.fasta, z.limit=z.limit, remove.always.WC=remove.always.WC, mask.regions=mask.regions, eval.ploidy=eval.ploidy, 
                 allow.contig.cuts=allow.contig.cuts, em.param.estim=em.param.estim)
  config <- c(config, params[setdiff(names(params), names(config))])
  
  ## Make a copy of the config file ##
  writeConfig(config = config, configfile=file.path(outputfolder, 'SaaRclust.config'))
  
  ## Get contigs/scaffolds names and sizes from BAM
  bamfile <- list.files(bamfolder, pattern = ".bam$", full.names = TRUE)[1]
  file.header <- Rsamtools::scanBamHeader(bamfile)[[1]]
  chrom.lengths <- file.header$targets
  
  ## Get contigs/scaffolds names and sizes from FASTA
  if (!is.null(config[['assembly.fasta']]) & is.character(config[['assembly.fasta']])) {
    ## Check if submitted fasta file is indexed
    assembly.fasta.idx <- paste0(config[['assembly.fasta']], ".fai")
    if (!file.exists(assembly.fasta.idx)) {
      ptm <- startTimedMessage("Fasta file is not indexed, indexing ")
      fa.idx <- Rsamtools::indexFa(file = config[['assembly.fasta']])
      stopTimedMessage(ptm)
    }
    fa.file <- open(Rsamtools::FaFile(config[['assembly.fasta']]))
    fa.idx <- Rsamtools::scanFaIndex(fa.file)
    ## Check if BAMs have been aligned to thr FASTA file submitted as 'assembly.fasta'
    if (!all(names(chrom.lengths) %in% GenomeInfoDb::seqlevels(fa.idx))) {
      warning("Not all sequence names in BAMs match those in submitted FASTA file.
              Final FASTA file cannot be exported, setting 'assembly.fasta=NULL'.")
      assembly.fasta <- NULL
    } else {
      ## Check if submitted BAM and FASTA files have the same seqlengths
      if(!all(chrom.lengths == GenomeInfoDb::seqlengths(fa.idx)[names(chrom.lengths)])) {
        warning("Not all sequence lengths in BAMs match those in submitted FASTA file.
                 Final FASTA file cannot be exported, setting 'assembly.fasta=NULL'.")
        assembly.fasta <- NULL
      }
    }
  }
  
  ## Keep only contigs/scaffolds of size defined by min.contig.size
  filt <- chrom.lengths >= config[['min.contig.size']]
  chroms.in.data <- names(chrom.lengths[filt])
  
  ## Process only user defined chromosomes/contigs
  if (!is.null(chromosomes) & is.character(chromosomes)) {
    if (any(chromosomes %in% chroms.in.data)) {
      chroms.in.data <- chroms.in.data[chroms.in.data %in% chromosomes]
      chrom.lengths <- chrom.lengths[names(chrom.lengths) %in% chromosomes]
    } else {
      warning("None of the sequence names defined in 'chromosomes' found, after 'min.contig.size' filtering !!!")
    }
  }
  
  ## Create a mask of regions with and excess of read coverage that appears always WC
  ## and bins that have very low read counts
  destination <- file.path(datapath, paste0("maskRegions.RData"))
  if (config[['mask.regions']]) {
    if (config[['reuse.data.obj']]) {
      if (file.exists(destination)) {
        message("Loading previously generated region mask ...\n", destination)
        blacklist <- get(load(destination))
        #blacklist.gr <- c(blacklist$alwaysWC, blacklist$alwaysZero)
      } else {
        bins.gr <- makeFixedBins(bamfile = bamfile, bin.size = 50000, step.size = 10000, chromosomes = chroms.in.data)
        blacklist <- suppressWarnings( 
          maskAlwaysWCandZeroBins(bamfolder = bamfolder, genomic.bins = bins.gr, min.mapq = config[['min.mapq']], pairedEndReads = config[['pairedEndReads']])
        )
      }  
    } else {
      bins.gr <- makeFixedBins(bamfile = bamfile, bin.size = 50000, step.size = 10000, chromosomes = chroms.in.data)
      blacklist <- suppressWarnings( 
        maskAlwaysWCandZeroBins(bamfolder = bamfolder, genomic.bins = bins.gr, min.mapq = config[['min.mapq']], pairedEndReads = config[['pairedEndReads']])
      )
    } 
    blacklist.gr <- c(blacklist$alwaysWC, blacklist$alwaysZero)
    ## Store data object
    if (config[['store.data.obj']]) {
      save(blacklist, file = destination)
    }
  }

  ## Get counts per genomic regions ##
  destination <- file.path(datapath, paste0("rawCounts_", config[['bin.size']],"bp_", config[['bin.method']], ".RData"))
  if (config[['reuse.data.obj']]) {
    if (file.exists(destination)) {
      message("Loading previously generated BAM read counts ...\n", destination)
      counts.l <- get(load(destination))
    } else {
      if (config[['mask.regions']]) {
        counts.l <- importBams(bamfolder = bamfolder, chromosomes = chroms.in.data, pairedEndReads = config[['pairedEndReads']], min.mapq = config[['min.mapq']], custom.bins = config[['custom.bins']], bin.size = config[['bin.size']], step.size = config[['step.size']], bin.method = config[['bin.method']], blacklist = blacklist.gr)
      } else {
        counts.l <- importBams(bamfolder = bamfolder, chromosomes = chroms.in.data, pairedEndReads = config[['pairedEndReads']], min.mapq = config[['min.mapq']], custom.bins = config[['custom.bins']], bin.size = config[['bin.size']], step.size = config[['step.size']], bin.method = config[['bin.method']])
      }
      ## Remove bins with zero counts across all cells
      counts.sums <- Reduce('+', counts.l)
      counts.sums <- rowSums(counts.sums)
      zero.bins <- which(counts.sums == 0)
      if (length(zero.bins) > 0) {
        counts.l <- lapply(counts.l, function(x) x[-zero.bins,])
        message("Removed zero count bins ", length(zero.bins), "/", length(counts.sums))
      }  
    }
  } else {
    if (config[['mask.regions']]) {
      counts.l <- importBams(bamfolder = bamfolder, chromosomes = chroms.in.data, pairedEndReads = config[['pairedEndReads']], min.mapq = config[['min.mapq']], custom.bins = config[['custom.bins']], bin.size = config[['bin.size']], step.size = config[['step.size']], bin.method = config[['bin.method']], blacklist = blacklist.gr)
    } else {
      counts.l <- importBams(bamfolder = bamfolder, chromosomes = chroms.in.data, pairedEndReads = config[['pairedEndReads']], min.mapq = config[['min.mapq']], custom.bins = config[['custom.bins']], bin.size = config[['bin.size']], step.size = config[['step.size']], bin.method = config[['bin.method']])
    }
    ## Remove bins with zero counts across all cells
    counts.sums <- Reduce('+', counts.l)
    counts.sums <- rowSums(counts.sums)
    zero.bins <- which(counts.sums == 0)
    if (length(zero.bins) > 0) {
      counts.l <- lapply(counts.l, function(x) x[-zero.bins,])
      message("Removed zero count bins ", length(zero.bins), "/", length(counts.sums))
    }
  }
  ## Store data object
  if (config[['store.data.obj']]) {
    save(counts.l, file = destination)
  }
  
  ## Check user input of 'em.param.estim' ##
  if (!is.null(config[['em.param.estim']])) {
    ## Load theta and pi parameter estimates from the user defined saarclust object in 'em.param.estim'
    saarclust.obj <- get(load(config[['em.param.estim']]))
    if (class(saarclust.obj) == 'saarclust') {
      theta.param <- saarclust.obj$theta.param
      pi.param <- saarclust.obj$pi.param
      if (length(theta.param) == length(counts.l)) {
        message("Loading previously generated 'theta' and 'pi' parameter estimates ...\n", config[['em.param.estim']])
      } else {
        warning("Submitted 'em.param.estim' contains ", length(theta.param), " cells while binned counts were created from ", length(counts.l), " cells, setting 'em.param.estim' to NULL." )
        config[['em.param.estim']] <- NULL
      }
    } else {
      warning("Submitted 'em.param.estim' is not a 'saarclust' object. Setting 'em.param.estim' to NULL.")
      config[['em.param.estim']] <- NULL
    }
  }
  
  ## Perform hard clustering ##
  destination <- file.path(datapath, paste0("hardClust_", config[['num.clusters']], "K_", config[['bin.size']], "bp_", config[['bin.method']], ".RData"))
  if (config[['reuse.data.obj']] & is.null(config[['em.param.estim']])) {
    if (file.exists(destination)) {
      message("Loading previously generated hard clustering results ...\n", destination)
      hardClust.ord <- get(load(destination))
    } else {
      set.seed(1000) ## to reproduce hard clustering results
      hardClust.ord <- hardClust(counts.l, num.clusters=config[['num.clusters']], nstart = 1000)
    }
  } else if (is.null(config[['em.param.estim']])) {
    set.seed(1000) ## to reproduce hard clustering results
    hardClust.ord <- hardClust(counts.l, num.clusters=config[['num.clusters']], nstart = 1000)
  }
  ## Store data object
  if (config[['store.data.obj']]) {
    if (exists('hardClust.ord')) {
      save(hardClust.ord, file = destination)
    }  
  }
  
  ## Estimate EM parameters
  if (is.null(config[['em.param.estim']])) {
    ## Estimate theta parameter from hard clustering results
    theta.param <- estimateTheta(counts.l, hard.clust=hardClust.ord, alpha=config[['alpha']])
    ## Estimate pi parameter from hard clustering results
    readsPerCluts <- BiocGenerics::table(hardClust.ord)
    pi.param <- readsPerCluts / sum(readsPerCluts)
  }

  ## RUN EM ##
  destination <- file.path(datapath, paste0("softClust_", config[['num.clusters']], "K_", config[['bin.size']], "bp_", config[['bin.method']], ".RData"))
  if (config[['reuse.data.obj']]) {
    if (file.exists(destination) & is.null(config[['em.param.estim']])) {
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
  
  ## Find clusters with WC state in majority of cells ##
  theta.sums <- Reduce("+", EM.obj$theta.param)
  theta.zscore <- (theta.sums[,3] - mean(theta.sums[,3])) / sd(theta.sums[,3])
  wc.clust.idx <- which(theta.zscore >= 2.576) ## 99% confdence level
  # if (length(wc.clust.idx) == 0) {
  #   wc.clust.idx <- which.max(theta.zscore)
  # }
  
  ## Get names of always WC regions
  if (length(wc.clust.idx) > 0) {
    ## Get probability table
    soft.prob <- EM.obj$soft.pVal
    ## Remove segments that do not reach required prob.th
    row.max.prob <- apply(soft.prob, 1, max)
    mask <- row.max.prob >= config[['prob.th']]
    soft.prob <- soft.prob[mask,]
    clust.ID <- apply(soft.prob, 1, which.max)
    always.wc.ctgs <- names(clust.ID)[clust.ID %in% wc.clust.idx]
    always.wc.ctgs.gr <- string2GRanges(always.wc.ctgs)
    always.wc.ctgs.gr <- always.wc.ctgs.gr[,0]
    ## Store data object
    destination <- file.path(datapath, paste0("alwaysWCregions", config[['bin.size']], "bp_", config[['bin.method']], ".RData"))
    if (store.data.obj) {
      save(always.wc.ctgs.gr, file = destination)
    }
  } else {
    always.wc.ctgs.gr <- NULL
  }
  
  ## Remove always WC clusters
  if (config[['remove.always.WC']]) {
    if (length(wc.clust.idx) > 0) {
      ## Remove cluster with the most WC states
      message("Removing cluster ", paste(wc.clust.idx, collapse = ", "), " with the most WC states !!!")
      EM.obj$theta.param <- lapply(EM.obj$theta.param, function(x) x[-wc.clust.idx,])
      EM.obj$soft.pVal <- EM.obj$soft.pVal[,-wc.clust.idx]
      EM.obj$pi.param <- EM.obj$pi.param[-wc.clust.idx]
    }
  }
  
  ## Assign contigs to clusters based on soft probabilities ##
  clustered.grl <- counts2ranges(counts.l, saarclust.obj=EM.obj, best.prob=1, prob.th=config[['prob.th']])
  ## Get regions assigned to each cluster
  regions.clustered <- clustered.grl[[1]][,1]
  
  ## Find clusters with WW and CC state in majority of cells [haploid clusters] ##
  theta.sums <- Reduce("+", EM.obj$theta.param)
  theta.zscore.hap <- (theta.sums[,3] - mean(theta.sums[,3])) / sd(theta.sums[,3])
  hap.clust.idx <- which(theta.zscore.hap <= -2)
  if (length(hap.clust.idx) > 0) {
    message("Haploid clusters detected ", paste(hap.clust.idx, collapse = ", "), " !!!")
    ## Get haploid contigs
    ctg2clusters <- clustered.grl[[1]][,1]
    hap.ctgs <- GenomicRanges::reduce(ctg2clusters[ctg2clusters$clust.ID %in% hap.clust.idx])
    ## Store data object
    destination <- file.path(datapath, paste0("haploid_contig_regions.RData"))
    save(hap.ctgs, file = destination)
  } else {
    message("NO haploid clusters detected !!!")
  }
  
  ## Get cluster IDs that belong to the same chromosome/scaffold ##
  nclust <- nrow(EM.obj$theta.param[[1]])
  if (!is.null(config[['desired.num.clusters']])) {
    if (nclust > 2 & config[['num.clusters']] > config[['desired.num.clusters']]) {
      split.pairs <- connectDividedClusters(theta.param=EM.obj$theta.param, 
                                            clustered.gr=regions.clustered, 
                                            z.limit=config[['z.limit']], 
                                            desired.num.clusters=config[['desired.num.clusters']], 
                                            max.cluster.length.mbp=config[['max.cluster.length.mbp']])
    } else {
      clusters <- BiocGenerics::as.list(c(1:nclust))
      names(clusters) <- c(1:nclust)
      split.pairs <- list(clusters=clusters, putative.HETs=NULL)
    }
  } else if (nclust > 2) {
    split.pairs <- connectDividedClusters(theta.param=EM.obj$theta.param,
                                          clustered.gr=regions.clustered,
                                          z.limit=config[['z.limit']], 
                                          desired.num.clusters=config[['desired.num.clusters']],
                                          max.cluster.length.mbp=config[['max.cluster.length.mbp']])
  } else {
    clusters <- as.list(c(1:nclust))
    names(clusters) <- c(1:nclust)
    split.pairs <- list(clusters=clusters, putative.HETs=NULL)
  }
  
  ## Store data object
  destination <- file.path(datapath, paste0("connectedClusters_", config[['bin.size']], "bp_", config[['bin.method']], ".RData"))
  if (config[['store.data.obj']]) {
    save(split.pairs, file = destination)
  }
  
  ## Assign contigs to clusters based on soft probabilities ##
  #clustered.grl <- counts2ranges(counts.l, saarclust.obj=EM.obj, best.prob=1, prob.th=config[['prob.th']])
  
  ## Order and orient contigs ##
  #destination <- file.path(asmpath, paste0("ordered&oriented_", config[['bin.size']], "bp_", config[['bin.method']], ".tsv"))
  destination <- NULL
  ordered.contigs.gr <- orderAndOrientClusters(clustered.grl=clustered.grl, split.pairs=split.pairs, ord.method=config[['ord.method']], alpha=config[['alpha']], min.region.to.order=config[['min.region.to.order']], filename=destination)
  GenomeInfoDb::seqlengths(ordered.contigs.gr) <- chrom.lengths[GenomeInfoDb::seqlevels(ordered.contigs.gr)]
  
  ## Extend gaps between ranges
  ordered.contigs.gr <- expandGaps(ordered.contigs.gr)
  
  ## Add contig ploidy information ##
  if (config[['eval.ploidy']]) {
    if (exists('hap.ctgs')) {
      if (length(hap.ctgs) > 0) {
        ordered.contigs.gr <- labelGenomicRegions(gr = ordered.contigs.gr, label.gr = hap.ctgs, label.gr.ID = '1n')
      }  
    }  
    if(!config[['remove.always.WC']] & length(always.wc.ctgs.gr) > 0) {
      ordered.contigs.gr <- labelGenomicRegions(gr = ordered.contigs.gr, label.gr = always.wc.ctgs.gr, label.gr.ID = '>2n')
    }  
  }
  
  ## Report contigs assigned to more than one cluster or with putative misorient ##
  putative.errors <- ordered.contigs.gr
  putative.errors.grl <- GenomicRanges::split(putative.errors, GenomicRanges::seqnames(putative.errors))
  putative.errors.grl <- putative.errors.grl[lengths(putative.errors.grl) > 1]
  mask <- lapply(putative.errors.grl, function(gr) length(unique(gr$ID)) > 1 | length(unique(gr$dir)) > 1)
  idx <- which(mask == TRUE)
  if (length(idx) > 0) {
    putative.errors.gr <- unlist(putative.errors.grl[idx], use.names = FALSE)
    putative.errors.gr  <- GenomeInfoDb::keepSeqlevels(putative.errors.gr, value = unique(GenomicRanges::seqnames(putative.errors.gr)))
    ## Collapse consecutive ranges
    putative.errors.gr <- GenomicRanges::sort(putative.errors.gr)
    putative.errors.gr$collapse.ID <- paste0(putative.errors.gr$dir, '_', putative.errors.gr$ID)
    id.field <- ncol(mcols(putative.errors.gr))
    putative.errors.gr <- SaaRclust::collapseBins(putative.errors.gr, id.field = id.field)
    ## Store data object
    destination <- file.path(datapath, paste0("putativeAsmErrors_", config[['bin.size']], "bp_", config[['bin.method']], ".RData"))
    if (store.data.obj) {
      save(putative.errors.gr, file = destination)
    }
    ## Export summary table of assembly errors
    putative.errors.report <- reportMisAsmCTGs(gr = putative.errors.gr)
    destination <- file.path(asmpath, paste0("asmErrorsReport_", config[['bin.size']], "bp_", config[['bin.method']], ".tsv"))
    utils::write.table(putative.errors.report, file = destination, quote = FALSE, row.names = FALSE, append = FALSE, sep = "\t")
  } else {
    putative.errors.report <- NULL
  }

  ## UNCut putative assembly errors ##
  if (!config[['allow.contig.cuts']]) {
    ## Concatenate divided contigs (assigned to more than one cluster or direction) into a single contig range
    grl <- putative.errors.grl
    concat.grl <- GenomicRanges::GRangesList()
    for (i in seq_along(grl)) {
      gr.sub <- grl[[i]]
      ## Get the largest chunk of a divided contig
      gr.sub.max <- gr.sub[which.max(width(gr.sub))]
      ## Get full contig size
      gr.full <- range(gr.sub)
      ## Add mcols from the largest chunk of a divided contig
      GenomicRanges::mcols(gr.full) <- GenomicRanges::mcols(gr.sub.max)
      concat.grl[[i]] <- gr.full
    }
    gr <- unlist(concat.grl, use.names = FALSE)
    ## Replace divided contigs with their continuous range
    ordered.contigs.gr <- ordered.contigs.gr[!as.character(GenomicRanges::seqnames(ordered.contigs.gr)) %in% as.character(GenomicRanges::seqnames(gr))]
    ordered.contigs.gr <- GenomicRanges::sort(c(ordered.contigs.gr, gr))
  }

  ## Store data object
  destination <- file.path(datapath, paste0("ordered&oriented_", config[['bin.size']], "bp_chunks.RData"))
  if (store.data.obj) {
    save(ordered.contigs.gr, file = destination)
  }
  
  ## Report contigs that have been filtered out
  filtered.ctgs <- chrom.lengths[!names(chrom.lengths) %in% GenomeInfoDb::seqlevels(ordered.contigs.gr)]
  if (length(filtered.ctgs) > 0) {
    filtered.ctgs.gr <- GenomicRanges::GRanges(seqnames = names(filtered.ctgs), ranges=IRanges::IRanges(start = 1, end = filtered.ctgs))
    GenomeInfoDb::seqlengths(filtered.ctgs.gr) <- filtered.ctgs
    if (length(filtered.ctgs.gr) > 0) {
      destination <- file.path(datapath, paste0("filteredContigs_", config[['bin.size']], "bp_", config[['bin.method']], ".RData"))
      if (store.data.obj) {
        save(filtered.ctgs.gr, file = destination)
      }
    }
  }  
  
  ## Get statistics on clustered contigs
  all.ctgs <- data.frame(ctg=names(chrom.lengths), ctg.len=chrom.lengths, stringsAsFactors = FALSE, row.names = NULL, index='all.ctgs')
  size.select.ctgs <- data.frame(ctg=names(chrom.lengths)[names(chrom.lengths) %in% chroms.in.data], 
                                 ctg.len=chrom.lengths[names(chrom.lengths) %in% chroms.in.data], 
                                 stringsAsFactors = FALSE, row.names = NULL, 
                                 index = paste0("min.ctg.len.", config[['min.contig.size']]))
  clustered.ctgs <- data.frame(ctg=names(chrom.lengths)[names(chrom.lengths) %in% seqlevels(ordered.contigs.gr)], 
                                ctg.len=chrom.lengths[names(chrom.lengths) %in% seqlevels(ordered.contigs.gr)], 
                                stringsAsFactors = FALSE, row.names = NULL, 
                                index = "clustered.ctgs")
  ctg.stat <- rbind(all.ctgs, size.select.ctgs, clustered.ctgs)
  
  ## Store data object
  destination <- file.path(datapath, paste0("ctgStat_minCtgSize_", config[['min.contig.size']], ".RData"))
  if (store.data.obj) {
    save(ctg.stat, file = destination)
  }
  
  ## Get contig report table ##
  ctgs.report <- data.frame(ctg=names(chrom.lengths), ctg.len=chrom.lengths, stringsAsFactors = FALSE, row.names = NULL)
  ctgs.report$sizeSelect <- FALSE
  ctgs.report$sizeSelect[ctgs.report$ctg %in% size.select.ctgs$ctg] <- TRUE
  ctgs.report$Clustered <- FALSE
  ctgs.report$Clustered[ctgs.report$ctg %in% clustered.ctgs$ctg] <- TRUE
  ctgs.report$putative.error <- FALSE
  if (!is.null(putative.errors.report)) {
    ctgs.report$putative.error[ctgs.report$ctg %in% putative.errors.report$seqnames] <- TRUE
  }  
  ctgs.report$Dir <- NA
  ctgs.report$Dir[match(seqlevels(ordered.contigs.gr), ctgs.report$ctg)] <- ordered.contigs.gr$dir
  ctgs.report$Cluster.ID <- NA
  ctgs.report$Cluster.ID[match(seqlevels(ordered.contigs.gr), ctgs.report$ctg)] <- ordered.contigs.gr$ID
  if (config[['eval.ploidy']]) {
    ctgs.report$Ploidy <- NA
    ctgs.report$Ploidy[match(seqlevels(ordered.contigs.gr), ctgs.report$ctg)] <- ordered.contigs.gr$ploidy
  }
  ctgs.report <- ctgs.report[order(ctgs.report$Cluster.ID, ctgs.report$ctg.len),]
  ## Export summary table of assembly errors
  destination <- file.path(asmpath, paste0("ctgReport_", config[['bin.size']], "bp_", config[['bin.method']], ".tsv"))
  utils::write.table(ctgs.report, file = destination, quote = FALSE, row.names = FALSE, append = FALSE, sep = "\t")
  
  ## Export clustered FASTA ##
  if (!is.null(config[['assembly.fasta']]) & is.character(config[['assembly.fasta']])) {
    if (file.exists(config[['assembly.fasta']])) {
      exportPseudoChromosomalScaffolds(clustered.gr=ordered.contigs.gr, assembly.fasta=config[['assembly.fasta']], outputfolder=asmpath, concat.fasta=config[['concat.fasta']])
    } else {
      message("Submitted 'assembly.fasta' does not exists !!!")
    }  
  }
  
  ## Plot theta parameter
  theta.plt <- plotThetaEstimates(theta.param = EM.obj$theta.param)
  suppressWarnings(
  ggplot2::ggsave(filename = file.path(pltpath, 'theta_param.pdf'), 
         plot = theta.plt, 
         width = config[['num.clusters']] / 10, 
         height = length(EM.obj$theta.param) / 10, 
         limitsize = FALSE)
  )
  
  ## Plot abundance of WC states per initial clusters (n=num.clusters)
  zscore.plt <- plotStrandStateZscore(zscores = theta.zscore)
  suppressWarnings(
  ggplot2::ggsave(filename = file.path(pltpath, 'zscores_wc_states.pdf'), 
                  plot = zscore.plt, 
                  width = 8, 
                  height = 4, 
                  limitsize = FALSE)
  )
  
  ## Plot probability distribution
  prob.plt <- plotEMprobs(em.prob = EM.obj$soft.pVal, prob.th = config[['prob.th']])
  suppressWarnings(
  ggplot2::ggsave(filename = file.path(pltpath, 'prob_dist.pdf'), 
         plot = prob.plt, 
         width = 6, 
         height = 4, 
         limitsize = FALSE)
  )
  
  ## Plot cluster size distribution
  clust.plt <- plotClusteredContigSizes(clustered.gr = ordered.contigs.gr)
  suppressWarnings(
  ggplot2::ggsave(filename = file.path(pltpath, 'cluster_sizes.pdf'), 
         plot = clust.plt, 
         width = 8, 
         height = 4, 
         limitsize = FALSE)
  )
  
  ## Plot statistics of clustered contigs
  ctg.stat.plt <- plotCTGstat(ctg.stat = ctg.stat)
  suppressWarnings(
  ggplot2::ggsave(filename = file.path(pltpath, 'clustered_ctgs.pdf'), 
        plot = ctg.stat.plt, 
        width = 7, 
        height = 5, 
        limitsize = FALSE)
  )
  
  ## Report total processing time
  time <- proc.time() - ptm
  message("\nTotal analysis time: ", round(time[3],2), "s")
}
