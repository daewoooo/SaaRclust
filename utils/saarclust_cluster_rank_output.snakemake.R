log <- file(snakemake@log[[1]], open='wt')
sink(file=log, type='message')
sink(file=log, type='output')

source("postProcessing.R")

v <- accuracyRanking(snakemake@input[[1]])

saveRDS(v[["ranks"]], file = snakemake@output[["ranks"]])
saveRDS(v[["trueClust"]], file = snakemake@output[["true_clust"]])
