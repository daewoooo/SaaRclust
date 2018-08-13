log <- file(snakemake@log[[1]], open='wt')
sink(file=log, type='message')
sink(file=log, type='output')

source("postProcessing.R")

clusters <- snakemake@input[["clusters"]]
raw_data <- snakemake@input[["raw_data"]]
v <- boxplotDistSSlibsPerPB(Clusters2process=clusters, Quals2process=raw_data)

saveRDS(v, file = snakemake@output[[1]])
