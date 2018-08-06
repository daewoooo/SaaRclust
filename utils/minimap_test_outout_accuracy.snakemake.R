log <- file(snakemake@log[[1]], open='wt')
sink(file=log, type='message')
sink(file=log, type='output')

source("processMinimapTest.R")

file <- snakemake@input[[1]]
v <- getAccuracy(file)

write.table(v, file = snakemake@output[[1]], quote = F, sep = "\t", row.names = F)
