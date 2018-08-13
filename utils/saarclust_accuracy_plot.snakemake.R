log <- file(snakemake@log[[1]], open='wt')
sink(file=log, type='message')
sink(file=log, type='output')

library(ggplot2)

results <- lapply(snakemake@input, function(x) read.table(x, header = T))
clust.acc.df <- Reduce("+", results)
  
#calcualte accuracy percentages
clust.acc.df$prob.th <- results[[1]]$prob.th
clust.acc.df$th.acc <- clust.acc.df$acc.th.match / clust.acc.df$acc.th.sum
clust.acc.df$th.clustReads <- clust.acc.df$acc.th.sum / clust.acc.df$allReads

#get genome size
suppressMessages( library("biovizBase") )
hg38Ideogram <- getIdeogram("hg38", cytoband = FALSE)
hg38Ideogram <- keepSeqlevels(hg38Ideogram, paste0('chr', c(1:22,'X')), pruning.mode="coarse")
genome.size <- sum(as.numeric(seqlengths(hg38Ideogram)))
clust.acc.df$depth <- ceiling(clust.acc.df$seq.bases/genome.size)

acc.plt <- ggplot(clust.acc.df) + geom_point(aes(x=th.acc, y=th.clustReads), color="deepskyblue4", size=10) + geom_linerange(aes(ymin=-Inf, x=th.acc, ymax=th.clustReads),color="deepskyblue4") + scale_y_continuous(limits = c(0.8, 1)) + scale_x_continuous(limits = c(0,1)) + scale_y_continuous(limits = c(0,1)) + ylab("(%) evaluated PB reads") + xlab("(%) correctly assigned PB reads") + geom_text(aes(x=th.acc, y=th.clustReads), label=c('all', clust.acc.df$prob.th[-1]), color="white") + theme_bw()  + theme_bw() + geom_text(aes(x=th.acc, y=th.clustReads+0.05), label=paste0(clust.acc.df$depth, "x"), color="black")
 
write.table(clust.acc.df, file=snakemake@output[["accuracy_table"]], quote=F, sep="\t", row.names = F)
ggsave(filename = snakemake@output[["accuracy_plot"]], plot = acc.plt)
