log <- file(snakemake@log[[1]], open='wt')
sink(file=log, type='message')
sink(file=log, type='output')

library('scales')
library('ggplot2')
source("postProcessing.R")

allCluster.ranks <- lapply(snakemake@input[["ranks"]], readRDS)
all.prob.trueClust <- lapply(snakemake@input[["true_clust"]], readRDS)

allCluster.ranks.sums <- Reduce("+", allCluster.ranks)
  
#plotting
  
table.ranks.df <- data.frame(rank.acc=names(allCluster.ranks.sums), Freq=allCluster.ranks.sums)
table.ranks.df$rank.acc <- factor( table.ranks.df$rank.acc, levels= table.ranks.df$rank.acc)
ranking.plt <- ggplot(table.ranks.df, aes(x=rank.acc, y=Freq)) + geom_bar(fill="chartreuse4", stat="identity") + xlab("Probability ranking of true cluster") + ylab("Frequency") + theme_bw() + scale_y_continuous(labels=comma)
  
#get accuracy measure (sum probabilities of true cluster divided by all PB reads)
all.prob.trueClust.v <- unlist(all.prob.trueClust)
overal.acc <- sum(all.prob.trueClust.v)/length(all.prob.trueClust.v)

#return(list(ranking.plt=ranking.plt, ranking.table=table.ranks.df,  overal.acc=overal.acc))
ggsave(filename = snakemake@output[["plot"]], plot = ranking.plt)
write.table(table.ranks.df, file = snakemake@output[["rank_table"]])
write.table(overal.acc, file = snakemake@output[["overal_acc"]])
