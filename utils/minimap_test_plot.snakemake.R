log <- file(snakemake@log[[1]], open='wt')
sink(file=log, type='message')
sink(file=log, type='output')

library("reshape")
library("ggplot2")
library("data.table")
library("RColorBrewer")

#source("processMinimapTest.R")

results <- lapply(snakemake@input, function(x) read.table(x, header = T))
final.tab <- do.call(rbind, results)
final.tab <- as.data.frame(final.tab)

final.tab.long <- melt(final.tab, id.vars=c("k","w","L","f"), measure.vars=c("num.minimizers","real.time","mean.SSreads.perPB","mean.SSlib.perPB","mean.SSreads.perlib.perPB","total.PBreads","total.SSreads","falses","trues"))

## Prepare plot
plt <- ggplot(final.tab.long, aes(x='', y=value, color=variable))
plt <- plt + geom_point(size=2) + facet_grid(variable ~ L + k + f, scales="free_y", labeller = labeller(.cols = label_both))
plt <- plt + scale_color_manual(values = brewer.pal(n=9, name="Set1"), guide="none")
plt <- plt + theme(strip.text.y = element_text(angle = 360))
plt <- plt + xlab("")

ggsave(filename = snakemake@output[[1]], plot = plt)


#read.table(v, file = snakemake@input[[1]], header = T)
