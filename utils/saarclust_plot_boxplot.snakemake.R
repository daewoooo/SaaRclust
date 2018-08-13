log <- file(snakemake@log[[1]], open='wt')
sink(file=log, type='message')
sink(file=log, type='output')

library(ggplot2)

allSSlib.perPB <- lapply(snakemake@input, readRDS)
print(class(allSSlib.perPB))
print(allSSlib.perPB)

unlist( sapply(allSSlib.perPB, function(x) x['trues']), use.names = F ) -> trues
unlist( sapply(allSSlib.perPB, function(x) x['falses']), use.names = F )-> falses
ID <- rep(c('Correct', 'Incorrect'), c(length(trues), length(falses)))
counts <- c(trues, falses)
allSSlib.perPB.df <- data.frame(counts=counts, ID=ID)
#Plot 
box.plt <- ggplot(allSSlib.perPB.df, aes(x=ID, y=counts, fill=ID)) + geom_boxplot(outlier.colour="red") + scale_fill_manual(values = c("darkolivegreen3" ,"darkgoldenrod1"), guide="none") + xlab("") + ylab("# of Strand-seq libraries per PB read") + theme_bw()

write.table(allSSlib.perPB.df, file=snakemake@output[["boxplot_table"]], quote=F, sep="\t", row.names = F)
ggsave(filename = snakemake@output[["boxplot"]], plot = box.plt)
