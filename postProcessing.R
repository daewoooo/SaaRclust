############################
### data post-processing ###
############################

#plot likelihood function
logl.df <- data.frame(log.l=soft.clust$log.l, iter=1:length(soft.clust$log.l))
logl.diff <- round(max(logl.df$log.l) - min(logl.df$log.l), digits = 0)
p1 <- ggplot(logl.df) + geom_line(aes(x=iter, y=log.l), color="red") + xlab("EM iterations") + ylab("Likelihood function") + annotate("text",  x=Inf, y = Inf, label = paste("Diff: ",logl.diff), vjust=1, hjust=1)

#plot cluster accuracy
p2 <- plotClustAccuracy(pVal.df = soft.clust.df, num.clusters = num.clusters)

p12 <- plot_grid(p1, p2, nrow = 1)

#plot pi values
pi.param.df <- data.frame(pi.param=soft.clust$pi.param, id=1:length(soft.clust$pi.param))
p3 <- ggplot(pi.param.df , aes(x=reorder(id, -pi.param), y=pi.param)) + geom_bar(stat='identity') + xlab("Estimated cluster sizes")
#NOTE: overlay with known chromosome sizes

#plot heatmap values
#p4 <- plotHeatmap(pVal.df = soft.clust.df, colOrder=NULL, num.clusters = num.clusters)

p4 <- plotThetaEstimates(theta.param = soft.clust$theta.param, title = "No WC constrain")
p5 <- plotThetaEstimates(theta.param = soft.clust.rescale$theta.param, title = "Applied WC constrain")

final.plot <- plot_grid(p12, p3, p4, p5, nrow = 4, rel_heights = c(1,1,4,4))


#NOTE check representation of SSreads in properly clustered PB reads
pVals <- soft.clust.df[,c(1:num.clusters)]
ord <- apply(pVals, 1, which.max)
soft.clust.df$clustID <- ord 
max.pVal <- apply(pVals, 1, max)
mask <- max.pVal >= 0.9
postpric.tab <- soft.clust.df[mask,]
PBreadsClust <- split(postpric.tab$PBreadNames, postpric.tab$clustID)

test <- tab.filt[tab.filt$PBreadNames %in% PBreadsClust[[1]],]

#get cluster sizes
#library("BSgenome.Hsapiens.UCSC.hg38")
#seq.len <- seqlengths(BSgenome.Hsapiens.UCSC.hg38)[paste0('chr', c(1:22))]
#seq.len <- as.numeric(seq.len)
#seq.len.ratio <- seq.len/sum(seq.len)
#seq.len.ratio <- rep(seq.len.ratio, each=2)
#cluster.sizes <- round((chunk.size/2)*seq.len.ratio)
#cluster.sizes <- rep(cluster.sizes, each=2)
