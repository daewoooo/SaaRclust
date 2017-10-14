#' Plot data quality measures.
#'
#' Takes imported data table and plots some relevant data quality measures.
#'
#' @param summary.tab Imported data table.
#' @author David Porubsky
#' @export

plotQualMeasure <- function(summary.tab) {

  mapp.dist.tab <- summary.tab$mapp.stat.counts
  mapp.gaps.tab <- summary.tab$mapp.gaps.stat
  SScov.stat.m <- summary.tab$SScov.stat          
  ord <- order(match(rownames(SScov.stat.m), mapp.dist.tab$PBreadNames))
  SScov.stat.m <- SScov.stat.m[ord,]
  #mapp.SScov.stat <- summary.tab$mapp.SScov.stat
  
  
  #plotting distribution of SSreads mapped to PB reads
  read.map.dist.mean <- round(mean(mapp.dist.tab$SSread.perPB))
  #read.dist.plt <- ggplot(mapp.dist.tab) + geom_linerange(aes(x=c(1:nrow(mapp.dist.tab)),ymin=0, ymax=SSread.perPB), stat="identity") + theme_bw() + geom_hline(yintercept = read.map.dist.mean, color="red") + geom_text(aes(nrow(mapp.dist.tab), read.map.dist.mean, label = read.map.dist.mean, vjust = -1), color="red") + xlab("Sorted PB reads") + ylab("# of mapped SS reads")
  #read.dist.plt.log <- ggplot(mapp.dist.tab) + geom_linerange(aes(x=c(1:nrow(mapp.dist.tab)),ymin=0, ymax=SSread.perPB), stat="identity") + theme_bw() + geom_hline(yintercept = read.map.dist.mean, color="red") + geom_text(aes(nrow(mapp.dist.tab), read.map.dist.mean, label = read.map.dist.mean, vjust = -1), color="red") + xlab("Sorted PB reads") + ylab("# of mapped SS reads (log10)") + scale_y_continuous(trans = "log10")
  read.dist.plt <- ggplot(mapp.dist.tab) + geom_line(aes(x=c(1:nrow(mapp.dist.tab)),y=SSread.perPB)) + theme_bw() + geom_hline(yintercept = read.map.dist.mean, color="red") + geom_text(aes(nrow(mapp.dist.tab), read.map.dist.mean, label = read.map.dist.mean, vjust = -1), color="red") + xlab("Sorted PB reads") + ylab("# of mapped SS reads")
  read.dist.plt.log <- ggplot(mapp.dist.tab) + geom_line(aes(x=c(1:nrow(mapp.dist.tab)),y=SSread.perPB)) + theme_bw() + geom_hline(yintercept = read.map.dist.mean, color="red") + geom_text(aes(nrow(mapp.dist.tab), read.map.dist.mean, label = read.map.dist.mean, vjust = -1), color="red") + xlab("Sorted PB reads") + ylab("# of mapped SS reads (log10)") + scale_y_continuous(trans = "log10")
  
  #plotting distribution of SSlibs represented per PB read
  SSlib.perPB.mean <- round(mean(mapp.dist.tab$SSlib.perPB))
  #SSlib.perPB.plt <- ggplot(mapp.dist.tab) + geom_linerange(aes(x=c(1:nrow(mapp.dist.tab)),ymin=0, ymax=SSlib.perPB), stat="identity") + theme_bw() + geom_hline(yintercept = SSlib.perPB.mean, color="red") + geom_text(aes(nrow(mapp.dist.tab), SSlib.perPB.mean, label = SSlib.perPB.mean, vjust = -1), color="red") + xlab("PB reads sorted by the number\nof SSreads per PBread") + ylab("# of SS libs per PB read")
  #SSlib.perPB.plt.log <- ggplot(mapp.dist.tab) + geom_linerange(aes(x=c(1:nrow(mapp.dist.tab)),ymin=0, ymax=SSlib.perPB), stat="identity") + theme_bw() + geom_hline(yintercept = SSlib.perPB.mean, color="red") + geom_text(aes(nrow(mapp.dist.tab), SSlib.perPB.mean, label = SSlib.perPB.mean, vjust = -1), color="red") + xlab("PB reads sorted by the number\nof SSreads per PBread") + ylab("# of SS libs per PB read (log10)") + scale_y_continuous(trans = "log10")
  SSlib.perPB.plt <- ggplot(mapp.dist.tab) + geom_line(aes(x=c(1:nrow(mapp.dist.tab)),y=SSlib.perPB)) + theme_bw() + geom_hline(yintercept = SSlib.perPB.mean, color="red") + geom_text(aes(nrow(mapp.dist.tab), SSlib.perPB.mean, label = SSlib.perPB.mean, vjust = -1), color="red") + xlab("PB reads sorted by the number\nof SSreads per PBread") + ylab("# of SS libs per PB read")
  SSlib.perPB.plt.log <- ggplot(mapp.dist.tab) + geom_line(aes(x=c(1:nrow(mapp.dist.tab)),y=SSlib.perPB)) + theme_bw() + geom_hline(yintercept = SSlib.perPB.mean, color="red") + geom_text(aes(nrow(mapp.dist.tab), SSlib.perPB.mean, label = SSlib.perPB.mean, vjust = -1), color="red") + xlab("PB reads sorted by the number\nof SSreads per PBread") + ylab("# of SS libs per PB read (log10)") + scale_y_continuous(trans = "log10")
  
  
  #plotting sum of gaps per PB read normalized by SSread counts mapped to a given PB read
  gaps.perPB.top1perc <- round(quantile(mapp.dist.tab$gaps.perPB.norm,prob=1-1/100)) #99th quantile
  #gaps.perPB.plt <- ggplot(mapp.dist.tab) + geom_linerange(aes(x=c(1:nrow(mapp.dist.tab)),ymin=0, ymax=gaps.perPB.norm), stat="identity") + theme_bw() + geom_hline(yintercept = gaps.perPB.top1perc, color="red") + geom_text(aes(nrow(mapp.dist.tab), gaps.perPB.top1perc, label = gaps.perPB.top1perc, vjust = -1), color="red") + xlab("PB reads sorted by the number\nof SSreads per PBread") + ylab("Normalized sum of gaps per PB read")
  gaps.perPB.plt <- ggplot(mapp.dist.tab) + geom_line(aes(x=c(1:nrow(mapp.dist.tab)),y=gaps.perPB.norm)) + theme_bw() + geom_hline(yintercept = gaps.perPB.top1perc, color="red") + geom_text(aes(nrow(mapp.dist.tab), gaps.perPB.top1perc, label = gaps.perPB.top1perc, vjust = -1), color="red") + xlab("PB reads sorted by the number\nof SSreads per PBread") + ylab("Normalized sum of gaps per PB read")
  
  #plotting counts of SS reads per lib per PB read
  mapp.SScov.stat.mean <- mean(SScov.stat.m[SScov.stat.m>0])
  SScov.stat.m[SScov.stat.m==0] <- NA
  SScov.stat.df <- data.frame(value=rowMeans(SScov.stat.m, na.rm = T))
  SScov.stat.plt <- ggplot(SScov.stat.df) + geom_line(aes(x=c(1:nrow(SScov.stat.df)),y=value)) + theme_bw() + geom_hline(yintercept = mapp.SScov.stat.mean, color="red") + geom_text(aes(nrow(SScov.stat.df), mapp.SScov.stat.mean, label = mapp.SScov.stat.mean, vjust = -1), color="red") + xlab("PB reads sorted by the mean number\nof SSreads per SSlib per PBread") + ylab("Mean counts of SSreads per SSlib per PB read")
  #hm <- Heatmap(SScov.stat.m, name = "SScounts", cluster_columns = F, cluster_rows = F, show_row_names = FALSE)
  #hm <- Heatmap(test, name = "SScounts", cluster_columns = F, cluster_rows = F, show_row_names = FALSE, col = colorRamp2(c(1,2,3,4,5), topo.colors(n=5)))
  #Heatmap(SScov.stat.m, name = "SScounts", cluster_columns = F, cluster_rows = F, show_row_names = FALSE, col = colorRamp2(c(0,1,2,3), c("white","blue", "green","red")))
  
  #mapp.SScov.stat.mean <- mean(mapp.SScov.stat)
  #mapp.SScov.stat.df <- data.frame(value=sort(mapp.SScov.stat, decreasing = T))
  #mapp.SScov.stat.plt <- ggplot(mapp.SScov.stat.df) + geom_line(aes(x=c(1:nrow(mapp.SScov.stat.df)),y=value)) + theme_bw() + geom_hline(yintercept = mapp.SScov.stat.mean, color="red") + geom_text(aes(nrow(mapp.SScov.stat.df), mapp.SScov.stat.mean, label = mapp.SScov.stat.mean, vjust = -1), color="red") + xlab("PB reads sorted by the mean number\nof SSreads per SSlib per PBread") + ylab("Counts of SSreads per SSlib per PB read")
  
  #plotting match with gaps per SS read flagged TRUE or FALSE based on mapping agreement
  gaps.perSS.mean <- round(mean(mapp.gaps.tab$matchWithgaps))
  accur.counts <- table(mapp.gaps.tab$mapp.accur)
  falses <- round((accur.counts[1]/sum(accur.counts))*100, digits = 2)
  trues <- round((accur.counts[2]/sum(accur.counts))*100, digits = 2)
  false.lab <- paste('MISS ', falses, '%', sep = "")
  true.lab <- paste('MATCH ', trues, '%', sep = "")
  gaps.perSS.plt <- ggplot(mapp.gaps.tab) + geom_linerange(aes(x=c(1:nrow(mapp.gaps.tab)),ymin=0, ymax=matchWithgaps, color=mapp.accur)) + theme_bw() + geom_hline(yintercept = gaps.perSS.mean, color="red") + geom_text(aes(nrow(mapp.gaps.tab), gaps.perSS.mean, label = gaps.perSS.mean, vjust = -1), color="red") + xlab("Sorted SS reads by the length\nof SS read alignment (bp)") + ylab("(bp) SS read alignment with gaps (log10)") + scale_y_continuous(trans = "log10") + scale_color_manual(labels = c(false.lab, true.lab), values = c("red", "green"))
 
  #gaps.rle <- rle(mapp.gaps.tab$matchWithgaps)
  #gaps.df <- data.frame(gaps=gaps.rle$values, mapp.accur=mapp.gaps.tab$mapp.accur[which(gaps.rle$values %in% mapp.gaps.tab$matchWithgaps)])
  #gaps.perSS.plt <- ggplot(gaps.df) + geom_linerange(aes(x=c(1:nrow(gaps.df)), ymin=0, ymax=gaps, color=mapp.accur)) + theme_bw() + geom_hline(yintercept = gaps.perSS.mean, color="red") + geom_text(aes(nrow(gaps.df), gaps.perSS.mean, label = gaps.perSS.mean, vjust = -1), color="red") + xlab("Sorted SS reads by the length\nof SS read alignment (bp)") + ylab("(bp) SS read alignment with gaps (log10)") + scale_y_continuous(trans = "log10") + scale_color_manual(labels = c(false.lab, true.lab), values = c("red", "green"))
  
  #suppressWarnings( plt <- plot_grid(read.dist.plt, read.dist.plt.log, SSlib.perPB.plt, SSlib.perPB.plt.log, gaps.perPB.plt, gaps.perSS.plt, ncol = 2) )
  suppressWarnings( plt <- plot_grid(read.dist.plt, read.dist.plt.log, SSlib.perPB.plt, SSlib.perPB.plt.log, SScov.stat.plt, gaps.perPB.plt, ncol = 2) )
  return(plt)
}