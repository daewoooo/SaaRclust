library("gsubfn")
library("data.table")
library("RColorBrewer")

path <-  "/media/daewoooo/WORK/SS2PacBio_alignment_HG00733/Test_minimap/minimap_test/" 
files2process <- list.files(path = path, pattern = '\\.log')

results <- list()
plots <- list()
for (file in files2process) {
  message("Working on ",file)
  options <- unlist(strsplit(file, split = "_"))[c(4:7)]

  filePath <- file.path(path, file)
  conn <- file(filePath, open="r")
  lines <-readLines(conn)
  lines.l <- list()
  for (i in 1:length(lines)){
    lines.l[[i]] <- lines[i]
  }
  close(conn)

  minimizers <- strapplyc(lines.l[[4]], "consider: (\\d+)", simplify = T)
  real.time <- strapplyc(lines.l[[7]], "Real time: (\\d+.\\d+)", simplify = T)
  cpu.time <- strapplyc(lines.l[[7]], "CPU: (\\d+.\\d+)", simplify = T)

  k <- strapplyc(options[1], "(\\d+)", simplify = T)
  w <- strapplyc(options[2], "(\\d+)", simplify = T)
  L <- strapplyc(options[3], "(\\d+)", simplify = T)
  f <- strapplyc(options[4], "(\\d+.\\d+)", simplify = T)

  minimap.file <- gsub(file, pattern = '\\.log', replacement = '')
  minimap.file.path <- file.path(path, minimap.file)
  tab.in <- importTestData(infile = minimap.file.path)

  #export data quality measures
  qual.measures <- getQualMeasure(tab.in)
  mapp.stat.counts <- qual.measures$mapp.stat.counts
  mapp.gaps.tab <- qual.measures$mapp.gaps.stat
  
  #qual.plt <- plotQualMeasure(qual.measures)
  #plt.file <- gsub(file, pattern = '\\.log', replacement = '\\_plot.pdf')
  #plt.file <- file.path(path, plt.file)
  #ggsave(filename = plt.file, plot = qual.plt, width = 10, height = 10)

  read.map.dist.mean <- round(mean(mapp.stat.counts$SSread.perPB))
  SSlib.perPB.mean <- round(mean(mapp.stat.counts$SSlib.perPB))
  gaps.perSS.mean <- round(mean(mapp.gaps.tab$matchWithgap))
  accur.counts <- table(mapp.gaps.tab$mapp.accur)
  falses <- round((accur.counts[1]/sum(accur.counts))*100, digits = 2)
  trues <- round((accur.counts[2]/sum(accur.counts))*100, digits = 2)

  total.PBreads <- length(unique(tab.in$PBreadNames))
  total.SSreads <- length(unique(tab.in$SSreadNames))

  v <- c(k=as.numeric(k), w=as.numeric(w), L=as.numeric(L), f=as.numeric(f), num.minimizers=as.numeric(minimizers), real.time=as.numeric(real.time), cpu.time=as.numeric(cpu.time), mean.SSread.perPB=read.map.dist.mean, mean.SSlib.perPB=SSlib.perPB.mean, gaps.perSS.mean=gaps.perSS.mean, total.PBreads=as.numeric(total.PBreads), total.SSreads=as.numeric(total.SSreads), falses=falses, trues=trues)
  results[[file]] <- v
  #plots[[file]] <- qual.plt
}

final.tab <- do.call(rbind,results)
final.tab <- as.data.frame(final.tab)
fileName <- file.path(path, 'results_minimap_test.txt')
write.table(x=final.tab, file = fileName, quote = F)

save(plots, file = "plots_minimap_test.RData")

#final.tab.long <- melt(final.tab, id.vars=c("k","w","L","f"), measure.vars=c("num.minimizers","real.time","cpu.time","mean.SSread.perPB","mean.SSlib.perPB","gaps.perSS.mean","total.PBreads","total.SSreads","falses.FALSE","trues.TRUE"))
final.tab.long <- melt(final.tab, id.vars=c("k","w","L","f"), measure.vars=c("num.minimizers","real.time","mean.SSread.perPB","mean.SSlib.perPB","gaps.perSS.mean","total.PBreads","total.SSreads","falses.FALSE","trues.TRUE"))

ggplot(final.tab, aes(x=f, y=trues.TRUE)) + geom_line() + geom_point() + facet_grid(w ~ .)

final.tab.long$value[is.na(final.tab.long$value)] <- 0
ggplot(final.tab.long, aes(x=f, y=value, group=variable, color=variable)) + geom_line() + geom_point() + facet_grid(variable ~ w + k, scales="free_y") + scale_color_manual(values = brewer.pal(n=9, name="Set1"))

ggplot(final.tab.long[final.tab.long$k != 11,], aes(x=f, y=value, group=variable, color=variable)) + geom_line() + geom_point() + facet_grid(variable ~ w + k, scales="free_y") + scale_color_manual(values = brewer.pal(n=9, name="Set1"))

ggplot(final.tab.long[final.tab.long$k != 11 & final.tab.long$f != 0.5,], aes(x='', y=value, group=variable, color=variable)) + geom_line() + geom_point() + facet_grid(variable ~ w + k + f, scales="free_y") + scale_color_manual(values = brewer.pal(n=9, name="Set1"))

