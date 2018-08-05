## Load required libraries
library("gsubfn")
library("data.table")
library("RColorBrewer")
library("SaaRclust")

## Path to the test minimap folder to be processed
path <- "/home/porubsky/WORK/PROJECTS/SaaRclust_project/TEST_minimap/aligns_50000reads" 
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
  tab.in <- importData(infile = minimap.file.path)

  ## Export data quality measures
  qual.measures <- getQualMeasure(tab.in)

  mean.SSreads.perPB <- round(mean(qual.measures$SSreads.perPB))
  mean.SSlib.perPB <- round(mean(qual.measures$SSlib.perPB$counts))
  mean.SSreads.perlib.perPB <- round(mean(qual.measures$SSreads.perlib.perPB))
  accur.counts <- table(tab.in$SSchrom == tab.in$PBchrom)
  falses <- round((accur.counts[1]/sum(accur.counts))*100, digits = 2)
  trues <- round((accur.counts[2]/sum(accur.counts))*100, digits = 2)

  total.PBreads <- length(unique(tab.in$PBreadNames))
  total.SSreads <- length(unique(tab.in$SSreadNames))

  v <- c(k=as.numeric(k), w=as.numeric(w), L=as.numeric(L), f=as.numeric(f), num.minimizers=as.numeric(minimizers), real.time=as.numeric(real.time), cpu.time=as.numeric(cpu.time), mean.SSreads.perPB=mean.SSreads.perPB, mean.SSlib.perPB=mean.SSlib.perPB, mean.SSreads.perlib.perPB=mean.SSreads.perlib.perPB, total.PBreads=as.numeric(total.PBreads), total.SSreads=as.numeric(total.SSreads), falses=falses, trues=trues)
  results[[file]] <- v
}

final.tab <- do.call(rbind,results)
final.tab <- as.data.frame(final.tab)
fileName <- file.path(path, 'results_minimap_test.txt')
write.table(x=final.tab, file = fileName, quote = F)

save(plots, file = "plots_minimap_test.RData")

## Prepare table for plotting
final.tab.long <- melt(final.tab, id.vars=c("k","w","L","f"), measure.vars=c("num.minimizers","real.time","mean.SSreads.perPB","mean.SSlib.perPB","mean.SSreads.perlib.perPB","total.PBreads","total.SSreads","falses.FALSE","trues.TRUE"))

## Prepare plot
plt <- ggplot(final.tab.long, aes(x='', y=value, color=variable))
plt <- plt + geom_point(size=2) + facet_grid(variable ~ L + k + f, scales="free_y", labeller = labeller(.cols = label_both))
plt <- plt + scale_color_manual(values = brewer.pal(n=9, name="Set1"), guide="none")
plt <- plt + theme(strip.text.y = element_text(angle = 360))
plt <- plt + xlab("")

