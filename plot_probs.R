library(ggplot2) #requires ggplot2 2.0
library(reshape2)

files_pp <- dir()[grep("postProbSingle", dir())]
ref_file <- dir()[grep("Ref", dir(), ignore.case=T)]

name <- read.table(ref_file, nrows=1, stringsAsFactors = FALSE, sep="\t")[1,]
name <- substr(name,2,nchar(name))
sequence <- unlist(strsplit(read.table(ref_file, nrows=1, skip=1, stringsAsFactors = FALSE)[1,1], split=""))
structure <- unlist(strsplit(read.table(ref_file, nrows=1, skip=2, stringsAsFactors = FALSE)[1,1], split=""))
structure_col <- structure
structure_col[which(structure %in% ".")] <- "unpaired"
structure_col[which(structure %in% c("(",")"))] <- "paired"

postProb_CMCT <- 1-read.table(files_pp[grep("-CMCT", files_pp, fixed=T)][1])
postProb_DMS <- 1-read.table(files_pp[grep("-DMS", files_pp, fixed=T)][1])
postProb_SeqOnly <- 1-read.table(files_pp[grep("Seq_only", files_pp, fixed=T)][1])
postProb_SHAPE <- 1-read.table(files_pp[grep("-SHAPE", files_pp, fixed=T)][2])
postProb_All <-  1-read.table(files_pp[grep("SHAPE-DMS-CMCT", files_pp, fixed=T)][1])

rectangles <- (1:length(structure))[structure_col=="paired"]
xmin_new <- rectangles[1]
xmax_new <- NULL
previous <- rectangles[1]
for(i in 2:length(rectangles)) {
  current <- rectangles[i]
  if (previous+1 != current) {
    xmax_new <- c(xmax_new,previous)
    xmin_new <- c(xmin_new,current)
  }
  previous <- current
  if (i == length(rectangles) && length(xmax_new)!=length(xmin_new)) xmax_new <- c(xmax_new,current)
}
rectangles <- data.frame(xmin = xmin_new-0.34, xmax = xmax_new+0.34)

plotter <- data.frame(Position=1:length(structure) , Structure=structure, Known_structure=structure_col, Sequence=sequence, postProb_SeqOnly=postProb_SeqOnly, postProb_CMCT=postProb_CMCT, postProb_DMS=postProb_DMS, postProb_SHAPE=postProb_SHAPE, postProb_All=postProb_All)
plotter <- melt(plotter, id.vars = c("Position", "Structure","Known_structure", "Sequence"))
plotter$variable <- factor(plotter$variable, labels=c("Sequence only", "Sequence and CMCT", "Sequence and DMS",  "Sequence and SHAPE", "Sequence and all probing sets"))
colour_palette <- c("darkblue","darkgrey")
pdf(file=paste(name,"paired_pp_barplots.pdf",sep="_"),width=11.7,height=8.27)
ggplot(plotter) + geom_bar(stat="identity", aes(x=Position, y=value, fill=Known_structure), alpha=0.9, width=0.75) + geom_hline(yintercept=0.5) + facet_wrap(~variable, ncol=1) + theme_bw() + theme(legend.position="none", panel.grid.minor=element_blank(), strip.background=element_rect(colour = "white")) + ylab("Posterior probability of pairing") + scale_x_discrete(breaks=1:length(structure), labels=structure, name="Known structure") + scale_y_continuous(limits=c(0,1)) + ggtitle(name) + scale_fill_manual(values=colour_palette) + geom_rect(data=rectangles, aes(ymin=-Inf, ymax=Inf, xmin=xmin, xmax=xmax), fill='gray80', alpha=0.25)
dev.off()
cat("done plotting")