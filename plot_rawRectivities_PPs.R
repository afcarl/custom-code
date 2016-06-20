library(ggplot2)
library(reshape2)

files_pp <- dir()[grep("postProbSingle", dir())]
files_predictions <- dir()[grep("dat", dir())]
ref_file <- dir()[grep("Ref", dir(), ignore.case=T)]
files_raw <- dir()[grep("raw", dir())]

name <- read.table(ref_file, nrows=1, stringsAsFactors = FALSE, sep="\t")[1,]
name <- substr(name,2,nchar(name))
sequence <- unlist(strsplit(read.table(ref_file, nrows=1, skip=1, stringsAsFactors = FALSE)[1,1], split=""))[6:90]
structure <- unlist(strsplit(read.table(ref_file, nrows=1, skip=2, stringsAsFactors = FALSE)[1,1], split=""))[6:90]
structure_col <- structure
structure_col[which(structure %in% ".")] <- "unpaired"
structure_col[which(structure %in% c("(",")"))] <- "paired"

raw_CMCT <- as.numeric(read.table(files_raw[grep("CMCT", files_raw)])[,2])[6:90]
raw_CMCT[is.na(raw_CMCT)] <- 0
raw_DMS <- as.numeric(read.table(files_raw[grep("DMS", files_raw)])[,2])[6:90]
raw_DMS[is.na(raw_DMS)] <- 0
raw_SHAPE <- as.numeric(read.table(files_raw[grep("SHAPE", files_raw)])[,2])[6:90]
raw_SHAPE[is.na(raw_SHAPE)] <- 0

plotter_raw <- data.frame(Position=1:length(structure) , Structure=structure, Known_structure=structure_col, raw_CMCT=raw_CMCT, raw_DMS=raw_DMS, raw_SHAPE=raw_SHAPE)
plotter_raw <- melt(plotter_raw, id.vars = c("Position", "Structure","Known_structure"))
plotter_raw$variable <- factor(plotter_raw$variable, labels=c("CMCT", "DMS", "SHAPE"))

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

colour_palette <- c("darkblue","darkgrey")
pdf(file=paste(name,"raw_reactivity_barplots_v2.pdf",sep="_"),width=7,height=5)
ggplot(plotter_raw) + geom_bar(stat="identity", aes(x=Position, y=value, fill=Known_structure), alpha=0.9, width=0.75) + facet_wrap(~variable, ncol=1, scales = "free_y") + theme_bw() + theme(legend.position="none", panel.grid=element_blank(), strip.text=element_blank()) + ylab("Raw reactivity") + scale_x_discrete(breaks=1:length(structure), labels=structure, name="Known structure") + ggtitle(name) + scale_fill_manual(values=colour_palette) + geom_rect(data=rectangles, aes(ymin=-Inf, ymax=Inf, xmin=xmin, xmax=xmax), fill='gray80', alpha=0.25)
dev.off()

postProb_CMCT <- 1-read.table(files_pp[grep("-CMCT", files_pp, fixed=T)][1])[6:90,1]
structure_CMCT <- unlist(strsplit(read.table(files_predictions[grep("CMCT.", files_predictions, fixed=T)][1], nrows=1, skip=2, stringsAsFactors = FALSE)[1,1], split=""))[6:90]
#structure_CMCT[which(structure_CMCT %in% ".")] <- "unpaired"
#structure_CMCT[which(structure_CMCT %in% c("(",")"))] <- "paired"

postProb_DMS <- 1-read.table(files_pp[grep("-DMS", files_pp, fixed=T)][1])[6:90,1]
structure_DMS <- unlist(strsplit(read.table(files_predictions[grep("DMS.", files_predictions, fixed=T)][1], nrows=1, skip=2, stringsAsFactors = FALSE)[1,1], split=""))[6:90]
#structure_DMS[which(structure_DMS %in% ".")] <- "unpaired"
#structure_DMS[which(structure_DMS %in% c("(",")"))] <- "paired"

postProb_SeqOnly <- 1-read.table(files_pp[grep("Seq_only", files_pp, fixed=T)][1])[6:90,1]
structure_SeqOnly <- unlist(strsplit(read.table(files_predictions[grep("SeqOnly.", files_predictions, fixed=T)][1], nrows=1, skip=2, stringsAsFactors = FALSE)[1,1], split=""))[6:90]
#structure_SeqOnly[which(structure_SeqOnly %in% ".")] <- "unpaired"
#structure_SeqOnly[which(structure_SeqOnly %in% c("(",")"))] <- "paired"

postProb_SHAPE <- 1-read.table(files_pp[grep("-SHAPE", files_pp, fixed=T)][2])[6:90,1]
structure_SHAPE <- unlist(strsplit(read.table(files_predictions[grep("SHAPE.", files_predictions, fixed=T)][1], nrows=1, skip=2, stringsAsFactors = FALSE)[1,1], split=""))[6:90]
#structure_SHAPE[which(structure_SHAPE %in% ".")] <- "unpaired"
#structure_SHAPE[which(structure_SHAPE %in% c("(",")"))] <- "paired"

postProb_All <-  1-read.table(files_pp[grep("SHAPE-DMS-CMCT", files_pp, fixed=T)][1])[6:90,1]
structure_All <- unlist(strsplit(read.table(files_predictions[grep("SHAPE-DMS-CMCT", files_predictions, fixed=T)][1], nrows=1, skip=2, stringsAsFactors = FALSE)[1,1], split=""))[6:90]
#structure_All[which(structure_All %in% ".")] <- "unpaired"
#structure_All[which(structure_All %in% c("(",")"))] <- "paired"

plotter <- data.frame(Position=1:length(structure) , Structure=structure, Known_structure=structure_col, postProb_SeqOnly=postProb_SeqOnly, postProb_CMCT=postProb_CMCT, postProb_DMS=postProb_DMS, postProb_SHAPE=postProb_SHAPE, postProb_All=postProb_All)
plotter <- melt(plotter, id.vars = c("Position", "Structure","Known_structure"))
plotter$variable <- factor(plotter$variable, labels=c("Sequence only", "Sequence and CMCT", "Sequence and DMS",  "Sequence and SHAPE", "Sequence and all probing sets"))
#plotter$Predicted_structure <- c(structure_SeqOnly, structure_CMCT, structure_DMS, structure_SHAPE, structure_All)
colour_palette <- c("darkblue","darkgrey")
pdf(file=paste(name,"pp_barplots_v2.pdf",sep="_"),width=7,height=8.27)
ggplot(plotter) + geom_bar(stat="identity", aes(x=Position, y=value, fill=Known_structure), alpha=0.9, width=0.75) + geom_hline(yintercept=0.5) + facet_wrap(~variable, ncol=1) + theme_bw() + theme(legend.position="none", panel.grid=element_blank(), strip.text=element_blank()) + ylab("Posterior probability of pairing") + scale_x_discrete(breaks=1:length(structure), labels=structure, name="Known structure") + scale_y_continuous(limits=c(0,1)) + ggtitle(name)+ scale_fill_manual(values=colour_palette) + geom_rect(data=rectangles, aes(ymin=-Inf, ymax=Inf, xmin=xmin, xmax=xmax), fill='gray80', alpha=0.25)
dev.off()
