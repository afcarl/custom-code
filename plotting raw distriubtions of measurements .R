require(ggplot2)
require(gridExtra)
as.numeric.factor <- function(x) {as.numeric(levels(x))[x]}
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  require(grid)

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)), ncol = cols, nrow = ceiling(numPlots/cols))
  }

 if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row, layout.pos.col = matchidx$col))
    }
  }
}

G1_pred <- G1[56:82]
G2_pred <- G2[487:730]
G1_train <- G1[1:55]
G2_train <- G2[1:487]

# training data alone
pdf(file="top10_sota_train.pdf",width=11.7,height=8.27)
colour_palette <- c("#56B4E9","#D55E00")
for (i in 1:10) {
	cand <- as.character(top100_sota[i])
	#pdf(file=paste(cand,".pdf",sep=""),width=11.7,height=8.27)
	
	temp <- data.frame(PROMOTER = mmatrix_BRCA_PROMOTER[cand,c(G1_train,G2_train)], BODY=mmatrix_BRCA_BODY[cand,c(G1_train,G2_train)], CPM=t(cpm[cand,c(G1_train,G2_train)]), sampleType=c(rep("AN",length(G1_train)),rep("T",length(G2_train))))
	colnames(temp)[3] <- "CPM"
	
	plot1 <- qplot(temp$sampleType,temp$CPM,geom="boxplot", colour=temp$sampleType) + theme_bw() + theme(panel.grid=element_line(linetype = 0), legend.position="none", plot.title=element_text(size=10)) + scale_y_log10() +xlab("") +ylab("Expression (reads per million)") + ggtitle(paste(cand,"; p-value=",signif(results_all[cand,1],digits=3),", edgeR",sep="")) + scale_colour_manual(values=colour_palette)
	plot2 <- qplot(temp$sampleType,temp$PROMOTER,geom="boxplot", colour=temp$sampleType) + theme_bw() + theme(panel.grid=element_line(linetype = 0), legend.position="none", plot.title=element_text(size=10)) +xlab("") +ylab("Pr. meth. (M-value)") + ggtitle(paste(cand,"; p-value=",signif(results_all[cand,2],digits=3),", Welch's t-test",sep="")) + scale_colour_manual(values=colour_palette)
	plot3 <- qplot(temp$sampleType,temp$BODY,geom="boxplot", colour=temp$sampleType) + theme_bw() + theme(panel.grid=element_line(linetype = 0), legend.position="none", plot.title=element_text(size=10)) +xlab("") +ylab("GB. meth. (M-value)") + ggtitle(paste(cand,"; p-value=",signif(results_all[cand,3],digits=3),", Welch's t-test",sep="")) + scale_colour_manual(values=colour_palette)
	
	correlation_an <- signif(cor(temp[G1_train,1], temp[G1_train,2], method="spearman"),digits=2)
	correlation_t <- signif(cor(temp[G2_train,1], temp[G2_train,2], method="spearman"),digits=2)
	plot4 <- ggplot(temp,aes(x=PROMOTER,y=BODY,colour=sampleType)) +theme_bw() + theme(panel.grid=element_line(linetype = 0), plot.title=element_text(size=8)) + geom_point(alpha=0.25) + ggtitle(paste(cand,"; cor(AN's)=", correlation_an,", cor(T's)=", correlation_t,", Spearman's rho",sep="")) + xlab("Pr. meth. (M-value)") + ylab("GB. meth. (M-value)") + stat_density2d(alpha=0.5) + theme(legend.position="bottom") + scale_colour_manual(values=colour_palette)
	
	correlation_an <- signif(cor(temp[G1_train,1], temp[G1_train,3], method="spearman"),digits=2)
	correlation_t <- signif(cor(temp[G2_train,1], temp[G2_train,3], method="spearman"),digits=2)
	plot5 <- ggplot(temp,aes(x=PROMOTER,y=CPM,colour=sampleType)) + theme_bw() + theme(panel.grid=element_line(linetype = 0), plot.title=element_text(size=8)) + scale_y_log10() + geom_point(alpha=0.25) + ggtitle(paste(cand,"; cor(AN's)=", correlation_an,", cor(T's)=", correlation_t,", Spearman's rho",sep="")) + xlab("Pr. meth. (M-value)") + ylab("Expression (reads per million)") + stat_density2d(alpha=0.5) + theme(legend.position="bottom") + scale_colour_manual(values=colour_palette)
	
	correlation_an <- signif(cor(temp[G1_train,2], temp[G1_train,3], method="spearman"),digits=2)
	correlation_t <- signif(cor(temp[G2_train,2], temp[G2_train,3], method="spearman"),digits=2)
	plot6 <- ggplot(temp,aes(x=BODY,y=CPM,colour=sampleType)) + scale_y_log10() + theme_bw() + theme(panel.grid=element_line(linetype = 0), plot.title=element_text(size=8)) + geom_point(alpha=0.25) + ggtitle(paste(cand,"; cor(AN's)=", correlation_an,", cor(T's)=", correlation_t,", Spearman's rho",sep="")) + xlab("GB. meth. (M-value)") + ylab("Expression (reads per million)") + stat_density2d(alpha=0.5) + theme(legend.position="bottom") + scale_colour_manual(values=colour_palette)
	print(multiplot(plot1,plot4,plot2,plot5,plot3,plot6,cols=3))
	#dev.off()
}
dev.off()

colour_palette <- c("green4","green4","orangered4","salmon3")
# training + prediction data using subseted matrices
for (i in 1:length(top20)) {
	cand <- top20[i]
	pdf(file=paste(cand,".pdf",sep=""),width=11.7,height=8.27)
	temp <- t(rbind(mmatrix_BRCA_PROMOTER_top20_train[cand,ANs],mmatrix_BRCA_BODY_top20_train[cand,ANs],cpm_BRCA_top20_train[cand,ANs],rep("AN",75)))
	temp <- rbind(temp,t(rbind(mmatrix_BRCA_PROMOTER_top20_train[cand,Ts],mmatrix_BRCA_BODY_top20_train[cand,Ts],cpm_BRCA_top20_train[cand,Ts],rep("T",75))))
	temp <- rbind(temp,t(rbind(mmatrix_BRCA_PROMOTER_top20_predict[cand,Ts_toPredict],mmatrix_BRCA_BODY_top20_predict[cand,Ts_toPredict],cpm_BRCA_top20_predict[cand,Ts_toPredict],rep("T_unseen",length(Ts_toPredict)))))
	temp <- rbind(temp,t(rbind(mmatrix_BRCA_PROMOTER_top20_predict[cand,ANs_toPredict],mmatrix_BRCA_BODY_top20_predict[cand,ANs_toPredict],cpm_BRCA_top20_predict[cand,ANs_toPredict],rep("AN_unseen",length(ANs_toPredict)))))
	colnames(temp) <- c("PROMOTER","BODY","CPM","sampleType")
	temp <- as.data.frame(temp)
	temp[,1] <- as.numeric.factor(temp[,1])
	temp[,2] <- as.numeric.factor(temp[,2])
	temp[,3] <- as.numeric.factor(temp[,3])
	
	plot1 <- ggplot(temp[1:150,],aes(x=PROMOTER,y=BODY,colour=sampleType)) +theme_bw() + geom_point(alpha=0.75) + stat_density2d(alpha=0.6) + geom_point(data=temp[151:936,],aes(x=PROMOTER,y=BODY,colour=sampleType),alpha=0.2) + ggtitle(cand) + xlab("Pr. meth.") + ylab("GB. meth.") + theme(legend.position="bottom") + scale_colour_manual(values=colour_palette)
	plot2 <- qplot(temp$sampleType,temp$CPM,geom="boxplot") + theme_bw() + scale_y_log10() +xlab("") +ylab("Expression") + ggtitle(cand)
	plot3 <- ggplot(temp[1:150,],aes(x=PROMOTER,y=CPM,colour=sampleType)) +theme_bw() + scale_y_log10() + geom_point(alpha=0.75) + stat_density2d(alpha=0.6) + geom_point(data=temp[151:936,],aes(x=PROMOTER,y=CPM,colour=sampleType),alpha=0.2) + ggtitle(cand) + xlab("Pr. meth.") + ylab("Expression") + theme(legend.position="bottom") + scale_colour_manual(values=colour_palette)
	plot4 <- ggplot(temp[1:150,],aes(x=BODY,y=CPM,colour=sampleType)) +theme_bw() + scale_y_log10() + geom_point(alpha=0.75) + stat_density2d(alpha=0.6) + geom_point(data=temp[151:936,],aes(x=BODY,y=CPM,colour=sampleType),alpha=0.2) + ggtitle(cand) + xlab("GB meth.") + ylab("Expression") + theme(legend.position="bottom") + scale_colour_manual(values=colour_palette)
	print(multiplot(plot1,plot3,plot2,plot4,cols=2))
	dev.off()
}

colour_palette <- c("red","green4","salmon3","yellowgreen","yellowgreen","salmon3")
names(colour_palette) <- c("FN","FP","TP","TN","AN","T")
# training contours + prediction data points labelled according to prediction
for (i in 1:length(top20)) {
	cand <- top20[i]
	pdf(file=paste(cand,"_predictions.pdf",sep=""),width=11.7,height=8.27)
	TPs <- rownames(predictions[[i]][which(predictions[[i]][Ts_toPredict,2] >= 0.5),])
	FNs <- rownames(predictions[[i]][which(predictions[[i]][Ts_toPredict,2] < 0.5),])
	TNs <- rownames(predictions[[i]][which(predictions[[i]][ANs_toPredict,2] < 0.5),])
	FPs <- rownames(predictions[[i]][which(predictions[[i]][ANs_toPredict,2] >= 0.5),])
	temp <- t(rbind(mmatrix_BRCA_PROMOTER_top20_train[cand,ANs],mmatrix_BRCA_BODY_top20_train[cand,ANs],cpm_BRCA_top20_train[cand,ANs],rep("AN",75)))
	temp <- rbind(temp,t(rbind(mmatrix_BRCA_PROMOTER_top20_train[cand,Ts],mmatrix_BRCA_BODY_top20_train[cand,Ts],cpm_BRCA_top20_train[cand,Ts],rep("T",75))))
	temp <- rbind(temp,t(rbind(mmatrix_BRCA_PROMOTER_top20_predict[cand,c(FPs,FNs)],mmatrix_BRCA_BODY_top20_predict[cand,c(FPs,FNs)],cpm_BRCA_top20_predict[cand,c(FPs,FNs)],c(rep("FP",length(FPs)),rep("FN",length(FNs))))))
	temp <- rbind(temp,t(rbind(mmatrix_BRCA_PROMOTER_top20_predict[cand,c(TNs,TPs)],mmatrix_BRCA_BODY_top20_predict[cand,c(TNs,TPs)],cpm_BRCA_top20_predict[cand,c(TNs,TPs)],c(rep("TN",length(TNs)),rep("TP",length(TPs))))))
	colnames(temp) <- c("PROMOTER","BODY","CPM","sampleType")
	temp <- as.data.frame(temp)
	temp[,1] <- as.numeric.factor(temp[,1])
	temp[,2] <- as.numeric.factor(temp[,2])
	temp[,3] <- as.numeric.factor(temp[,3])

	from_f <- 151
	to_f <- 151+length(c(FNs,FPs))
	from_t <- 151+length(c(FNs,FPs))+1
	to_t <- 936
	plot1 <- ggplot(temp[1:150,],aes(x=PROMOTER,y=CPM,colour=sampleType)) +theme_bw() + scale_y_log10() + stat_density2d(alpha=0.3) + geom_point(data=temp[from_t:to_t,],aes(x=PROMOTER,y=CPM,colour=sampleType),alpha=0.2) + geom_point(data=temp[from_f:to_f,],aes(x=PROMOTER,y=CPM,colour=sampleType),alpha=0.65) + ggtitle(cand) + xlab("Pr. meth.") + ylab("Expression") + theme(legend.position="bottom") + scale_colour_manual(breaks=c("FN","FP","TN","AN","TP","T"),values=colour_palette)
	plot2 <- ggplot(temp[1:150,],aes(x=BODY,y=CPM,colour=sampleType)) +theme_bw() + scale_y_log10() + stat_density2d(alpha=0.3) + geom_point(data=temp[from_t:to_t,],aes(x=BODY,y=CPM,colour=sampleType),alpha=0.2) + geom_point(data=temp[from_f:to_f,],aes(x=BODY,y=CPM,colour=sampleType),alpha=0.65) + ggtitle(cand) + xlab("GB meth.") + ylab("Expression") + theme(legend.position="bottom") + scale_colour_manual(breaks=c("FN","FP","TN","AN","TP","T"),values=colour_palette)
	plot3 <- ggplot(temp[1:150,],aes(x=PROMOTER,y=BODY,colour=sampleType)) +theme_bw() + stat_density2d(alpha=0.3) + geom_point(data=temp[from_t:to_t,],aes(x=PROMOTER,y=BODY,colour=sampleType),alpha=0.2) + geom_point(data=temp[from_f:to_f,],aes(x=PROMOTER,y=BODY,colour=sampleType),alpha=0.65) + ggtitle(cand) + xlab("Pr. meth.") + ylab("GB. meth.") + theme(legend.position="bottom") + scale_colour_manual(breaks=c("FN","FP","TN","AN","TP","T"),values=colour_palette)
	print(multiplot(plot1,plot3,plot2,cols=2))
	dev.off()
}


colour_palette <- c("orange","red")
pdf(file="top10_progression.pdf",width=11.7,height=8.27)

colour_palette <- c("orange","red")
pdf(file="serpins_progression.pdf",width=11.7,height=8.27)

colour_palette <- c("orange","red")
pdf(file="top10_sota_progression.pdf",width=11.7,height=8.27)

# progression data
for (i in 1:length(top10genes_sota)) {
	cand <- top10genes_sota[i]
	temp <- t(rbind(mmatrix_BRCA_PROMOTER[cand,G1],mmatrix_BRCA_BODY[cand,G1],cpm[cand,G1],rep("progresssed",length(G1))))
	temp <- rbind(temp,t(rbind(mmatrix_BRCA_PROMOTER[cand,G2],mmatrix_BRCA_BODY[cand,G2],cpm[cand,G2],rep("nonProgressed",length(G2)))))
	colnames(temp) <- c("PROMOTER","BODY","CPM","sampleType")
	temp <- as.data.frame(temp)
	temp[,1] <- as.numeric.factor(temp[,1])
	temp[,2] <- as.numeric.factor(temp[,2])
	temp[,3] <- as.numeric.factor(temp[,3])
	
	
	plot1 <- qplot(temp$sampleType,temp$CPM,geom="boxplot") + theme_bw() + scale_y_log10() +xlab("") +ylab("Expression") + ggtitle(cand)
	plot2 <- qplot(temp$sampleType,temp$PROMOTER,geom="boxplot") + theme_bw() +xlab("") +ylab("Pr. meth.") + ggtitle(cand)
	plot3 <- qplot(temp$sampleType,temp$BODY,geom="boxplot") + theme_bw()+xlab("") +ylab("GB. meth.") + ggtitle(cand)
	plot4 <- ggplot(temp,aes(x=PROMOTER,y=BODY,colour=sampleType)) +theme_bw() + geom_point(alpha=0.75) + ggtitle(cand) + xlab("Pr. meth.") + ylab("GB. meth.") + stat_density2d(alpha=0.5) + theme(legend.position="bottom") + scale_colour_manual(values=colour_palette)
	plot5 <- ggplot(temp,aes(x=PROMOTER,y=CPM,colour=sampleType)) + theme_bw() + scale_y_log10() + geom_point(alpha=0.75) +ggtitle(cand) + xlab("Pr. meth.") + ylab("Expression") + stat_density2d(alpha=0.5) + theme(legend.position="bottom") + scale_colour_manual(values=colour_palette)
	plot6 <- ggplot(temp,aes(x=BODY,y=CPM,colour=sampleType)) + scale_y_log10() + theme_bw() + geom_point(alpha=0.75) +ggtitle(cand) + xlab("GB. meth.") + ylab("Expression") + stat_density2d(alpha=0.5) + theme(legend.position="bottom") + scale_colour_manual(values=colour_palette)
	print(multiplot(plot1,plot4,plot2,plot5,plot3,plot6,cols=3))
}
dev.off()

# new progression plots
colour_palette <- c("#CC79A7","#009E73")
pdf(file="top10_progression_preSel.pdf",width=11.7,height=8.27)
for (i in 1:10) {
	cand <- as.character(top10_sel[i])
	
	temp <- data.frame(PROMOTER = mmatrix_BRCA_PROMOTER[cand,c(G1,G2)], BODY=mmatrix_BRCA_BODY[cand,c(G1,G2)], CPM=t(cpm[cand,c(G1,G2)]), sampleType=c(rep("progresssed",length(G1)),rep("nonProgressed",length(G2))))
	colnames(temp)[3] <- "CPM"
	
	plot1 <- qplot(temp$sampleType,temp$CPM,geom="boxplot", colour=temp$sampleType) + theme_bw() + theme(panel.grid=element_line(linetype = 0), legend.position="none", plot.title=element_text(size=10)) + scale_y_log10() +xlab("") +ylab("Expression (reads per million)") + ggtitle(paste(cand,"; p-value=",signif(results_progression[cand,2],digits=3),", edgeR",sep="")) + scale_colour_manual(values=colour_palette)
	plot2 <- qplot(temp$sampleType,temp$PROMOTER,geom="boxplot", colour=temp$sampleType) + theme_bw() + theme(panel.grid=element_line(linetype = 0), legend.position="none", plot.title=element_text(size=10)) +xlab("") +ylab("Pr. meth. (M-value)") + ggtitle(paste(cand,"; p-value=",signif(results_progression[cand,3],digits=3),", Welch's t-test",sep="")) + scale_colour_manual(values=colour_palette)
	plot3 <- qplot(temp$sampleType,temp$BODY,geom="boxplot", colour=temp$sampleType) + theme_bw() + theme(panel.grid=element_line(linetype = 0), legend.position="none", plot.title=element_text(size=10)) +xlab("") +ylab("GB. meth. (M-value)") + ggtitle(paste(cand,"; p-value=",signif(results_progression[cand,4],digits=3),", Welch's t-test",sep="")) + scale_colour_manual(values=colour_palette)
	
	correlation_an <- signif(cor(temp[G1,1], temp[G1,2], method="spearman"),digits=2)
	correlation_t <- signif(cor(temp[G2,1], temp[G2,2], method="spearman"),digits=2)
	plot4 <- ggplot(temp,aes(x=PROMOTER,y=BODY,colour=sampleType)) +theme_bw() + theme(panel.grid=element_line(linetype = 0), plot.title=element_text(size=7)) + geom_point(alpha=0.5) + ggtitle(paste(cand,"; cor(progr.)=", correlation_an,", cor(nonProgr.)=", correlation_t,", Spearman's rho",sep="")) + xlab("Pr. meth. (M-value)") + ylab("GB. meth. (M-value)") + stat_density2d(alpha=0.15) + theme(legend.position="bottom") + scale_colour_manual(values=colour_palette)
	
	correlation_an <- signif(cor(temp[G1,1], temp[G1,3], method="spearman"),digits=2)
	correlation_t <- signif(cor(temp[G2,1], temp[G2,3], method="spearman"),digits=2)
	plot5 <- ggplot(temp,aes(x=PROMOTER,y=CPM,colour=sampleType)) + theme_bw() + theme(panel.grid=element_line(linetype = 0), plot.title=element_text(size=7)) + scale_y_log10() + geom_point(alpha=0.5) + ggtitle(paste(cand,"; cor(progr.)=", correlation_an,", cor(nonProgr.)=", correlation_t,", Spearman's rho",sep="")) + xlab("Pr. meth. (M-value)") + ylab("Expression (reads per million)") + stat_density2d(alpha=0.15) + theme(legend.position="bottom") + scale_colour_manual(values=colour_palette)
	
	correlation_an <- signif(cor(temp[G1,2], temp[G1,3], method="spearman"),digits=2)
	correlation_t <- signif(cor(temp[G2,2], temp[G2,3], method="spearman"),digits=2)
	plot6 <- ggplot(temp,aes(x=BODY,y=CPM,colour=sampleType)) + scale_y_log10() + theme_bw() + theme(panel.grid=element_line(linetype = 0), plot.title=element_text(size=7)) + geom_point(alpha=0.5) + ggtitle(paste(cand,"; cor(progr.)=", correlation_an,", cor(nonProgr.)=", correlation_t,", Spearman's rho",sep="")) + xlab("GB. meth. (M-value)") + ylab("Expression (reads per million)") + stat_density2d(alpha=0.15) + theme(legend.position="bottom") + scale_colour_manual(values=colour_palette)
	print(multiplot(plot1,plot4,plot2,plot5,plot3,plot6,cols=3))
}
dev.off()




# updated raw variables' distributions 
colour_palette <- c("green4","green4","orangered4","salmon3")
pdf(file="top10.pdf",width=11.7,height=8.27)
# progression data
for (i in 1:10) {
	cand <- top100[i]
	pdf(file=paste(cand,"_updated.pdf",sep=""),width=11.7,height=8.27)
	temp <- t(rbind(mmatrix_BRCA_PROMOTER[cand,ANs],mmatrix_BRCA_BODY[cand,ANs],cpm[cand,ANs],rep("Adjacent normal samples",length(ANs))))
	temp <- rbind(temp,t(rbind(mmatrix_BRCA_PROMOTER[cand,Ts],mmatrix_BRCA_BODY[cand,Ts],cpm[cand,Ts],rep("Tumours",length(Ts)))))
	colnames(temp) <- c("PROMOTER","BODY","CPM","sampleType")
	temp <- as.data.frame(temp)
	temp[,1] <- as.numeric.factor(temp[,1])
	temp[,2] <- as.numeric.factor(temp[,2])
	temp[,3] <- as.numeric.factor(temp[,3])
	
	plot1 <- qplot(temp$sampleType,temp$CPM,geom="boxplot") + theme_bw() + scale_y_log10() +xlab("") +ylab("Expression") + ggtitle(cand)
	plot2 <- qplot(temp$sampleType,temp$PROMOTER,geom="boxplot") + theme_bw() +xlab("") +ylab("Pr. meth.") + ggtitle(cand)
	plot3 <- qplot(temp$sampleType,temp$BODY,geom="boxplot") + theme_bw()+xlab("") +ylab("GB. meth.") + ggtitle(cand)
	plot4 <- ggplot(temp,aes(x=PROMOTER,y=BODY,colour=sampleType)) +theme_bw() + geom_point(alpha=0.75) + ggtitle(cand) + xlab("Pr. meth.") + ylab("GB. meth.") + stat_density2d(alpha=0.5) + theme(legend.position="bottom") + scale_colour_manual(values=colour_palette)
	plot5 <- ggplot(temp,aes(x=PROMOTER,y=CPM,colour=sampleType)) + theme_bw() + scale_y_log10() + geom_point(alpha=0.75) +ggtitle(cand) + xlab("Pr. meth.") + ylab("Expression") + stat_density2d(alpha=0.5) + theme(legend.position="bottom") + scale_colour_manual(values=colour_palette)
	plot6 <- ggplot(temp,aes(x=BODY,y=CPM,colour=sampleType)) + scale_y_log10() + theme_bw() + geom_point(alpha=0.75) +ggtitle(cand) + xlab("GB. meth.") + ylab("Expression") + stat_density2d(alpha=0.5) + theme(legend.position="bottom") + scale_colour_manual(values=colour_palette)
	print(multiplot(plot1,plot4,plot2,plot5,plot3,plot6,cols=3))
}
dev.off()
