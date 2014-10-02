require(ggplot2)
require(gridExtra)

for (i in 1:length(top20)) {
	cand <- top20[i]
	pdf(file=paste(cand,".pdf",sep=""),paper="a4r")
	par(mfrow=c(2,2))
	temp <- t(rbind(bmatrix_BRCA_PROMOTER[cand,ANs],bmatrix_BRCA_BODY[cand,ANs],cpm_BRCA_plusOne[cand,ANs],rep("AN",75)))
	temp <- rbind(temp,t(rbind(bmatrix_BRCA_PROMOTER[cand,Ts],bmatrix_BRCA_BODY[cand,Ts],cpm_BRCA_plusOne[cand,Ts],rep("T",75))))
	colnames(temp) <- c("PROMOTER","BODY","CPM","sampleType")
	temp <- as.data.frame(temp)
	temp[,1] <- as.numeric.factor(temp[,1])
	temp[,2] <- as.numeric.factor(temp[,2])
	temp[,3] <- as.numeric.factor(temp[,3])
	
	plot1 <- ggplot(temp,aes(x=PROMOTER,y=BODY,colour=sampleType)) +theme_bw() + geom_point(alpha=0.75) + ggtitle(cand) + xlab("Pr. meth.") + ylab("GB. meth.") + stat_density2d(alpha=0.5) + theme(legend.position="bottom")
	plot2 <- qplot(temp$sampleType,temp$CPM,geom="boxplot") + theme_bw() + scale_y_log10() +xlab("") +ylab("Expression") + ggtitle(cand)
	plot3 <- ggplot(temp,aes(x=PROMOTER,y=CPM,colour=sampleType)) + theme_bw() + geom_point(alpha=0.75) +ggtitle(cand) + xlab("Pr. meth.")  + ylab("Expression") + stat_density2d(alpha=0.5) + theme(legend.position="bottom")
	plot4 <- ggplot(temp,aes(x=BODY,y=CPM,colour=sampleType)) + theme_bw() + geom_point(alpha=0.75) +ggtitle(cand) + xlab("GB. meth.")  + ylab("Expression") + stat_density2d(alpha=0.5) + theme(legend.position="bottom")
	print(multiplot(plot1,plot3,plot2,plot4,cols=2))
	dev.off()
}

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
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                    ncol = cols, nrow = ceiling(numPlots/cols))
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

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}
