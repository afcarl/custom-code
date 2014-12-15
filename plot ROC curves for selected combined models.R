require(ggplot2)
require(pROC)
require(Brobdingnag)
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

prod_loglik_g1_g1 <- rep(1,length(G1_pred))
prod_loglik_g1_g2 <- rep(1,length(G1_pred))
prod_loglik_g2_g1 <- rep(1,length(G2_pred))
prod_loglik_g2_g2 <- rep(1,length(G2_pred))
for (i in 1:5) {
	# naive Bayes combination of results
	prod_loglik_g2_g2 <- prod_loglik_g2_g2 * as.brob(exp(-mlogliks[[i]][G2_pred,3]))
	prod_loglik_g2_g1 <- prod_loglik_g2_g1 * as.brob(exp(-mlogliks[[i]][G2_pred,4]))
	
	prod_loglik_g1_g1 <- prod_loglik_g1_g1 * as.brob(exp(-mlogliks[[i]][G1_pred,4]))
	prod_loglik_g1_g2 <- prod_loglik_g1_g2 * as.brob(exp(-mlogliks[[i]][G1_pred,3]))
}
roc_pgm <- roc(cases=as.double(prod_loglik_g2_g2 / (prod_loglik_g2_g2+prod_loglik_g2_g1)), controls=as.double(prod_loglik_g1_g2 / (prod_loglik_g1_g2+prod_loglik_g1_g1)),ci=T)
roc_logistic <- roc(cases=predict(models_comb[[i]],dfs_predict[G2_pred,],"response"), controls=predict(models_comb[[i]],dfs_predict[G1_pred,],"response"),ci=T)
plotter_5c <- data.frame(Sensitivity=c(roc_pgm$sensitivities,roc_logistic$sensitivities),Specificity=c(roc_pgm$specificities,roc_logistic$specificities),Model=c(rep("top5 PGMs combined",length(roc_pgm$specificities)),rep("top5 logistic regression combined",length(roc_logistic$specificities))))
plot1 <- ggplot(plotter_5c,aes(x=Specificity,y=Sensitivity,colour=Model)) + geom_polygon(alpha=0) + scale_x_reverse() + theme_bw() + theme(legend.position="bottom") + geom_abline(intercept=1,slope=1,size=0.75) + scale_colour_brewer(palette="Set1",labels=paste(c("top20 LR-combined: AUC=","top20 PGM-combined: AUC="),c(signif(roc_logistic$auc,4),signif(roc_pgm$auc,4))," (", c(capture.output(roc_logistic$ci),capture.output(roc_pgm$ci)),")", sep = "")) + ggtitle("ROC plot") + theme(plot.title = element_text(face="bold", size=14),  axis.title.x = element_text(face="bold", size=12), axis.title.y = element_text(face="bold", size=12, angle=90), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.justification=c(1,0),  legend.position=c(1,0), legend.title=element_blank(), legend.key = element_blank(), legend.text = element_text(size=7))

for (i in 6:10) {
	# naive Bayes combination of results
	prod_loglik_g2_g2 <- prod_loglik_g2_g2 * as.brob(exp(-mlogliks[[i]][G2_pred,3]))
	prod_loglik_g2_g1 <- prod_loglik_g2_g1 * as.brob(exp(-mlogliks[[i]][G2_pred,4]))
	
	prod_loglik_g1_g1 <- prod_loglik_g1_g1 * as.brob(exp(-mlogliks[[i]][G1_pred,4]))
	prod_loglik_g1_g2 <- prod_loglik_g1_g2 * as.brob(exp(-mlogliks[[i]][G1_pred,3]))
}
roc_pgm <- roc(cases=as.double(prod_loglik_g2_g2 / (prod_loglik_g2_g2+prod_loglik_g2_g1)), controls=as.double(prod_loglik_g1_g2 / (prod_loglik_g1_g2+prod_loglik_g1_g1)),ci=T)
roc_logistic <- roc(cases=predict(models_comb[[i]],dfs_predict[G2_pred,],"response"), controls=predict(models_comb[[i]],dfs_predict[G1_pred,],"response"),ci=T)
plotter_10c <- data.frame(Sensitivity=c(roc_pgm$sensitivities,roc_logistic$sensitivities),Specificity=c(roc_pgm$specificities,roc_logistic$specificities),Model=c(rep("top10 PGMs combined",length(roc_pgm$specificities)),rep("top10 logistic regression combined",length(roc_logistic$specificities))))
plot2 <- ggplot(plotter_10c,aes(x=Specificity,y=Sensitivity,colour=Model)) + geom_polygon(alpha=0) + scale_x_reverse() + theme_bw() + theme(legend.position="bottom") + geom_abline(intercept=1,slope=1,size=0.75) + scale_colour_brewer(palette="Set1",labels=paste(c("top20 LR-combined: AUC=","top20 PGM-combined: AUC="),c(signif(roc_logistic$auc,4),signif(roc_pgm$auc,4))," (", c(capture.output(roc_logistic$ci),capture.output(roc_pgm$ci)),")", sep = "")) + ggtitle("ROC plot") + theme(plot.title = element_text(face="bold", size=14),  axis.title.x = element_text(face="bold", size=12), axis.title.y = element_text(face="bold", size=12, angle=90), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.justification=c(1,0),  legend.position=c(1,0), legend.title=element_blank(), legend.key = element_blank(), legend.text = element_text(size=7))

for (i in 11:15) {
	# naive Bayes combination of results
	prod_loglik_g2_g2 <- prod_loglik_g2_g2 * as.brob(exp(-mlogliks[[i]][G2_pred,3]))
	prod_loglik_g2_g1 <- prod_loglik_g2_g1 * as.brob(exp(-mlogliks[[i]][G2_pred,4]))
	
	prod_loglik_g1_g1 <- prod_loglik_g1_g1 * as.brob(exp(-mlogliks[[i]][G1_pred,4]))
	prod_loglik_g1_g2 <- prod_loglik_g1_g2 * as.brob(exp(-mlogliks[[i]][G1_pred,3]))
}
roc_pgm <- roc(cases=as.double(prod_loglik_g2_g2 / (prod_loglik_g2_g2+prod_loglik_g2_g1)), controls=as.double(prod_loglik_g1_g2 / (prod_loglik_g1_g2+prod_loglik_g1_g1)),ci=T)
roc_logistic <- roc(cases=predict(models_comb[[i]],dfs_predict[G2_pred,],"response"), controls=predict(models_comb[[i]],dfs_predict[G1_pred,],"response"),ci=T)
plotter_15c <- data.frame(Sensitivity=c(roc_pgm$sensitivities,roc_logistic$sensitivities),Specificity=c(roc_pgm$specificities,roc_logistic$specificities),Model=c(rep("top15 PGMs combined",length(roc_pgm$specificities)),rep("top15 logistic regression combined",length(roc_logistic$specificities))))
plot3 <- ggplot(plotter_15c,aes(x=Specificity,y=Sensitivity,colour=Model)) + geom_polygon(alpha=0) + scale_x_reverse() + theme_bw() + theme(legend.position="bottom") + geom_abline(intercept=1,slope=1,size=0.75) + scale_colour_brewer(palette="Set1",labels=paste(c("top20 LR-combined: AUC=","top20 PGM-combined: AUC="),c(signif(roc_logistic$auc,4),signif(roc_pgm$auc,4))," (", c(capture.output(roc_logistic$ci),capture.output(roc_pgm$ci)),")", sep = "")) + ggtitle("ROC plot") + theme(plot.title = element_text(face="bold", size=14),  axis.title.x = element_text(face="bold", size=12), axis.title.y = element_text(face="bold", size=12, angle=90), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.justification=c(1,0),  legend.position=c(1,0), legend.title=element_blank(), legend.key = element_blank(), legend.text = element_text(size=7))

for (i in 16:20) {
	# naive Bayes combination of results
	prod_loglik_g2_g2 <- prod_loglik_g2_g2 * as.brob(exp(-mlogliks[[i]][G2_pred,3]))
	prod_loglik_g2_g1 <- prod_loglik_g2_g1 * as.brob(exp(-mlogliks[[i]][G2_pred,4]))
	
	prod_loglik_g1_g1 <- prod_loglik_g1_g1 * as.brob(exp(-mlogliks[[i]][G1_pred,4]))
	prod_loglik_g1_g2 <- prod_loglik_g1_g2 * as.brob(exp(-mlogliks[[i]][G1_pred,3]))
}
roc_pgm <- roc(cases=as.double(prod_loglik_g2_g2 / (prod_loglik_g2_g2+prod_loglik_g2_g1)), controls=as.double(prod_loglik_g1_g2 / (prod_loglik_g1_g2+prod_loglik_g1_g1)),ci=T)
roc_logistic <- roc(cases=predict(models_comb[[i]],dfs_predict[G2_pred,],"response"), controls=predict(models_comb[[i]],dfs_predict[G1_pred,],"response"),ci=T)
plotter_20c <- data.frame(Sensitivity=c(roc_pgm$sensitivities,roc_logistic$sensitivities),Specificity=c(roc_pgm$specificities,roc_logistic$specificities),Model=c(rep("top20 PGMs combined",length(roc_pgm$specificities)),rep("top20 logistic regression combined",length(roc_logistic$specificities))))
plot4 <- ggplot(plotter_20c,aes(x=Specificity,y=Sensitivity,colour=Model)) + geom_polygon(alpha=0) + scale_x_reverse() + theme_bw() + theme(legend.position="bottom") + geom_abline(intercept=1,slope=1,size=0.75) + scale_colour_brewer(palette="Set1",labels=paste(c("top20 LR-combined: AUC=","top20 PGM-combined: AUC="),c(signif(roc_logistic$auc,4),signif(roc_pgm$auc,4))," (", c(capture.output(roc_logistic$ci),capture.output(roc_pgm$ci)),")", sep = "")) + ggtitle("ROC plot") + theme(plot.title = element_text(face="bold", size=14),  axis.title.x = element_text(face="bold", size=12), axis.title.y = element_text(face="bold", size=12, angle=90), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.justification=c(1,0),  legend.position=c(1,0), legend.title=element_blank(), legend.key = element_blank(), legend.text = element_text(size=7))

pdf(file="4ROCs_comparison.pdf",width=11.7,height=8.27)
print(multiplot(plot1,plot3,plot2,plot4,cols=2))
dev.off()