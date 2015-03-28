library(ggplot2)
cb_palette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00","#CC79A7")

# calculating NULLs
index <- sample(1:812,82)
deviation_body_null <- cbind(apply(mmatrix_BRCA_BODY[,c(G1,G2)[index]],1,sd,na.rm=TRUE),apply(mmatrix_BRCA_BODY[,c(G1,G2)[-index]],1,sd,na.rm=TRUE))
deviation_promoter_null <- cbind(apply(mmatrix_BRCA_PROMOTER[,c(G1,G2)[index]],1,sd,na.rm=TRUE),apply(mmatrix_BRCA_PROMOTER[,c(G1,G2)[-index]],1,sd,na.rm=TRUE))
deviation_expr_null <- cbind(apply(cpm[,c(G1,G2)[index]],1,sd,na.rm=TRUE),apply(cpm[,c(G1,G2)[-index]],1,sd,na.rm=TRUE))

for (i in 2:10) {
	index <- sample(1:812,82)
	deviation_body_null <- rbind(deviation_body_null,cbind(apply(mmatrix_BRCA_BODY[,c(G1,G2)[index]],1,sd,na.rm=TRUE),apply(mmatrix_BRCA_BODY[,c(G1,G2)[-index]],1,sd,na.rm=TRUE)))
	deviation_promoter_null <- rbind(deviation_promoter_null,cbind(apply(mmatrix_BRCA_PROMOTER[,c(G1,G2)[index]],1,sd,na.rm=TRUE),apply(mmatrix_BRCA_PROMOTER[,c(G1,G2)[-index]],1,sd,na.rm=TRUE)))
	deviation_expr_null <- rbind(deviation_expr_null,cbind(apply(cpm[,c(G1,G2)[index]],1,sd,na.rm=TRUE),apply(cpm[,c(G1,G2)[-index]],1,sd,na.rm=TRUE)))
}

deviation_body_null <- as.data.frame(deviation_body_null)
deviation_body_null[,3] <- deviation_body_null[,2]/deviation_body_null[,1]
deviation_promoter_null <- as.data.frame(deviation_promoter_null)
deviation_promoter_null[,3] <- deviation_promoter_null[,2]/deviation_promoter_null[,1]
deviation_expr_null <- as.data.frame(deviation_expr_null)
deviation_expr_null[,3] <- deviation_expr_null[,2]/deviation_expr_null[,1]

# methylation split into Pr and GB
deviation_body <- cbind(apply(mmatrix_BRCA_BODY[,G1],1,sd,na.rm=TRUE),apply(mmatrix_BRCA_BODY[,G2],1,sd,na.rm=TRUE))
deviation_body <- as.data.frame(deviation_body)
deviation_body[,3] <- deviation_body[,2]/deviation_body[,1]
deviation_promoter <- cbind(apply(mmatrix_BRCA_PROMOTER[,G1],1,sd,na.rm=TRUE),apply(mmatrix_BRCA_PROMOTER[,G2],1,sd,na.rm=TRUE))
deviation_promoter <- as.data.frame(deviation_promoter)
deviation_promoter[,3] <- deviation_promoter[,2]/deviation_promoter[,1]
deviation_expr <- cbind(apply(cpm[,G1],1,sd,na.rm=TRUE),apply(cpm[,G2],1,sd,na.rm=TRUE))
deviation_expr <- as.data.frame(deviation_expr)
deviation_expr[,3] <- deviation_expr[,2]/deviation_expr[,1]

density_actual_gb <- density(deviation_body[,3],from=0,to=6)
density_null_gb <- density(deviation_body_null[,3],from=0,to=6)
density_actual_pr <- density(deviation_promoter[,3],from=0,to=6)
density_null_pr <- density(deviation_promoter_null[,3],from=0,to=6)
density_actual_expr <- density(deviation_expr[,3],from=0,to=6,na.rm=TRUE)
density_null_expr <- density(deviation_expr_null[,3],from=0,to=6,na.rm=TRUE)

# plotter dfs
plotter_gb <- data.frame(SD_ratio=rep(density_actual_gb$x,2),density=c(density_actual_gb$y,density_null_gb$y),set=c(rep("TCGA Breast Cancer dataset",length(density_actual_gb$x)),rep("expected by random",length(density_null_gb$x))))
plotter_pr <- data.frame(SD_ratio=rep(density_actual_pr$x,2),density=c(density_actual_pr$y,density_null_pr$y),set=c(rep("TCGA Breast Cancer dataset",length(density_actual_pr$x)),rep("expected by random",length(density_null_pr$x))))
plotter_expr <- data.frame(SD_ratio=rep(density_actual_expr$x,2),density=c(density_actual_expr$y,density_null_expr$y),set=c(rep("TCGA Breast Cancer dataset",length(density_actual_expr$x)),rep("expected by random",length(density_null_expr$x))))

plotter <- rbind(plotter_expr,plotter_pr,plotter_gb)
plotter$variable <- c(rep("Gene expression",nrow(plotter_expr)),rep("Gene body methylation",nrow(plotter_gb)),rep("Promoter methylation",nrow(plotter_pr)))
plotter$variable <- factor(plotter$variable,levels=c("Gene expression","Promoter methylation","Gene body methylation"))
plotter$set <- factor(plotter$set,levels=c("expected by random","TCGA Breast Cancer dataset"))

# finding significance thresholds
density_null_gb_sd <- density(deviation_body_null[,3],from=-0,to=6,na.rm=T,n=10240)
density_null_gb_sd$y <- cumsum(density_null_gb_sd$y)/max(cumsum(density_null_gb_sd$y))
density_null_pr_sd <- density(deviation_promoter_null[,3],from=0,to=6,na.rm=T,n=10240)
density_null_pr_sd$y <- cumsum(density_null_pr_sd$y)/max(cumsum(density_null_pr_sd$y))
density_null_expr_sd <- density(deviation_expr_null[,3],from=0,to=6,na.rm=T,n=1024)
density_null_expr_sd$y <- cumsum(density_null_expr_sd$y)/max(cumsum(density_null_expr_sd$y))
signif_gb_sd <- c(density_null_gb_sd$x[which(density_null_gb_sd$y >= 0.025)[1]], density_null_gb_sd$x[which(density_null_gb_sd$y >= 0.975)[1]])
signif_pr_sd <- c(density_null_pr_sd$x[which(density_null_pr_sd$y >= 0.025)[1]], density_null_pr_sd$x[which(density_null_pr_sd$y >= 0.975)[1]])
signif_expr_sd <- c(density_null_expr_sd$x[which(density_null_expr_sd$y >= 0.025)[1]], density_null_expr_sd$x[which(density_null_expr_sd$y >= 0.975)[1]])

# finding significant fractions of empirical distributions
density_actual_gb$y <- cumsum(density_actual_gb$y)/max(cumsum(density_actual_gb$y))
density_actual_pr$y <- cumsum(density_actual_pr$y)/max(cumsum(density_actual_pr$y))
density_actual_expr$y <- cumsum(density_actual_expr$y)/max(cumsum(density_actual_expr$y))

1-density_actual_gb$y[which(density_actual_gb$x >= signif_gb_sd[2])[1]]
1-density_actual_pr$y[which(density_actual_pr$x >= signif_pr_sd[2])[1]]
1-density_actual_expr$y[which(density_actual_expr$x >= signif_expr_sd[2])[1]]

density_actual_gb$y[which(density_actual_gb$x >= signif_gb_sd[1])[1]]
density_actual_pr$y[which(density_actual_pr$x >= signif_pr_sd[1])[1]]
density_actual_expr$y[which(density_actual_expr$x >= signif_expr_sd[1])[1]]


df.signif_sd <- data.frame(thresholds=c(signif_expr_sd,signif_pr_sd,signif_gb_sd),labels=rep("5% significance level",6),variable=c("Gene expression","Gene expression","Promoter methylation","Promoter methylation","Gene body methylation","Gene body methylation"))
palette <- c("#2FAC66","#000000","#2FAC66")
plot1 <- ggplot(plotter,aes(x=SD_ratio,y=density)) + geom_line(aes(colour=set),lwd=0.5) + scale_colour_manual(values=palette,labels=c("5% significance level","expected by random","TCGA Breast Cancer\ndataset")) + facet_wrap(facets=~variable) + theme_bw() + theme(legend.position="right", panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.background = element_rect(fill="white")) + xlab("SD ratio = SD(T's) / SD(AN's)") + ylab("density across genes") + geom_vline(data=df.signif_sd, aes(xintercept = thresholds, colour=labels),linetype=3) + geom_hline(yintercept=0,size=0.15)

# bimodal datas
bimodals <- ("cg21428001","cg03363248","cg10029496","cg18287241","cg04450606","cg21936464","cg25743221","cg06805940","cg05484770","cg12476579","cg22549268")

for (i in 1:length(bimodals)){
	max <- max(mmatrix[bimodals[i],])
	min <- min(mmatrix[bimodals[i],])
	plot(density(mmatrix[bimodals[i],ANs]),col="green",xlab="M-value",xlim=c(min,max),type="l",main=c("ANs and Ts distribution of ",paste(bimodals[i])))
	lines(density(mmatrix[bimodals[i],Ts]),col="red")
}

deviation_expr <- deviation_expr[order(-deviation_expr[,3]),]
for (i in 1:50){
	max <- max(log(cpm_BRCA_plusOne[rownames(deviation_expr[i,]),]))
	min <- min(log(cpm_BRCA_plusOne[rownames(deviation_expr[i,]),]))
	dens_an <- density(log(cpm_BRCA_plusOne[rownames(deviation_expr[i,]),ANs]),from=min,to=max)
	dens_t <- density(log(cpm_BRCA_plusOne[rownames(deviation_expr[i,]),Ts]),from=min,to=max)
	plot(dens_an,col="green",xlab="log(RPM)",ylab="density",type="l",main=c("ANs and Ts distribution of ",paste(rownames(deviation_expr[i,]))," expression"))
	lines(dens_t,col="red")
}
