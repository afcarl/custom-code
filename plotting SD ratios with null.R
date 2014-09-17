library(ggplot2)
cb_palette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00","#CC79A7")

# methylation
# actual
mmatrix_capped <- mmatrix
for (i in 1:ncol(mmatrix_capped)) mmatrix_capped[which(mmatrix_capped[,i] < -7),i] <- -7
for (i in 1:ncol(mmatrix_capped)) mmatrix_capped[which(mmatrix_capped[,i] > 7),i] <- 7
deviation_meth <- cbind(apply(mmatrix_capped[,ANs],1,sd),apply(mmatrix_capped[,Ts],1,sd))
deviation_meth <- as.data.frame(deviation_meth)
deviation_meth[,3] <- deviation_meth[,2]/deviation_meth[,1]

# null
index <- sample(1:150,75)
deviation_meth_null <- cbind(apply(mmatrix_capped[,c(ANs,Ts)[index]],1,sd),apply(mmatrix_capped[,c(ANs,Ts)[-index]],1,sd))
deviation_meth_null <- as.data.frame(deviation_meth_null)
deviation_meth_null[,3] <- deviation_meth_null[,2]/deviation_meth_null[,1]

density_actual <- density(deviation_meth[,3],from=0,to=6)
density_null <- density(deviation_meth_null[,3],from=0,to=6)

plotter1 <- data.frame(SD_ratio=rep(density_actual$x,2),density=c(density_actual$y,density_null$y),type=c(rep("BRCA",length(density_actual$x)),rep("null",length(density_null$x))))

ggplot(plotter1,aes(x=SD_ratio,y=density)) + geom_line(aes(colour=type),lwd=0.5) + ggtitle("Ratio of Ts and ANs probe-wise standard deviations across genome based off 450k")+ theme(legend.position="bottom")  + scale_colour_manual(values=cb_palette)

# expression
cpm_BRCA_plusOne <- counts_BRCA_plusOne
for (i in 1:ncol(counts_BRCA_plusOne)) cpm_BRCA_plusOne[,i] <- cpm_BRCA_plusOne[,i]/factors_ls[i]

deviation_expr <- cbind(apply(cpm_BRCA_plusOne[,ANs],1,sd),apply(cpm_BRCA_plusOne[,Ts],1,sd))
deviation_expr <- as.data.frame(deviation_expr)
deviation_expr[,3] <- deviation_expr[,2]/deviation_expr[,1]

deviation_expr_null <- cbind(apply(counts_BRCA_plusOne[,c(ANs,Ts)[index]],1,sd),apply(counts_BRCA_plusOne[,c(ANs,Ts)[-index]],1,sd))
deviation_expr_null <- as.data.frame(deviation_expr_null)
deviation_expr_null[,3] <- deviation_expr_null[,2]/deviation_expr_null[,1]

density_actual <- density(deviation_expr[,3],from=0,to=6)
density_null <- density(deviation_expr_null[,3],from=0,to=6)


plotter2 <- data.frame(SD_ratio=rep(density_actual$x,2),density=c(density_actual$y,density_null$y),type=c(rep("BRCA",length(density_actual$x)),rep("null",length(density_null$x))))

ggplot(plotter2,aes(x=SD_ratio,y=density)) + geom_line(aes(colour=type),lwd=0.5) + ggtitle("Ratio of Ts and ANs gene-wise standard deviations across genome based off RNA-seq")+ theme(legend.position="bottom")  + scale_colour_manual(values=cb_palette)

plotter3 <- rbind(plotter1,plotter2)
plotter3$variable <- c(rep("Methylation",nrow(plotter1)),rep("Expression",nrow(plotter2)))

ggplot(plotter3,aes(x=SD_ratio,y=density)) + geom_line(aes(colour=type),lwd=0.5) + ggtitle("Ratio of Ts and ANs standard deviations across genome")+ theme(legend.position="bottom") + scale_colour_manual(values=cb_palette) + facet_wrap(facets=~variable)


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
