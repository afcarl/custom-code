library(ggplot2)
cb_palette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00","#CC79A7")

# methylation split into Pr and GB
# GB
bmatrix_BRCA_BODY2_cap <- bmatrix_BRCA_BODY2
bmatrix_BRCA_BODY2_cap[which(bmatrix_BRCA_BODY2_cap < -7)] <- -7
deviation_body <- cbind(apply(bmatrix_BRCA_BODY2_cap[,ANs],1,sd),apply(bmatrix_BRCA_BODY2_cap[,Ts],1,sd))
deviation_body <- as.data.frame(deviation_body)
deviation_body[,3] <- deviation_body[,2]/deviation_body[,1]
# GB null
index <- sample(1:150,75)
deviation_body_null <- cbind(apply(bmatrix_BRCA_BODY2_cap[,c(ANs,Ts)[index]],1,sd),apply(bmatrix_BRCA_BODY2_cap[,c(ANs,Ts)[-index]],1,sd))
deviation_body_null <- as.data.frame(deviation_body_null)
deviation_body_null[,3] <- deviation_body_null[,2]/deviation_body_null[,1]

density_actual <- density(deviation_body[,3],from=0,to=6)
density_null <- density(deviation_body_null[,3],from=0,to=6)

plotter_gb <- data.frame(SD_ratio=rep(density_actual$x,2),density=c(density_actual$y,density_null$y),type=c(rep("BRCA",length(density_actual$x)),rep("null",length(density_null$x))))

# Pr
bmatrix_BRCA_PROMOTER2_cap <- bmatrix_BRCA_PROMOTER2
bmatrix_BRCA_PROMOTER2_cap[which(bmatrix_BRCA_PROMOTER2_cap < -7)] <- -7
deviation_promoter <- cbind(apply(bmatrix_BRCA_PROMOTER2_cap[,ANs],1,sd),apply(bmatrix_BRCA_PROMOTER2_cap[,Ts],1,sd))
deviation_promoter <- as.data.frame(deviation_promoter)
deviation_promoter[,3] <- deviation_promoter[,2]/deviation_promoter[,1]
# Pr null
deviation_promoter_null <- cbind(apply(bmatrix_BRCA_PROMOTER2_cap[,c(ANs,Ts)[index]],1,sd),apply(bmatrix_BRCA_PROMOTER2_cap[,c(ANs,Ts)[-index]],1,sd))
deviation_promoter_null <- as.data.frame(deviation_promoter_null)
deviation_promoter_null[,3] <- deviation_promoter_null[,2]/deviation_promoter_null[,1]

density_actual <- density(deviation_promoter[,3],from=0,to=6)
density_null <- density(deviation_promoter_null[,3],from=0,to=6)

plotter_pr <- data.frame(SD_ratio=rep(density_actual$x,2),density=c(density_actual$y,density_null$y),type=c(rep("BRCA",length(density_actual$x)),rep("null",length(density_null$x))))

plotter_both <- 

ggplot(plotter_pr,aes(x=SD_ratio,y=density)) + geom_line(aes(colour=type),lwd=0.5) + theme(legend.position="bottom") + scale_colour_manual(values=cb_palette)


# expression
cpm_BRCA_plusOne <- counts_BRCA_plusOne[workingList_BRCA,]
for (i in 1:ncol(cpm_BRCA_plusOne)) cpm_BRCA_plusOne[,i] <- cpm_BRCA_plusOne[,i]/factors_ls[i]

deviation_expr <- cbind(apply(cpm_BRCA_plusOne[,ANs],1,sd),apply(cpm_BRCA_plusOne[,Ts],1,sd))
deviation_expr <- as.data.frame(deviation_expr)
deviation_expr[,3] <- deviation_expr[,2]/deviation_expr[,1]

deviation_expr_null <- cbind(apply(cpm_BRCA_plusOne[,c(ANs,Ts)[index]],1,sd),apply(cpm_BRCA_plusOne[,c(ANs,Ts)[-index]],1,sd))
deviation_expr_null <- as.data.frame(deviation_expr_null)
deviation_expr_null[,3] <- deviation_expr_null[,2]/deviation_expr_null[,1]

density_actual <- density(deviation_expr[,3],from=0,to=6)
density_null <- density(deviation_expr_null[,3],from=0,to=6)

plotter_expr <- data.frame(SD_ratio=rep(density_actual$x,2),density=c(density_actual$y,density_null$y),type=c(rep("BRCA",length(density_actual$x)),rep("null",length(density_null$x))))

ggplot(plotter_expr,aes(x=SD_ratio,y=density)) + geom_line(aes(colour=type),lwd=0.5) + theme(legend.position="bottom")  + scale_colour_manual(values=cb_palette)

plotter <- rbind(plotter_expr,plotter_pr,plotter_gb)
plotter$variable <- c(rep("Gene expression",nrow(plotter_expr)),rep("Gene body methylation",nrow(plotter_gb)),rep("Promoter methylation",nrow(plotter_pr)))
plotter$variable <- factor(plotter$variable,levels=c("Gene expression","Promoter methylation","Gene body methylation"))
plotter$type <- factor(plotter$type,levels=c("null","BRCA"))

ggplot(plotter,aes(x=SD_ratio,y=density)) + geom_line(aes(colour=type),lwd=0.5) + ggtitle("Group-wise ratios of standard deviations across genome") + scale_colour_manual(values=cb_palette) + facet_wrap(facets=~variable) + theme_bw() + theme(legend.position="bottom") + xlab("SD ratio")



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
