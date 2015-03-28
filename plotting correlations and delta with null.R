findDiff <- function(x,y) {
	if(is.na(x) | is.na(y)) return(NA)
	out <- 0
	if (x < y) out <- y - x else if (x > y) out <- -(x - y)
	return(out)
}
# actual
gb_expr_spear_AN <- gb_expr_spear_T <- pr_expr_spear_AN <- pr_expr_spear_T <- vector(mode="numeric",length=17728)
for (i in 1:17728) {
	gb_expr_spear_AN[i] <- cor(mmatrix_BRCA_BODY[workingList_BRCA[i],G1],t(cpm[workingList_BRCA[cur],G1])[,1],method="spearman")
	gb_expr_spear_T[i] <- cor(mmatrix_BRCA_BODY[workingList_BRCA[i],G2],t(cpm[workingList_BRCA[cur],G2])[,1],method="spearman")
	pr_expr_spear_AN[i] <- cor(mmatrix_BRCA_PROMOTER[workingList_BRCA[i],G1],t(cpm[workingList_BRCA[cur],G1])[,1],method="spearman")
	pr_expr_spear_T[i] <- cor(mmatrix_BRCA_PROMOTER[workingList_BRCA[i],G2],t(cpm[workingList_BRCA[cur],G2])[,1],method="spearman")
}
spear_pr <- pr_expr_spear_T-pr_expr_spear_AN
spear_gb <- gb_expr_spear_T-gb_expr_spear_AN

# NULL
null_gb_an <- null_gb_t <- null_pr_t <- null_pr_an <- vector(mode="numeric",length=10*17728)
for (outer in 1:17728) {
	for (inner in 0:9) {
		cur <- outer
		index <- inner*17728+outer
		null_gb_an[index] <- cor(mmatrix_BRCA_BODY[workingList_BRCA[cur],G1], sample(t(cpm[workingList_BRCA[cur],G1])[,1]),method="spearman")
		null_gb_t[index] <- cor(mmatrix_BRCA_BODY[workingList_BRCA[cur],G2], sample(t(cpm[workingList_BRCA[cur],G2])[,1]),method="spearman")
		null_pr_an[index] <- cor(mmatrix_BRCA_PROMOTER[workingList_BRCA[cur],G1], sample(t(cpm[workingList_BRCA[cur],G1])[,1]),method="spearman")
		null_pr_t[index] <- cor(mmatrix_BRCA_PROMOTER[workingList_BRCA[cur],G2], sample(t(cpm[workingList_BRCA[cur],G2])[,1]),method="spearman")
	}
}
null_spear_pr <- null_pr_t-null_pr_an
null_spear_gb <- null_gb_t-null_gb_an

# finding significance thresholds
density_null_gb_t <- density(null_gb_t,from=-1,to=1,na.rm=T,n=10240)
density_null_gb_t$y <- cumsum(density_null_gb_t$y)/max(cumsum(density_null_gb_t$y))
density_null_gb_an <- density(null_gb_an,from=-1,to=1,na.rm=T,n=10240)
density_null_gb_an$y <- cumsum(density_null_gb_an$y)/max(cumsum(density_null_gb_an$y))
density_null_pr_t <- density(null_pr_t,from=-1,to=1,na.rm=T,n=10240)
density_null_pr_t$y <- cumsum(density_null_pr_t$y)/max(cumsum(density_null_pr_t$y))
density_null_pr_an <- density(null_pr_an,from=-1,to=1,na.rm=T,n=10240)
density_null_pr_an$y <- cumsum(density_null_pr_an$y)/max(cumsum(density_null_pr_an$y))
density_null_spear_pr <- density(null_spear_pr,from=-1,to=1,na.rm=T,n=10240)
density_null_spear_pr$y <- cumsum(density_null_spear_pr$y)/max(cumsum(density_null_spear_pr$y))
density_null_spear_gb <- density(null_spear_gb,from=-1,to=1,na.rm=T,n=10240)
density_null_spear_gb$y <- cumsum(density_null_spear_gb$y)/max(cumsum(density_null_spear_gb$y))

signif_gb_t <- c(density_null_gb_t$x[which(density_null_gb_t$y >= 0.025)[1]], density_null_gb_t$x[which(density_null_gb_t$y >= 0.975)[1]])
signif_gb_an <- c(density_null_gb_an$x[which(density_null_gb_an$y >= 0.025)[1]], density_null_gb_an$x[which(density_null_gb_an$y >= 0.975)[1]])
signif_pr_t <- c(density_null_pr_t$x[which(density_null_pr_t$y >= 0.025)[1]], density_null_pr_t$x[which(density_null_pr_t$y >= 0.975)[1]])
signif_pr_an <- c(density_null_pr_an$x[which(density_null_pr_an$y >= 0.025)[1]], density_null_pr_an$x[which(density_null_pr_an$y >= 0.975)[1]])
signif_spear_pr <- c(density_null_spear_pr$x[which(density_null_spear_pr$y >= 0.025)[1]],density_null_spear_pr$x[which(density_null_spear_pr$y >= 0.975)[1]])
signif_spear_gb <- c(density_null_spear_gb$x[which(density_null_spear_gb$y >= 0.025)[1]],density_null_spear_gb$x[which(density_null_spear_gb$y >= 0.975)[1]])

# finding significant fractions of empirical distributions
d_pr_an <- density(pr_expr_spear_AN,from=-1,to=1,na.rm=T)
d_pr_t <- density(pr_expr_spear_T,from=-1,to=1,na.rm=T)
d_gb_an <- density(gb_expr_spear_AN,from=-1,to=1,na.rm=T)
d_gb_t <- density(gb_expr_spear_T,from=-1,to=1,na.rm=T)

d_pr_an$y <- cumsum(d_pr_an$y)/max(cumsum(d_pr_an$y))
d_pr_t$y <- cumsum(d_pr_t$y)/max(cumsum(d_pr_t$y))
d_gb_an$y <- cumsum(d_gb_an$y)/max(cumsum(d_gb_an$y))
d_gb_t$y <- cumsum(d_gb_t$y)/max(cumsum(d_gb_t$y))

1-d_pr_an$y[which(d_pr_an$x >= signif_pr_an[2])[1]]
1-d_pr_t$y[which(d_pr_t$x >= signif_pr_t[2])[1]]
1-d_gb_an$y[which(d_gb_an$x >= signif_gb_an[2])[1]]
1-d_gb_t$y[which(d_gb_t$x >= signif_gb_t[2])[1]]

d_pr_an$y[which(d_pr_an$x >= signif_pr_an[1])[1]]
d_pr_t$y[which(d_pr_t$x >= signif_pr_t[1])[1]]
d_gb_an$y[which(d_gb_an$x >= signif_gb_an[1])[1]]
d_gb_t$y[which(d_gb_t$x >= signif_gb_t[1])[1]]

d_spear_gb <- density(spear_gb,from=-1,to=1,na.rm=T)
d_spear_pr <- density(spear_pr,from=-1,to=1,na.rm=T)
d_spear_gb$y <- cumsum(d_spear_gb$y)/max(cumsum(d_spear_gb$y))
d_spear_pr$y <- cumsum(d_spear_pr$y)/max(cumsum(d_spear_pr$y))

1-d_spear_gb$y[which(d_spear_gb$x >= signif_spear_gb[2])[1]]
1-d_spear_pr$y[which(d_spear_pr$x >= signif_spear_pr[2])[1]]
d_spear_gb$y[which(d_spear_gb$x >= signif_spear_gb[1])[1]]
d_spear_pr$y[which(d_spear_pr$x >= signif_spear_pr[1])[1]]

# plotting
require(ggplot2)
require(gridExtra)
cb_palette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00","#CC79A7")

plotter1 <- data.frame(correlation=rep(density(null_gb_an,from=-1,to=1,na.rm=T)$x,4),density=c(
density(null_gb_an,from=-1,to=1,na.rm=T)$y,
density(null_gb_t,from=-1,to=1,na.rm=T)$y,
density(gb_expr_spear_AN,from=-1,to=1,na.rm=T)$y,
density(gb_expr_spear_T,from=-1,to=1,na.rm=T)$y),
set=c(rep("AN's expected by random",512),rep("T's expected by random",512),rep("AN's",512),rep("T's",512)))

plotter2 <- data.frame(correlation=rep(density(null_pr_an,from=-1,to=1,na.rm=T)$x,4),density=c(
density(null_pr_an,from=-1,to=1,na.rm=T)$y,
density(null_pr_t,from=-1,to=1,na.rm=T)$y,
density(pr_expr_spear_AN,from=-1,to=1,na.rm=T)$y,
density(pr_expr_spear_T,from=-1,to=1,na.rm=T)$y),
set=c(rep("AN's expected by random",512),rep("T's expected by random",512),rep("AN's",512),rep("T's",512)))

plotter3 <- rbind(plotter1,plotter2)
plotter3$element <- c(rep("Gene body meth. vs expr.",2048),rep("Promoter meth. vs expr.",2048))
plotter3$set <- factor(plotter3$set,levels=c("AN's expected by random","AN's","AN's 5% significance level","T's expected by random","T's","T's 5% significance level"))

df.signif_prgb <- data.frame(thresholds=c(signif_gb_t,signif_pr_t,signif_gb_an,signif_pr_an),set=c(rep("T's 5% significance level",4),rep("AN's 5% significance level",4)),element=c("Gene body meth. vs expr.","Gene body meth. vs expr.","Promoter meth. vs expr.","Promoter meth. vs expr.","Gene body meth. vs expr.","Gene body meth. vs expr.","Promoter meth. vs expr.","Promoter meth. vs expr."))
df.signif_prgb$set <- factor(df.signif_prgb$set,levels=c("AN's expected by random","AN's","AN's 5% significance level","T's expected by random","T's","T's 5% significance level"))

palette <- c("#E69F00","#E69F00", "#000000", "#56B4E9", "#56B4E9","darkgrey")
plot2 <- ggplot(plotter3,aes(x=correlation,y=density)) + geom_line(aes(colour=set),lwd=0.5) +facet_wrap(facets=~element) + theme_bw() + theme(legend.position="right", panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.background = element_rect(fill="white")) + scale_colour_manual(values=palette) + ylab("density across genes") + geom_vline(data=df.signif_prgb, aes(xintercept = thresholds, colour=set),linetype=3) + geom_hline(yintercept=0,size=0.15)

plotter4 <- data.frame(delta_cor=rep(density(null_spear_gb,from=-1,to=1,na.rm=T)$x,2),density=c(
    density(null_spear_gb,from=-1,to=1,na.rm=T)$y,
    density(spear_gb,from=-1,to=1,na.rm=T)$y),
    set=c(rep("expected by random",512),rep("TCGA Breast Cancer dataset",512)))

plotter5 <- data.frame(delta_cor=rep(density(null_spear_pr,from=-1,to=1,na.rm=T)$x,2),density=c(
    density(null_spear_pr,from=-1,to=1,na.rm=T)$y,
    density(spear_pr,from=-1,to=1,na.rm=T)$y),
    set=c(rep("expected by random",512),rep("TCGA Breast Cancer dataset",512)))

plotter6 <- rbind(plotter4,plotter5)
plotter6$element <- c(rep("Gene body meth. vs expr.",1024),rep("Promoter meth. vs expr.",1024))
plotter6$set <- factor(plotter6$set,levels=c("expected by random","TCGA Breast Cancer dataset"))

df.signif_spear <- data.frame(thresholds=c(signif_spear_gb,signif_spear_pr),labels=c("5% significance level","5% significance level","5% significance level","5% significance level"),element=c("Gene body meth. vs expr.","Gene body meth. vs expr.","Promoter meth. vs expr.","Promoter meth. vs expr."))

palette <- c("#2FAC66","#000000","#2FAC66")
plot3 <- ggplot(plotter6,aes(x=delta_cor,y=density)) + geom_line(aes(colour=set),lwd=0.5) + facet_wrap(facets=~element) + theme_bw() + theme(legend.position="right", panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.background = element_rect(fill="white")) + scale_colour_manual(values=palette,labels=c("5% significance level","expected by random","TCGA Breast Cancer\ndataset")) + ylab("density across genes") + xlab("delta correlation = cor(T’s) - cor(AN’s)") + geom_vline(data=df.signif_spear, aes(xintercept = thresholds, colour=labels),linetype=3) + geom_hline(yintercept=0,size=0.15)

pdf(file="test.pdf",height=11.7,width=8.27)
multiplot(plot1,plot2,plot3)
dev.off()