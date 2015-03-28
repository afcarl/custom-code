plotter_expr <- data.frame(SD_ratio=density_actual_expr$x,density=density_actual_expr$y,set=rep("TCGA Breast Cancer dataset",length(density_actual_expr$x)))
plotter_pr <- data.frame(SD_ratio=density_actual_pr$x,density=density_actual_pr$y,set=rep("TCGA Breast Cancer dataset",length(density_actual_pr$x)))
plotter_gb <- data.frame(SD_ratio=density_actual_gb$x,density=density_actual_gb$y,set=rep("TCGA Breast Cancer dataset",length(density_actual_gb$x)))
plotter <- rbind(plotter_expr,plotter_pr,plotter_gb)
plotter$variable <- c(rep("Gene expression",nrow(plotter_expr)),rep("Gene body methylation",nrow(plotter_gb)),rep("Promoter methylation",nrow(plotter_pr)))
plotter$variable <- factor(plotter$variable,levels=c("Gene expression","Promoter methylation","Gene body methylation"))
plotter$set <- factor(plotter$set,levels=c("expected by random","TCGA Breast Cancer dataset"))

palette <- c("#2FAC66","#2FAC66")
plot1 <- ggplot(plotter,aes(x=SD_ratio,y=density)) + geom_line(aes(colour=set),lwd=0.5) + scale_colour_manual(values=palette,labels=c("5% significance level","TCGA Breast Cancer\ndataset")) + facet_wrap(facets=~variable) + theme_bw() + theme(legend.position="right", panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.background = element_rect(fill="white")) + xlab("SD ratio = SD(T's) / SD(AN's)") + ylab("density across genes") + geom_vline(data=df.signif_sd, aes(xintercept = thresholds, colour=labels),linetype=3) + geom_hline(yintercept=0,size=0.15)


plotter1 <- data.frame(correlation=rep(density(null_gb_an,from=-1,to=1,na.rm=T)$x,2),density=c(
density(gb_expr_spear_AN,from=-1,to=1,na.rm=T)$y,
density(gb_expr_spear_T,from=-1,to=1,na.rm=T)$y),
set=c(rep("AN's",512),rep("T's",512)))

plotter2 <- data.frame(correlation=rep(density(null_pr_an,from=-1,to=1,na.rm=T)$x,2),density=c(
density(pr_expr_spear_AN,from=-1,to=1,na.rm=T)$y,
density(pr_expr_spear_T,from=-1,to=1,na.rm=T)$y),
set=c(rep("AN's",512),rep("T's",512)))

plotter3 <- rbind(plotter1,plotter2)
plotter3$element <- c(rep("Gene body meth. vs expr.",1024),rep("Promoter meth. vs expr.",1024))

plotter4 <- data.frame(delta_cor=density(null_spear_gb,from=-1,to=1,na.rm=T)$x,density=
    density(spear_gb,from=-1,to=1,na.rm=T)$y,
    set=rep("TCGA Breast Cancer dataset",512))

plotter5 <- data.frame(delta_cor=density(null_spear_pr,from=-1,to=1,na.rm=T)$x,density=
    density(spear_pr,from=-1,to=1,na.rm=T)$y,
    set=c(rep("TCGA Breast Cancer dataset",512)))

plotter6 <- rbind(plotter4,plotter5)
plotter6$element <- c(rep("Gene body meth. vs expr.",512),rep("Promoter meth. vs expr.",512))

palette <- c("#E69F00","#E69F00", "#56B4E9", "#56B4E9")
plot2 <- ggplot(plotter3,aes(x=correlation,y=density)) + geom_line(aes(colour=set),lwd=0.5) +facet_wrap(facets=~element) + theme_bw() + theme(legend.position="right", panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.background = element_rect(fill="white")) + scale_colour_manual(values=palette) + ylab("density across genes") + geom_vline(data=df.signif_prgb, aes(xintercept = thresholds, colour=set),linetype=3) + geom_hline(yintercept=0,size=0.15)

palette <- c("#2FAC66","#2FAC66")
plot3 <- ggplot(plotter6,aes(x=delta_cor,y=density)) + geom_line(aes(colour=set),lwd=0.5) + facet_wrap(facets=~element) + theme_bw() + theme(legend.position="right", panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.background = element_rect(fill="white")) + scale_colour_manual(values=palette,labels=c("5% significance level","TCGA Breast Cancer\ndataset")) + ylab("density across genes") + xlab("delta correlation = cor(T’s) - cor(AN’s)") + geom_vline(data=df.signif_spear, aes(xintercept = thresholds, colour=labels),linetype=3) + geom_hline(yintercept=0,size=0.15)

pdf(file="test.pdf",height=11.7,width=8.27)
multiplot(plot1,plot2,plot3)
dev.off()