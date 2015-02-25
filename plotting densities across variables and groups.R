require(ggplot2)
cb_palette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00","#CC79A7")

# Expression
density_expr_g1 <- density(log10(unlist(cpm[,G1])+1),from=0,to=3.3,n=2048)
density_expr_g2 <- density(log10(unlist(cpm[,G2])+1),from=0,to=3.3,n=2048)
density_expr_g1$y <- density_expr_g1$y/sum(density_expr_g1$y)*102.4
density_expr_g2$y <- density_expr_g2$y/sum(density_expr_g2$y)*102.4

plotter_dens_expr <- data.frame(metric=rep(density_expr_g1$x,2),density=c(density_expr_g1$y,density_expr_g2$y),set=c(rep("ANs",length(density_expr_g1$x)),rep("Ts",length(density_expr_g2$x))))

plot_expr <- ggplot(plotter_dens_expr,aes(x=metric,y=density,colour=set)) + geom_line() + theme_bw() + scale_colour_manual(values=cb_palette) + xlab("log10(RPM+1)")

# Promoter
density_promoter_g1 <- density(unlist(mmatrix_BRCA_PROMOTER[,G1]),from=-7,to=7,n=2048)
density_promoter_g2 <- density(unlist(mmatrix_BRCA_PROMOTER[,G2]),from=-7,to=7,n=2048)
density_promoter_g1$y <- density_promoter_g1$y/sum(density_promoter_g1$y)*102.4
density_promoter_g2$y <- density_promoter_g2$y/sum(density_promoter_g2$y)*102.4

plotter_dens_promoter <- data.frame(metric=rep(density_promoter_g1$x,2),density=c(density_promoter_g1$y,density_promoter_g2$y),set=c(rep("ANs",length(density_promoter_g1$x)),rep("Ts",length(density_promoter_g2$x)))) 

plot_promoter <- ggplot(plotter_dens_promoter,aes(x=metric,y=density,colour=set)) + geom_line() + theme_bw() + scale_colour_manual(values=cb_palette) + ylab("density across genes") + xlab("M-value")
# GB
density_body_g1 <- density(unlist(mmatrix_BRCA_BODY[,G1]),from=-7,to=7,n=2048)
density_body_g2 <- density(unlist(mmatrix_BRCA_BODY[,G2]),from=-7,to=7,n=2048)
density_body_g1$y <- density_body_g1$y/sum(density_body_g1$y)*102.4
density_body_g2$y <- density_body_g2$y/sum(density_body_g2$y)*102.4

plotter_dens_body <- data.frame(metric=rep(density_body_g1$x,2),density=c(density_body_g1$y,density_body_g2$y),set=c(rep("ANs",length(density_body_g1$x)),rep("Ts",length(density_body_g2$x)))) 

plot_body <- ggplot(plotter_dens_body,aes(x=metric,y=density,colour=set)) + geom_line() + theme_bw() + scale_colour_manual(values=cb_palette) + ylab("density across genes") + xlab("M-value")


pdf(file="distributions of variables across genes.pdf",width=11.7,height=8.27)
pdf(file="distributions of variables across genes.pdf")
print(multiplot(plot_expr,plot_promoter,plot_body,cols=3))
dev.off()


plotter <- rbind(plotter_dens_expr,plotter_dens_body,plotter_dens_promoter)
plotter$variable <- c(rep("Gene expression",nrow(plotter_dens_expr)),rep("Gene body methylation",nrow(plotter_dens_body)),rep("Promoter methylation",nrow(plotter_dens_promoter)))
plotter$variable <- factor(plotter$variable,levels=c("Gene expression","Promoter methylation","Gene body methylation"))

pdf(file="temp1.pdf",width=12,height=5)
ggplot(plotter,aes(x=metric,y=density,colour=set)) + geom_line() + theme_bw() + scale_colour_manual(values=cb_palette) + ylab("density across genes") + xlab("M-value") + facet_wrap(facets=~variable) + theme(legend.position="bottom") + xlim(c(0,3.3))
dev.off()

pdf(file="temp2.pdf",width=12,height=5)
ggplot(plotter,aes(x=metric,y=density,colour=set)) + geom_line() + theme_bw() + scale_colour_manual(values=cb_palette) + ylab("density across genes") + xlab("M-value") + facet_wrap(facets=~variable) + theme(legend.position="bottom") + xlim(c(-7,7))
dev.off()