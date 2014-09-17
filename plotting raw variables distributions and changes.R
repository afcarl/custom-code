findDiff <- function(x,y) {
	out <- 0
	if (x < y) out <- y - x
	if (x > y) out <- -(x - y)
	return(out)
}

require(ggplot2)
require(gridExtra)
cb_palette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00","#CC79A7")



plotter1 <- data.frame(metric=rep(density(bmatrix_BRCA_PROMOTER,from=0,to=1)$x,2),density=c(
density(bmatrix_BRCA_PROMOTER[,ANs2],from=0,to=1)$y,
density(bmatrix_BRCA_PROMOTER[,Ts2],from=0,to=1)$y), type=c(rep("ANs_all",512),rep("Ts_all",512)))
#ggplot(plotter1,aes(x=B_value,y=density)) + geom_line(aes(colour=type),lwd=0.5) + scale_colour_manual(values=cb_palette)

plotter2 <- data.frame(metric=rep(density(bmatrix_BRCA_BODY,from=0,to=1)$x,2),density=c(
density(bmatrix_BRCA_BODY[,ANs2],from=0,to=1)$y,
density(bmatrix_BRCA_BODY[,Ts2],from=0,to=1)$y), type=c(rep("ANs_all",512),rep("Ts_all",512)))
#ggplot(plotter2,aes(x=B_value,y=density)) + geom_line(aes(colour=type),lwd=0.5) + scale_colour_manual(values=cb_palette)

plotter3 <- data.frame(metric=rep(density(cpm_BRCA_plusOne[,ANs2],from=0,to=250)$x,2),density=c(
density(cpm_BRCA_plusOne[,ANs2],from=0,to=250)$y,
density(cpm_BRCA_plusOne[,Ts2],from=0,to=250)$y), type=c(rep("ANs_all",512),rep("Ts_all",512)))
ggplot(plotter3,aes(x=metric,y=density)) + geom_line(aes(colour=type),lwd=0.5) + scale_colour_manual(values=cb_palette) + scale_x_log10()

plotter4 <- rbind(plotter1,plotter2,plotter3)
plotter4$element <- c(rep("Promoter methylation",1024),rep("Gene body methylation",1024),rep("Expression",1024))
plotter4$element <- factor(plotter4$element,levels=c("Expression","Promoter methylation","Gene body methylation"))
ggplot(plotter4,aes(x=metric,y=density)) + geom_line(aes(colour=type),lwd=0.5) +facet_wrap(facets=~element,scales="free") +ggtitle("Global distributions of variables")+ theme(legend.position="bottom") + scale_colour_manual(values=cb_palette)

diff_expression <- Vectorize(findDiff)(cpm_BRCA_plusOne[,ANs2],cpm_BRCA_plusOne[,Ts2])
diff_expression_signif <- Vectorize(findDiff)(cpm_BRCA_plusOne[which(expr_BRCA[,4] <= 0.05),ANs2],cpm_BRCA_plusOne[which(expr_BRCA[,4] <= 0.05),Ts2])

diff_promoter_signif <- Vectorize(findDiff)(bmatrix_BRCA_PROMOTER[,ANs2],bmatrix_BRCA_PROMOTER[,Ts2])
diff_promoter <- Vectorize(findDiff)(bmatrix_BRCA_PROMOTER[which(BRCA_PROMOTERtest[,2] <= 0.05),ANs2],bmatrix_BRCA_PROMOTER[which(BRCA_PROMOTERtest[,2] <= 0.05),Ts2])

diff_body <- Vectorize(findDiff)(bmatrix_BRCA_BODY[,ANs2],bmatrix_BRCA_BODY[,Ts2])
diff_body_signif <- Vectorize(findDiff)(bmatrix_BRCA_BODY[which(BRCA_BODYtest[,2] <= 0.05),ANs2],bmatrix_BRCA_BODY[which(BRCA_BODYtest[,2] <= 0.05),Ts2])




plotter4 <- data.frame(delta_cor=rep(density(null_gb,from=-1,to=1)$x,3),density=c(density(null_spear_gb,from=-1,to=1)$y,density(Vectorize(findDiff)(gb_expr_spear_AN,gb_expr_spear_T),from=-1,to=1)$y,density(Vectorize(findDiff)(gb_expr_spear_AN[top1000_pgm],gb_expr_spear_T[top1000_pgm]),from=-1,to=1)$y),type=c(rep("null",512),rep("allGenes",512),rep("top1000",512)))

plotter5 <- data.frame(delta_cor=rep(density(null_pr,from=-1,to=1)$x,3),density=c(density(null_spear_pr,from=-1,to=1)$y,density(Vectorize(findDiff)(pr_expr_spear_AN,pr_expr_spear_T),from=-1,to=1)$y,density(Vectorize(findDiff)(pr_expr_spear_AN[top1000_pgm],pr_expr_spear_T[top1000_pgm]),from=-1,to=1)$y),type=c(rep("null",512),rep("allGenes",512),rep("top1000",512)))

plotter6 <- rbind(plotter4,plotter5)
plotter6$element <- c(rep("Gene body meth. vs expr.",1536),rep("Promoter meth. vs expr.",1536))

ggplot(plotter6,aes(x=delta_cor,y=density)) + geom_line(aes(colour=type),lwd=0.5) + ggtitle("change of gene-wise correlations of methylation vs expression") + facet_wrap(facets=~element) + theme(legend.position="bottom")  + scale_colour_manual(values=cb_palette)