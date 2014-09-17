findDiff <- function(x,y) {
	out <- 0
	if (x < y) out <- y - x
	if (x > y) out <- -(x - y)
	return(out)
}

null_spear_pr <- vector(mode="numeric",length=(17728))
null_spear_gb <- vector(mode="numeric",length=(17728))

null_gb <- vector(mode="numeric",length=(17728))
null_pr <- vector(mode="numeric",length=(17728))


for (i in 1:17728) {
	null_gb_an <- cor.test(bmatrix_BRCA_BODY2[workingList_BRCA4[i],ANs2],sample(cpm_BRCA_plusOne[workingList_BRCA4[i],ANs2]),method="spearman")$estimate
	null_gb[i] <- null_gb_an
	null_gb_t <- cor.test(bmatrix_BRCA_BODY2[workingList_BRCA4[i],Ts2],sample(cpm_BRCA_plusOne[workingList_BRCA4[i],Ts2]),method="spearman")$estimate
	null_spear_gb[i] <- findDiff(null_gb_an,null_gb_t)
	null_pr_an <- cor.test(bmatrix_BRCA_PROMOTER2[workingList_BRCA4[i],ANs2],sample(cpm_BRCA_plusOne[workingList_BRCA4[i],ANs2]),method="spearman")$estimate
	null_pr[i] <- null_pr_an
	null_pr_t <- cor.test(bmatrix_BRCA_PROMOTER2[workingList_BRCA4[i],Ts2],sample(cpm_BRCA_plusOne[workingList_BRCA4[i],Ts2]),method="spearman")$estimate
	null_spear_pr[i] <- findDiff(null_pr_an,null_pr_t)
}
require(ggplot2)
require(gridExtra)
cb_palette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00","#CC79A7")

plotter1 <- data.frame(correlation=rep(density(null_gb,from=-1,to=1)$x,5),density=c(
density(null_gb,from=-1,to=1)$y,
density(gb_expr_spear_AN,from=-1,to=1)$y,
density(gb_expr_spear_AN[top1000_pgm],from=-1,to=1)$y,
density(gb_expr_spear_T,from=-1,to=1)$y,
density(gb_expr_spear_T[top1000_pgm],from=-1,to=1)$y),type=c(rep("null",512),rep("AN_all",512),rep("AN_top1000",512),rep("T_all",512),rep("T_top1000",512)))

plotter2 <- data.frame(correlation=rep(density(null_pr,from=-1,to=1)$x,5),density=c(
density(null_pr,from=-1,to=1)$y,
density(pr_expr_spear_AN,from=-1,to=1)$y,
density(pr_expr_spear_AN[top1000_pgm],from=-1,to=1)$y,
density(pr_expr_spear_T,from=-1,to=1)$y,
density(pr_expr_spear_T[top1000_pgm],from=-1,to=1)$y),type=c(rep("null",512),rep("AN_all",512),rep("AN_top1000",512),rep("T_all",512),rep("T_top1000",512)))

plotter3 <- rbind(plotter1,plotter2)
plotter3$element <- c(rep("Gene body meth. vs expr.",2560),rep("Promoter meth. vs expr.",2560))
ggplot(plotter3,aes(x=correlation,y=density)) + geom_line(aes(colour=type),lwd=0.5) +facet_wrap(facets=~element) +ggtitle("correlations of methylation with expression")+ theme(legend.position="bottom")  + scale_colour_manual(values=cb_palette)

plotter4 <- data.frame(delta_cor=rep(density(null_gb,from=-1,to=1)$x,3),density=c(density(null_spear_gb,from=-1,to=1)$y,density(Vectorize(findDiff)(gb_expr_spear_AN,gb_expr_spear_T),from=-1,to=1)$y,density(Vectorize(findDiff)(gb_expr_spear_AN[top1000_pgm],gb_expr_spear_T[top1000_pgm]),from=-1,to=1)$y),type=c(rep("null",512),rep("allGenes",512),rep("top1000",512)))

plotter5 <- data.frame(delta_cor=rep(density(null_pr,from=-1,to=1)$x,3),density=c(density(null_spear_pr,from=-1,to=1)$y,density(Vectorize(findDiff)(pr_expr_spear_AN,pr_expr_spear_T),from=-1,to=1)$y,density(Vectorize(findDiff)(pr_expr_spear_AN[top1000_pgm],pr_expr_spear_T[top1000_pgm]),from=-1,to=1)$y),type=c(rep("null",512),rep("allGenes",512),rep("top1000",512)))

plotter6 <- rbind(plotter4,plotter5)
plotter6$element <- c(rep("Gene body meth. vs expr.",1536),rep("Promoter meth. vs expr.",1536))

ggplot(plotter6,aes(x=delta_cor,y=density)) + geom_line(aes(colour=type),lwd=0.5) + ggtitle("change of gene-wise correlations of methylation vs expression") + facet_wrap(facets=~element) + theme(legend.position="bottom")  + scale_colour_manual(values=cb_palette)