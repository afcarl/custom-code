require(ggplot2)
require(gridExtra)

which <- which(r %in% genes_underConsideration)
cand <- c("STEAP4")
for (i in 1:length(cand)) {
	temp <- t(rbind(bmatrix_BRCA_PROMOTER2[cand[i],ANs2],bmatrix_BRCA_BODY2[cand[i],ANs2],cpm_BRCA_plusOne[cand[i],ANs2],rep("AN",75)))
	temp <- rbind(temp,t(rbind(bmatrix_BRCA_PROMOTER2[cand[i],Ts2],bmatrix_BRCA_BODY2[cand[i],Ts],cpm_BRCA_plusOne[cand[i],Ts2],rep("T",75))))
	colnames(temp) <- c("PROMOTER","BODY","CPM","sampleType")
	temp <- as.data.frame(temp)
	temp[,1] <- as.numeric(levels(temp[,1])[as.integer(temp[,1])])
	temp[,2] <- as.numeric(levels(temp[,2])[as.integer(temp[,2])])
	temp[,3] <- as.numeric(levels(temp[,3])[as.integer(temp[,3])])
	
	ggplot(temp,aes(x=PROMOTER,y=BODY,colour=sampleType)) + geom_point(alpha=0.75) + ggtitle(cand[i]) + xlab("Pr. meth.")  + ylab("GB. meth.")
	qplot(temp$sampleType,temp$CPM,geom="boxplot") + scale_y_log10() +xlab("") +ylab("EXPRESSION") + ggtitle(cand[i]) + ylab("EXPRESSION")
	ggplot(temp,aes(x=PROMOTER,y=CPM,colour=sampleType)) + geom_point(alpha=0.75) +ggtitle(cand[i]) + xlab("Pr. meth.")  + ylab("EXPRESSION")
	ggplot(temp,aes(x=BODY,y=CPM,colour=sampleType)) + geom_point(alpha=0.75) +ggtitle(cand[i]) + xlab("GB. meth.")  + ylab("EXPRESSION")
	
	
	
}
dev.off()

logliks_temp <- cbind(c(ANs_logliks,Ts_logliks),rep("alt. model",150))
logliks_temp <- rbind(logliks_temp,cbind(all_logliks,rep("null",150)))
logliks_temp <- as.data.frame(logliks_temp)
logliks_temp[,1] <- as.numeric(levels(logliks_temp[,1])[as.integer(logliks_temp[,1])])
colnames(logliks_temp) <- c("loglik","model")