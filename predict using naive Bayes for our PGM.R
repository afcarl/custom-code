ids <- results_all[sort(rank(-results_all$Zs_2way),index.return=TRUE)$ix[1:20],14]

# read predictions from files produced by 2way_predict.R
predictions <- NULL
for (i in 1:length(ids)) {
	predictions[[i]] <- read.table(paste(ids[i],".predicted",sep="",collapse=""),header=TRUE)
	rownames(predictions[[i]]) <- predictions[[i]][,1]
}


library(Brobdingnag)
library(pROC)
performances <- matrix(ncol=6,nrow=20)
colnames(performances) <- c("sensitivity","specificity","AUC","runningNaiveBayes_sen","runningNaiveBayes_spe","AUC")

prod_loglik_an_an <- rep(1,length(ANs_toPredict))
prod_loglik_an_t <- rep(1,length(ANs_toPredict))

prod_loglik_t_an <- rep(1,length(Ts_toPredict))
prod_loglik_t_t <- rep(1,length(Ts_toPredict))


for (i in 1:length(predictions)) {
	performances[i,1] <- length(which(predictions[[i]][Ts_toPredict,2] >= 0.5))/length(Ts)
	performances[i,2] <- length(which(predictions[[i]][ANs_toPredict,2] < 0.5))/length(ANs)
	
	performances[i,3] <- auc(predictor=c(predictions[[i]][Ts_toPredict,2],predictions[[i]][ANs_toPredict,2]),response=c(rep("Pos",length(Ts_toPredict)),rep("Neg",length(ANs_toPredict))))
	
	# naive Bayes combination of results
	prod_loglik_t_t <- prod_loglik_t_t * as.brob(exp(-predictions[[i]][Ts_toPredict,3]))
	prod_loglik_t_an <- prod_loglik_t_an * as.brob(exp(-predictions[[i]][Ts_toPredict,4]))
	
	prod_loglik_an_an <- prod_loglik_an_an * as.brob(exp(-predictions[[i]][ANs_toPredict,4]))
	prod_loglik_an_t <- prod_loglik_an_t * as.brob(exp(-predictions[[i]][ANs_toPredict,3]))
	
	performances[i,5] <- length(which(as.double(prod_loglik_an_t / (prod_loglik_an_t+prod_loglik_an_an)) < 0.5))/length(ANs_toPredict)
	performances[i,4] <- length(which(as.double(prod_loglik_t_t / (prod_loglik_t_t+prod_loglik_t_an)) >= 0.5))/length(Ts_toPredict)
	
	performances[i,6] <- auc(predictor=c(as.double(prod_loglik_t_t / (prod_loglik_t_t+prod_loglik_t_an)),as.double(prod_loglik_an_t / (prod_loglik_an_t+prod_loglik_an_an))),response=c(rep("Pos",length(Ts_toPredict)),rep("Neg",length(ANs_toPredict))))
}
