prod_loglik_G1_pred_G2_pred <- rep(1,length(G1_pred))
prod_loglik_G1_pred_G1_pred <- rep(1,length(G1_pred))
prod_loglik_G2_pred_G2_pred <- rep(1,length(G2_pred))
prod_loglik_G2_pred_G1_pred <- rep(1,length(G2_pred))
which_top10 <- which(workingList_BRCA %in% top10)
which_top10 <- which_top10[c(6,7,1,3,2,9,5,10,4,8)]
performances_PGM_top10_noncv <- matrix(ncol=2,nrow=length(which_top10))
colnames(performances_PGM_top10_noncv) <- c("AUC","AUC_r")
for (j in 1:length(which_top10)) {
	i <- which_top10[j]
	rownames(mlogliks[[i]]) <- mlogliks[[i]][,1]
	
	performances_PGM_top10_noncv[j,1] <- auc(predictor=as.double(as.brob(exp(1))^-mlogliks[[i]][,2] / (as.brob(exp(1))^-mlogliks[[i]][,2] + as.brob(exp(1))^-mlogliks[[i]][,3])),response=c(rep("Pos",length(G1_pred)),rep("Neg",length(G2_pred))))
	
	# naive Bayes combination of results
	prod_loglik_G1_pred_G1_pred <- prod_loglik_G1_pred_G1_pred * as.brob(exp(1))^-mlogliks[[i]][G1_pred,2]
	prod_loglik_G1_pred_G2_pred <- prod_loglik_G1_pred_G2_pred * as.brob(exp(1))^-mlogliks[[i]][G1_pred,3]
	
	prod_loglik_G2_pred_G2_pred <- prod_loglik_G2_pred_G2_pred * as.brob(exp(1))^-mlogliks[[i]][G2_pred,3]
	prod_loglik_G2_pred_G1_pred <- prod_loglik_G2_pred_G1_pred * as.brob(exp(1))^-mlogliks[[i]][G2_pred,2]
	
	performances_PGM_top10_noncv[j,2] <- auc(predictor=c(as.double(prod_loglik_G1_pred_G1_pred / (prod_loglik_G1_pred_G1_pred+prod_loglik_G1_pred_G2_pred)),as.double(prod_loglik_G2_pred_G1_pred / (prod_loglik_G2_pred_G1_pred+prod_loglik_G2_pred_G2_pred))),response=c(rep("Pos",length(G1_pred)),rep("Neg",length(G2_pred))))
}