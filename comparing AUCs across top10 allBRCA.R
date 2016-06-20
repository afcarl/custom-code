require(pROC)
library(Brobdingnag)

pvals_aucTest <- matrix(ncol=2,nrow=10,data=0)
colnames(pvals_aucTest) <- c("single","running")

prod_loglik_G1_pred_G2_pred <- rep(1,length(G1_pred))
prod_loglik_G1_pred_G1_pred <- rep(1,length(G1_pred))
prod_loglik_G2_pred_G2_pred <- rep(1,length(G2_pred))
prod_loglik_G2_pred_G1_pred <- rep(1,length(G2_pred))

for (i in 1:10) {
	
	roc_pgm_s <- roc(predictor=as.double(as.brob(exp(1))^-mlogliks[[i]][,3] / (as.brob(exp(1))^-mlogliks[[i]][,3] + as.brob(exp(1))^-mlogliks[[i]][,4])),response=c(rep("Pos",length(G1_pred)),rep("Neg",length(G2_pred))))
	# naive Bayes combination of results
	prod_loglik_G1_pred_G1_pred <- prod_loglik_G1_pred_G1_pred * as.brob(exp(1))^-mlogliks[[i]][G1_pred,3]
	prod_loglik_G1_pred_G2_pred <- prod_loglik_G1_pred_G2_pred * as.brob(exp(1))^-mlogliks[[i]][G1_pred,4]
	
	prod_loglik_G2_pred_G2_pred <- prod_loglik_G2_pred_G2_pred * as.brob(exp(1))^-mlogliks[[i]][G2_pred,4]
	prod_loglik_G2_pred_G1_pred <- prod_loglik_G2_pred_G1_pred * as.brob(exp(1))^-mlogliks[[i]][G2_pred,3]
	roc_pgm_r <- roc(predictor=c(as.double(prod_loglik_G1_pred_G1_pred / (prod_loglik_G1_pred_G1_pred+prod_loglik_G1_pred_G2_pred)),as.double(prod_loglik_G2_pred_G1_pred / (prod_loglik_G2_pred_G1_pred+prod_loglik_G2_pred_G2_pred))),response=c(rep("Pos",length(G1_pred)),rep("Neg",length(G2_pred))))
	
	roc_lr_s <- roc(predictor=indiv_predictions_pgmBased[i,],response=c(rep("Pos",length(G1_pred)),rep("Neg",length(G2_pred))))
	roc_lr_r <- roc(predictor=indiv_predictions_pgmBased_r[i,],response=c(rep("Pos",length(G1_pred)),rep("Neg",length(G2_pred))))
	
	# test and record p-vals
	pvals_aucTest[i,1] <- roc.test(roc_pgm_s, roc_lr_s, method="delong")$p.value
	pvals_aucTest[i,2] <- roc.test(roc_pgm_r, roc_lr_r, method="delong")$p.value
}