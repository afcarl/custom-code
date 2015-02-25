# predict using all significant genes in the progression dataset
library(Brobdingnag)
library(pROC)
performances_PGM <- matrix(ncol=6,nrow=length(mlogliks))
colnames(performances_PGM) <- c("sensitivity","specificity","AUC","runningNaiveBayes_sen","runningNaiveBayes_spe","AUC")

prod_loglik_G1_pred_G2_pred <- rep(1,length(G1_pred))
prod_loglik_G1_pred_G1_pred <- rep(1,length(G1_pred))

prod_loglik_G2_pred_G2_pred <- rep(1,length(G2_pred))
prod_loglik_G2_pred_G1_pred <- rep(1,length(G2_pred))

for (i in 1:length(mlogliks)) {
	rownames(mlogliks[[i]]) <- mlogliks[[i]][,1]
	performances_PGM[i,1] <- length(which(as.brob(exp(1))^-mlogliks[[i]][G1_pred,3] / (as.brob(exp(1))^-mlogliks[[i]][G1_pred,3] + as.brob(exp(1))^-mlogliks[[i]][G1_pred,4]) >= 0.5))/length(G1_pred)
	performances_PGM[i,2] <- length(which(as.brob(exp(1))^-mlogliks[[i]][G2_pred,3] / (as.brob(exp(1))^-mlogliks[[i]][G2_pred,3] + as.brob(exp(1))^-mlogliks[[i]][G2_pred,4]) < 0.5))/length(G2_pred)
	
	performances_PGM[i,3] <- auc(predictor=as.double(as.brob(exp(1))^-mlogliks[[i]][,3] / (as.brob(exp(1))^-mlogliks[[i]][,3] + as.brob(exp(1))^-mlogliks[[i]][,4])),response=c(rep("Pos",length(G1_pred)),rep("Neg",length(G2_pred))))
	
	# naive Bayes combination of results
	prod_loglik_G1_pred_G1_pred <- prod_loglik_G1_pred_G1_pred * as.brob(exp(1))^-mlogliks[[i]][G1_pred,3]
	prod_loglik_G1_pred_G2_pred <- prod_loglik_G1_pred_G2_pred * as.brob(exp(1))^-mlogliks[[i]][G1_pred,4]
	
	prod_loglik_G2_pred_G2_pred <- prod_loglik_G2_pred_G2_pred * as.brob(exp(1))^-mlogliks[[i]][G2_pred,4]
	prod_loglik_G2_pred_G1_pred <- prod_loglik_G2_pred_G1_pred * as.brob(exp(1))^-mlogliks[[i]][G2_pred,3]
	
	performances_PGM[i,4] <- length(which(as.double(prod_loglik_G1_pred_G1_pred / (prod_loglik_G1_pred_G1_pred+prod_loglik_G1_pred_G2_pred)) >= 0.5))/length(G1_pred)
	performances_PGM[i,5] <- length(which(as.double(prod_loglik_G2_pred_G1_pred / (prod_loglik_G2_pred_G1_pred+prod_loglik_G2_pred_G2_pred)) < 0.5))/length(G2_pred)
	
	performances_PGM[i,6] <- auc(predictor=c(as.double(prod_loglik_G1_pred_G1_pred / (prod_loglik_G1_pred_G1_pred+prod_loglik_G1_pred_G2_pred)),as.double(prod_loglik_G2_pred_G1_pred / (prod_loglik_G2_pred_G1_pred+prod_loglik_G2_pred_G2_pred))),response=c(rep("Pos",length(G1_pred)),rep("Neg",length(G2_pred))))
}




# build logistic regression models using PGM-found genes at each fold

temp_train <- as.data.frame(t(cpm[top100,c(G1_train,G2_train)]))
temp_train <- cbind(temp_train,as.data.frame(t(mmatrix_BRCA_PROMOTER[top100,c(G1_train,G2_train)])))
temp_train <- cbind(temp_train,as.data.frame(t(mmatrix_BRCA_BODY[top100,c(G1_train,G2_train)])))
colnames(temp_train)[1:100] <- paste(colnames(temp_train)[1:100],"_E",sep="")
colnames(temp_train)[101:200] <- paste(colnames(temp_train)[101:200],"_PR",sep="")
colnames(temp_train)[201:300] <- paste(colnames(temp_train)[201:300],"_GB",sep="")
temp_train$class <- factor(c(rep("G1",length(G1_train)),rep("G2",length(G2_train))))

temp_eval <- as.data.frame(t(cpm[top100,c(G1_pred,G2_pred)]))
temp_eval <- cbind(temp_eval,as.data.frame(t(mmatrix_BRCA_PROMOTER[top100,c(G1_pred,G2_pred)])))
temp_eval <- cbind(temp_eval,as.data.frame(t(mmatrix_BRCA_BODY[top100,c(G1_pred,G2_pred)])))
colnames(temp_eval)[1:100] <- paste(colnames(temp_eval)[1:100],"_E",sep="")
colnames(temp_eval)[101:200] <- paste(colnames(temp_eval)[101:200],"_PR",sep="")
colnames(temp_eval)[201:300] <- paste(colnames(temp_eval)[201:300],"_GB",sep="")
temp_eval$class <- factor(c(rep("G1",length(G1_pred)),rep("G2",length(G2_pred))))

models <- models_r <- NULL
indiv_predictions_pgmBased <- matrix(ncol=length(c(G1_pred,G2_pred)), nrow=100)
colnames(indiv_predictions_pgmBased) <- c(G1_pred,G2_pred)
indiv_predictions_pgmBased_r <- matrix(ncol=length(c(G1_pred,G2_pred)), nrow=100)
colnames(indiv_predictions_pgmBased_r) <- c(G1_pred,G2_pred)

for (i in 1:100) {
	models[[i]] <- glm(as.formula(paste("class ~ ",paste(paste(top100[i],"_E",sep=""),paste(top100[i],"_PR",sep=""),paste(top100[i],"_GB",sep=""),sep=" + "))), data=temp_train, family="binomial")
	indiv_predictions_pgmBased[i,c(G1_pred,G2_pred)] <- predict(models[[i]],temp_eval,type="response")

	models_r[[i]] <- glm(as.formula(paste("class ~ ",paste(paste(top100[1:i],"_E",sep="", collapse=" + "), paste(top100[1:i],"_PR",sep="", collapse=" + "), paste(top100[1:i],"_GB",sep="", collapse=" + "),sep=" + "))), data=temp_train, family="binomial")
	indiv_predictions_pgmBased_r[i,c(G1_pred,G2_pred)] <- predict(models_r[[i]],temp_eval,type="response")
}

performances_pgmBased_LR <- matrix(ncol=6,nrow=length(top100))
colnames(performances_pgmBased_LR) <- c("sensitivity","specificity","AUC","sensitivity_r","specificity_r","AUC_r")

for (i in 1:length(top100)) {
	performances_pgmBased_LR[i,1] <- length(which(indiv_predictions_pgmBased[i,1:length(G1_pred)] < 0.5))/length(G1_pred)
	performances_pgmBased_LR[i,2] <- length(which(indiv_predictions_pgmBased[i,(1+length(G1_pred)):(length(G1_pred)+length(G2_pred))] >= 0.5))/length(G2_pred)
	performances_pgmBased_LR[i,3] <- auc(predictor=indiv_predictions_pgmBased[i,],response=c(rep("Pos",length(G1_pred)),rep("Neg",length(G2_pred))))
	performances_pgmBased_LR[i,4] <- length(which(indiv_predictions_pgmBased_r[i,1:length(G1_pred)] < 0.5))/length(G1_pred)
	performances_pgmBased_LR[i,5] <- length(which(indiv_predictions_pgmBased_r[i,(1+length(G1_pred)):(length(G1_pred)+length(G2_pred))] >= 0.5))/length(G2_pred)
	performances_pgmBased_LR[i,6] <- auc(predictor=indiv_predictions_pgmBased_r[i,],response=c(rep("Pos",length(G1_pred)),rep("Neg",length(G2_pred))))
}



# build models using top100 according to SotA methods
top100_sota <- rownames(results_all[order(results_all$sota_fish),])[1:100]

temp_train <- as.data.frame(t(cpm[top100_sota,c(G1_train,G2_train)]))
temp_train <- cbind(temp_train,as.data.frame(t(mmatrix_BRCA_PROMOTER[top100_sota,c(G1_train,G2_train)])))
temp_train <- cbind(temp_train,as.data.frame(t(mmatrix_BRCA_BODY[top100_sota,c(G1_train,G2_train)])))
colnames(temp_train)[1:100] <- paste(colnames(temp_train)[1:100],"_E",sep="")
colnames(temp_train)[101:200] <- paste(colnames(temp_train)[101:200],"_PR",sep="")
colnames(temp_train)[201:300] <- paste(colnames(temp_train)[201:300],"_GB",sep="")
temp_train$class <- factor(c(rep("G1",length(G1_train)),rep("G2",length(G2_train))))

temp_eval <- as.data.frame(t(cpm[top100_sota,c(G1_pred,G2_pred)]))
temp_eval <- cbind(temp_eval,as.data.frame(t(mmatrix_BRCA_PROMOTER[top100_sota,c(G1_pred,G2_pred)])))
temp_eval <- cbind(temp_eval,as.data.frame(t(mmatrix_BRCA_BODY[top100_sota,c(G1_pred,G2_pred)])))
colnames(temp_eval)[1:100] <- paste(colnames(temp_eval)[1:100],"_E",sep="")
colnames(temp_eval)[101:200] <- paste(colnames(temp_eval)[101:200],"_PR",sep="")
colnames(temp_eval)[201:300] <- paste(colnames(temp_eval)[201:300],"_GB",sep="")
temp_eval$class <- factor(c(rep("G1",length(G1_pred)),rep("G2",length(G2_pred))))

models <- models_r <- NULL
indiv_predictions_sotaBased <- matrix(ncol=length(c(G1_pred,G2_pred)), nrow=100)
colnames(indiv_predictions_sotaBased) <- c(G1_pred,G2_pred)
indiv_predictions_sotaBased_r <- matrix(ncol=length(c(G1_pred,G2_pred)), nrow=100)
colnames(indiv_predictions_sotaBased_r) <- c(G1_pred,G2_pred)

for (i in 1:100) {
	models[[i]] <- glm(as.formula(paste("class ~ ",paste(paste(top100_sota[i],"_E",sep=""),paste(top100_sota[i],"_PR",sep=""),paste(top100_sota[i],"_GB",sep=""),sep=" + "))), data=temp_train, family="binomial")
	indiv_predictions_sotaBased[i,c(G1_pred,G2_pred)] <- predict(models[[i]],temp_eval,type="response")

	models_r[[i]] <- glm(as.formula(paste("class ~ ",paste(paste(top100_sota[1:i],"_E",sep="", collapse=" + "), paste(top100_sota[1:i],"_PR",sep="", collapse=" + "), paste(top100_sota[1:i],"_GB",sep="", collapse=" + "),sep=" + "))), data=temp_train, family="binomial")
	indiv_predictions_sotaBased_r[i,c(G1_pred,G2_pred)] <- predict(models_r[[i]],temp_eval,type="response")
}

performances_sotaBased_LR <- matrix(ncol=6,nrow=length(top100_sota))
colnames(performances_sotaBased_LR) <- c("sensitivity","specificity","AUC","sensitivity_r","specificity_r","AUC_r")

for (i in 1:length(top100_sota)) {
	performances_sotaBased_LR[i,1] <- length(which(indiv_predictions_sotaBased[i,1:length(G1_pred)] < 0.5))/length(G1_pred)
	performances_sotaBased_LR[i,2] <- length(which(indiv_predictions_sotaBased[i,(1+length(G1_pred)):(length(G1_pred)+length(G2_pred))] >= 0.5))/length(G2_pred)
	performances_sotaBased_LR[i,3] <- auc(predictor=indiv_predictions_sotaBased[i,],response=c(rep("Pos",length(G1_pred)),rep("Neg",length(G2_pred))))
	performances_sotaBased_LR[i,4] <- length(which(indiv_predictions_sotaBased_r[i,1:length(G1_pred)] < 0.5))/length(G1_pred)
	performances_sotaBased_LR[i,5] <- length(which(indiv_predictions_sotaBased_r[i,(1+length(G1_pred)):(length(G1_pred)+length(G2_pred))] >= 0.5))/length(G2_pred)
	performances_sotaBased_LR[i,6] <- auc(predictor=indiv_predictions_sotaBased_r[i,],response=c(rep("Pos",length(G1_pred)),rep("Neg",length(G2_pred))))
}