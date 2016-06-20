library(caret)
library(pROC)

# Random Forest
indiv_predictions_pgmBased_cv <- matrix(ncol=71,nrow=10)
colnames(indiv_predictions_pgmBased_cv) <- c(G1,G2)
indiv_predictions_pgmBased_cv_r <- matrix(ncol=71,nrow=10)
colnames(indiv_predictions_pgmBased_cv_r) <- c(G1,G2)
for (fold in 0:13) {
	index_G1_fold <- (fold+1):(fold+1)
	index_G2_fold <- (fold*4+1):((fold+1)*4)
	if (fold == 13) index_G2_fold <- c(index_G2_fold,57)
	current <- top20_genes_perFold_14[1:10,fold+1]
	temp_train <- as.data.frame(t(cpm[current,c(G1[-index_G1_fold],G2[-index_G2_fold])]))
	temp_train <- cbind(temp_train,as.data.frame(t(mmatrix_BRCA_PROMOTER_noncap[current,c(G1[-index_G1_fold],G2[-index_G2_fold])])))
	temp_train <- cbind(temp_train,as.data.frame(t(mmatrix_BRCA_BODY_noncap[current,c(G1[-index_G1_fold],G2[-index_G2_fold])])))
	colnames(temp_train)[1:10] <- paste(colnames(temp_train)[1:10],"_E",sep="")
	colnames(temp_train)[11:20] <- paste(colnames(temp_train)[11:20],"_PR",sep="")
	colnames(temp_train)[21:30] <- paste(colnames(temp_train)[21:30],"_GB",sep="")
	temp_train$class <- factor(c(rep("G1",length(G1[-index_G1_fold])),rep("G2",length(G2[-index_G2_fold]))))
	
	temp_eval <- as.data.frame(t(cpm[current,c(G1[index_G1_fold],G2[index_G2_fold])]))
	temp_eval <- cbind(temp_eval,as.data.frame(t(mmatrix_BRCA_PROMOTER_noncap[current,c(G1[index_G1_fold],G2[index_G2_fold])])))
	temp_eval <- cbind(temp_eval,as.data.frame(t(mmatrix_BRCA_BODY_noncap[current,c(G1[index_G1_fold],G2[index_G2_fold])])))
	colnames(temp_eval)[1:10] <- paste(colnames(temp_eval)[1:10],"_E",sep="")
	colnames(temp_eval)[11:20] <- paste(colnames(temp_eval)[11:20],"_PR",sep="")
	colnames(temp_eval)[21:30] <- paste(colnames(temp_eval)[21:30],"_GB",sep="")
	temp_eval$class <- factor(c(rep("G1",length(G1[index_G1_fold])),rep("G2",length(G2[index_G2_fold]))))
	
	for (i in 0:9) {
		model <- train(temp_train[,((i*3)+1):(3+(i*3))], y=temp_train$class, method="rf")
		indiv_predictions_pgmBased_cv[i+1,c(G1[index_G1_fold],G2[index_G2_fold])] <- predict(model,temp_eval,type="prob")[,1]
		
		model <- train(temp_train[,1:(3+(i*3))], y=temp_train$class, method="rf")
		indiv_predictions_pgmBased_cv_r[i+1,c(G1[index_G1_fold],G2[index_G2_fold])] <- predict(model,temp_eval,type="prob")[,1]
	}
}
indiv_predictions_pgmBased_cv_rf <- indiv_predictions_pgmBased_cv
indiv_predictions_pgmBased_cv_r_rf <- indiv_predictions_pgmBased_cv_r

performances_pgmBased_RF_cv <- matrix(ncol=2,nrow=10)
colnames(performances_pgmBased_RF_cv) <- c("AUC","AUC_r")
for (i in 1:10) {
	performances_pgmBased_RF_cv[i,1] <- auc(predictor=indiv_predictions_pgmBased_cv_rf[i,],response=c(rep("Pos",length(G1)),rep("Neg",length(G2))))
	performances_pgmBased_RF_cv[i,2] <- auc(predictor=indiv_predictions_pgmBased_cv_r_rf[i,],response=c(rep("Pos",length(G1)),rep("Neg",length(G2))))
}


# Bayesian logistic regression
indiv_predictions_pgmBased_cv <- matrix(ncol=71,nrow=10)
colnames(indiv_predictions_pgmBased_cv) <- c(G1,G2)
indiv_predictions_pgmBased_cv_r <- matrix(ncol=71,nrow=10)
colnames(indiv_predictions_pgmBased_cv_r) <- c(G1,G2)
for (fold in 0:13) {
	index_G1_fold <- (fold+1):(fold+1)
	index_G2_fold <- (fold*4+1):((fold+1)*4)
	if (fold == 13) index_G2_fold <- c(index_G2_fold,57)
	current <- top20_genes_perFold_14[1:10,fold+1]
	temp_train <- as.data.frame(t(cpm[current,c(G1[-index_G1_fold],G2[-index_G2_fold])]))
	temp_train <- cbind(temp_train,as.data.frame(t(mmatrix_BRCA_PROMOTER_noncap[current,c(G1[-index_G1_fold],G2[-index_G2_fold])])))
	temp_train <- cbind(temp_train,as.data.frame(t(mmatrix_BRCA_BODY_noncap[current,c(G1[-index_G1_fold],G2[-index_G2_fold])])))
	colnames(temp_train)[1:10] <- paste(colnames(temp_train)[1:10],"_E",sep="")
	colnames(temp_train)[11:20] <- paste(colnames(temp_train)[11:20],"_PR",sep="")
	colnames(temp_train)[21:30] <- paste(colnames(temp_train)[21:30],"_GB",sep="")
	temp_train$class <- factor(c(rep("G1",length(G1[-index_G1_fold])),rep("G2",length(G2[-index_G2_fold]))))
	
	temp_eval <- as.data.frame(t(cpm[current,c(G1[index_G1_fold],G2[index_G2_fold])]))
	temp_eval <- cbind(temp_eval,as.data.frame(t(mmatrix_BRCA_PROMOTER_noncap[current,c(G1[index_G1_fold],G2[index_G2_fold])])))
	temp_eval <- cbind(temp_eval,as.data.frame(t(mmatrix_BRCA_BODY_noncap[current,c(G1[index_G1_fold],G2[index_G2_fold])])))
	colnames(temp_eval)[1:10] <- paste(colnames(temp_eval)[1:10],"_E",sep="")
	colnames(temp_eval)[11:20] <- paste(colnames(temp_eval)[11:20],"_PR",sep="")
	colnames(temp_eval)[21:30] <- paste(colnames(temp_eval)[21:30],"_GB",sep="")
	temp_eval$class <- factor(c(rep("G1",length(G1[index_G1_fold])),rep("G2",length(G2[index_G2_fold]))))
	
	for (i in 0:9) {
		model <- train(temp_train[,((i*3)+1):(3+(i*3))], y=temp_train$class, method="bayesglm")
		indiv_predictions_pgmBased_cv[i+1,c(G1[index_G1_fold],G2[index_G2_fold])] <- predict(model,temp_eval,type="prob")[,1]
		
		model <- train(temp_train[,1:(3+(i*3))], y=temp_train$class, method="bayesglm")
		indiv_predictions_pgmBased_cv_r[i+1,c(G1[index_G1_fold],G2[index_G2_fold])] <- predict(model,temp_eval,type="prob")[,1]
	}
}
indiv_predictions_pgmBased_cv_bayesglm <- indiv_predictions_pgmBased_cv
indiv_predictions_pgmBased_cv_r_bayesglm <- indiv_predictions_pgmBased_cv_r

performances_pgmBased_bayesglm_cv <- matrix(ncol=2,nrow=10)
colnames(performances_pgmBased_bayesglm_cv) <- c("AUC","AUC_r")
for (i in 1:10) {
	performances_pgmBased_bayesglm_cv[i,1] <- auc(predictor=indiv_predictions_pgmBased_cv_bayesglm[i,],response=c(rep("Pos",length(G1)),rep("Neg",length(G2))))
	performances_pgmBased_bayesglm_cv[i,2] <- auc(predictor=indiv_predictions_pgmBased_cv_r_bayesglm[i,],response=c(rep("Pos",length(G1)),rep("Neg",length(G2))))
}


# SVMs with linear kernel
indiv_predictions_pgmBased_cv <- matrix(ncol=71,nrow=10)
colnames(indiv_predictions_pgmBased_cv) <- c(G1,G2)
indiv_predictions_pgmBased_cv_r <- matrix(ncol=71,nrow=10)
colnames(indiv_predictions_pgmBased_cv_r) <- c(G1,G2)
for (fold in 0:13) {
	index_G1_fold <- (fold+1):(fold+1)
	index_G2_fold <- (fold*4+1):((fold+1)*4)
	if (fold == 13) index_G2_fold <- c(index_G2_fold,57)
	current <- top20_genes_perFold_14[1:10,fold+1]
	temp <- as.data.frame(t(cpm[current,c(G1,G2)]))
	temp <- cbind(temp,as.data.frame(t(mmatrix_BRCA_PROMOTER_noncap[current,c(G1,G2)])))
	temp <- cbind(temp,as.data.frame(t(mmatrix_BRCA_BODY_noncap[current,c(G1,G2)])))
	colnames(temp)[1:10] <- paste(colnames(temp)[1:10],"_E",sep="")
	colnames(temp)[11:20] <- paste(colnames(temp)[11:20],"_PR",sep="")
	colnames(temp)[21:30] <- paste(colnames(temp)[21:30],"_GB",sep="")
	preProcValues <- preProcess(temp, method = c("center", "scale"))
	
	temp_train <- as.data.frame(t(cpm[current,c(G1[-index_G1_fold],G2[-index_G2_fold])]))
	temp_train <- cbind(temp_train,as.data.frame(t(mmatrix_BRCA_PROMOTER_noncap[current,c(G1[-index_G1_fold],G2[-index_G2_fold])])))
	temp_train <- cbind(temp_train,as.data.frame(t(mmatrix_BRCA_BODY_noncap[current,c(G1[-index_G1_fold],G2[-index_G2_fold])])))
	colnames(temp_train)[1:10] <- paste(colnames(temp_train)[1:10],"_E",sep="")
	colnames(temp_train)[11:20] <- paste(colnames(temp_train)[11:20],"_PR",sep="")
	colnames(temp_train)[21:30] <- paste(colnames(temp_train)[21:30],"_GB",sep="")
	temp_train <- predict(preProcValues, temp_train)
	temp_train$class <- factor(c(rep("G1",length(G1[-index_G1_fold])),rep("G2",length(G2[-index_G2_fold]))))
	
	temp_eval <- as.data.frame(t(cpm[current,c(G1[index_G1_fold],G2[index_G2_fold])]))
	temp_eval <- cbind(temp_eval,as.data.frame(t(mmatrix_BRCA_PROMOTER_noncap[current,c(G1[index_G1_fold],G2[index_G2_fold])])))
	temp_eval <- cbind(temp_eval,as.data.frame(t(mmatrix_BRCA_BODY_noncap[current,c(G1[index_G1_fold],G2[index_G2_fold])])))
	colnames(temp_eval)[1:10] <- paste(colnames(temp_eval)[1:10],"_E",sep="")
	colnames(temp_eval)[11:20] <- paste(colnames(temp_eval)[11:20],"_PR",sep="")
	colnames(temp_eval)[21:30] <- paste(colnames(temp_eval)[21:30],"_GB",sep="")
	temp_eval <- predict(preProcValues, temp_eval)
	temp_eval$class <- factor(c(rep("G1",length(G1[index_G1_fold])),rep("G2",length(G2[index_G2_fold]))))
	
	for (i in 0:9) {
		model <- train(temp_train[,((i*3)+1):(3+(i*3))], y=temp_train$class, method="svmLinear2", probability = TRUE)
		indiv_predictions_pgmBased_cv[i+1,c(G1[index_G1_fold],G2[index_G2_fold])] <- predict(model,temp_eval[,((i*3)+1):(3+(i*3))],type="prob")[,1]
		
		model <- train(temp_train[,1:(3+(i*3))], y=temp_train$class, method="svmLinear2", probability = TRUE)
		indiv_predictions_pgmBased_cv_r[i+1,c(G1[index_G1_fold],G2[index_G2_fold])] <- predict(model,temp_eval[,1:(3+(i*3))],type="prob")[,1]
	}
}
indiv_predictions_pgmBased_cv_svmLinear <- indiv_predictions_pgmBased_cv
indiv_predictions_pgmBased_cv_r_svmLinear <- indiv_predictions_pgmBased_cv_r

performances_pgmBased_svmLinear_cv <- matrix(ncol=2,nrow=10)
colnames(performances_pgmBased_svmLinear_cv) <- c("AUC","AUC_r")
for (i in 1:10) {
	performances_pgmBased_svmLinear_cv[i,1] <- auc(predictor=indiv_predictions_pgmBased_cv_svmLinear[i,],response=c(rep("Pos",length(G1)),rep("Neg",length(G2))))
	performances_pgmBased_svmLinear_cv[i,2] <- auc(predictor=indiv_predictions_pgmBased_cv_r_svmLinear[i,],response=c(rep("Pos",length(G1)),rep("Neg",length(G2))))
}


# SVMs with polynomial kernel
indiv_predictions_pgmBased_cv <- matrix(ncol=71,nrow=10)
colnames(indiv_predictions_pgmBased_cv) <- c(G1,G2)
indiv_predictions_pgmBased_cv_r <- matrix(ncol=71,nrow=10)
colnames(indiv_predictions_pgmBased_cv_r) <- c(G1,G2)
for (fold in 0:13) {
	index_G1_fold <- (fold+1):(fold+1)
	index_G2_fold <- (fold*4+1):((fold+1)*4)
	if (fold == 13) index_G2_fold <- c(index_G2_fold,57)
	current <- top20_genes_perFold_14[1:10,fold+1]
	temp <- as.data.frame(t(cpm[current,c(G1,G2)]))
	temp <- cbind(temp,as.data.frame(t(mmatrix_BRCA_PROMOTER_noncap[current,c(G1,G2)])))
	temp <- cbind(temp,as.data.frame(t(mmatrix_BRCA_BODY_noncap[current,c(G1,G2)])))
	colnames(temp)[1:10] <- paste(colnames(temp)[1:10],"_E",sep="")
	colnames(temp)[11:20] <- paste(colnames(temp)[11:20],"_PR",sep="")
	colnames(temp)[21:30] <- paste(colnames(temp)[21:30],"_GB",sep="")
	preProcValues <- preProcess(temp, method = c("center", "scale"))
	
	temp_train <- as.data.frame(t(cpm[current,c(G1[-index_G1_fold],G2[-index_G2_fold])]))
	temp_train <- cbind(temp_train,as.data.frame(t(mmatrix_BRCA_PROMOTER_noncap[current,c(G1[-index_G1_fold],G2[-index_G2_fold])])))
	temp_train <- cbind(temp_train,as.data.frame(t(mmatrix_BRCA_BODY_noncap[current,c(G1[-index_G1_fold],G2[-index_G2_fold])])))
	colnames(temp_train)[1:10] <- paste(colnames(temp_train)[1:10],"_E",sep="")
	colnames(temp_train)[11:20] <- paste(colnames(temp_train)[11:20],"_PR",sep="")
	colnames(temp_train)[21:30] <- paste(colnames(temp_train)[21:30],"_GB",sep="")
	temp_train <- predict(preProcValues, temp_train)
	temp_train$class <- factor(c(rep("G1",length(G1[-index_G1_fold])),rep("G2",length(G2[-index_G2_fold]))))
	
	temp_eval <- as.data.frame(t(cpm[current,c(G1[index_G1_fold],G2[index_G2_fold])]))
	temp_eval <- cbind(temp_eval,as.data.frame(t(mmatrix_BRCA_PROMOTER_noncap[current,c(G1[index_G1_fold],G2[index_G2_fold])])))
	temp_eval <- cbind(temp_eval,as.data.frame(t(mmatrix_BRCA_BODY_noncap[current,c(G1[index_G1_fold],G2[index_G2_fold])])))
	colnames(temp_eval)[1:10] <- paste(colnames(temp_eval)[1:10],"_E",sep="")
	colnames(temp_eval)[11:20] <- paste(colnames(temp_eval)[11:20],"_PR",sep="")
	colnames(temp_eval)[21:30] <- paste(colnames(temp_eval)[21:30],"_GB",sep="")
	temp_eval <- predict(preProcValues, temp_eval)
	temp_eval$class <- factor(c(rep("G1",length(G1[index_G1_fold])),rep("G2",length(G2[index_G2_fold]))))
	
	for (i in 0:9) {
		model <- train(temp_train[,((i*3)+1):(3+(i*3))], y=temp_train$class, method="svmPoly", prob.model = TRUE)
		indiv_predictions_pgmBased_cv[i+1,c(G1[index_G1_fold],G2[index_G2_fold])] <- predict(model,temp_eval[,((i*3)+1):(3+(i*3))],type="prob")[,1]
		
		model <- train(temp_train[,1:(3+(i*3))], y=temp_train$class, method="svmPoly", prob.model = TRUE)
		indiv_predictions_pgmBased_cv_r[i+1,c(G1[index_G1_fold],G2[index_G2_fold])] <- predict(model,temp_eval[,1:(3+(i*3))],type="prob")[,1]
	}
}
indiv_predictions_pgmBased_cv_svmPoly <- indiv_predictions_pgmBased_cv
indiv_predictions_pgmBased_cv_r_svmPoly <- indiv_predictions_pgmBased_cv_r

performances_pgmBased_svmPoly_cv <- matrix(ncol=2,nrow=10)
colnames(performances_pgmBased_svmPoly_cv) <- c("AUC","AUC_r")
for (i in 1:10) {
	performances_pgmBased_svmPoly_cv[i,1] <- auc(predictor=indiv_predictions_pgmBased_cv_svmPoly[i,],response=c(rep("Pos",length(G1)),rep("Neg",length(G2))))
	performances_pgmBased_svmPoly_cv[i,2] <- auc(predictor=indiv_predictions_pgmBased_cv_r_svmPoly[i,],response=c(rep("Pos",length(G1)),rep("Neg",length(G2))))
}
