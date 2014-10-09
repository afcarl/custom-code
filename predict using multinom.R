library(pROC)
cpm_BRCA_top20 <- as.matrix(cpm_BRCA_top20)

# using multinom of nnet
library(nnet)
performances_multinomial <- matrix(ncol=3,nrow=20)
colnames(performances_multinomial) <- c("sensitivity","specificity","AUC")
for (i in 1:length(top20)) {
	df_train <- data.frame(promoter_meth = mmatrix_BRCA_PROMOTER_top20_train[top20[i],c(ANs,Ts)], gb_meth = mmatrix_BRCA_BODY_top20_train[top20[i],c(ANs,Ts)], expression = cpm_BRCA_top20_train[top20[i],c(ANs,Ts)], class=c(rep("AN",length(ANs)),rep("T",length(Ts))))
	df_predict <- data.frame(promoter_meth = mmatrix_BRCA_PROMOTER_top20[top20[i],c(ANs_toPredict,Ts_toPredict)], gb_meth = mmatrix_BRCA_BODY_top20[top20[i],c(ANs_toPredict,Ts_toPredict)], expression = cpm_BRCA_top20[top20[i],c(ANs_toPredict,Ts_toPredict)], class=c(rep("AN",length(ANs_toPredict)),rep("T",length(Ts_toPredict))))

	model <- multinom(class ~ promoter_meth + gb_meth + expression, df_train)
	
	predictions_t <- predict(model,df_predict[Ts_toPredict,],"probs")
	predictions_an <- predict(model,df_predict[ANs_toPredict,],"probs")
	
	performances_multinomial[i,1] <- length(which(predictions_t >= 0.5))/length(Ts_toPredict)
	performances_multinomial[i,2] <- length(which(predictions_an < 0.5))/length(ANs_toPredict)
	performances_multinomial[i,3] <- auc(predictor=c(predictions_t,predictions_an),response=c(rep("Pos",length(Ts_toPredict)),rep("Neg",length(ANs_toPredict))))
}

# using glm's binomial
performances_binomial <- matrix(ncol=3,nrow=20)
colnames(performances_binomial) <- c("sensitivity","specificity","AUC")
models <- NULL
for (i in 1:length(top20)) {
	df_train <- data.frame(promoter_meth = mmatrix_BRCA_PROMOTER_top20_train[top20[i],c(ANs,Ts)], gb_meth = mmatrix_BRCA_BODY_top20_train[top20[i],c(ANs,Ts)], expression = t(cpm_BRCA_top20_train[top20[i],c(ANs,Ts)]), class=c(rep("AN",length(ANs)),rep("T",length(Ts))))
	df_predict <- data.frame(promoter_meth = mmatrix_BRCA_PROMOTER_top20_predict[top20[i],c(ANs_toPredict,Ts_toPredict)], gb_meth = mmatrix_BRCA_BODY_top20_predict[top20[i],c(ANs_toPredict,Ts_toPredict)], expression = t(cpm_BRCA_top20_predict[top20[i],c(ANs_toPredict,Ts_toPredict)]), class=c(rep("AN",length(ANs_toPredict)),rep("T",length(Ts_toPredict))))
	colnames(df_predict)[3] <- "expression"
	colnames(df_train)[3] <- "expression"
	
	# independent model
	models[[i]] <- glm(class ~ promoter_meth + gb_meth + expression, data=df_train,family="binomial")
	# interaction terms model
	#models[[i]] <- glm(class ~ promoter_meth * expression + gb_meth * expression, data=df_train,family="binomial")
	
	predictions_t <- predict(models[[i]],df_predict[Ts_toPredict,],"response")
	predictions_an <- predict(models[[i]],df_predict[ANs_toPredict,],"response")
	
	performances_binomial[i,1] <- length(which(predictions_t >= 0.5))/length(Ts_toPredict)
	performances_binomial[i,2] <- length(which(predictions_an < 0.5))/length(ANs_toPredict)
	performances_binomial[i,3] <- auc(predictor=c(predictions_t,predictions_an),response=c(rep("Pos",length(Ts_toPredict)),rep("Neg",length(ANs_toPredict))))
}
