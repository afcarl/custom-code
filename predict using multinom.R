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
performances_binomial <- matrix(ncol=6,nrow=nrow(top100))
colnames(performances_binomial) <- c("sensitivity","specificity","AUC","sensitivity_comb","specificity_comb","AUC_comb")

G1_pred <- G1[-sample_G1]
G2_pred <- G2[-sample_G2]
G1_train <- G1[sample_G1]
G2_train <- G2[sample_G1]

models <- NULL
models_comb <- NULL
dfs_train <- data.frame(class=c(rep("AN",length(G1_train)),rep("T",length(G2_train))))
rownames(dfs_train) <- c(G1_train,G2_train)
dfs_predict <- data.frame(class=c(rep("AN",length(G1_pred)),rep("T",length(G2_pred))))
rownames(dfs_predict) <- c(G1_pred,G2_pred)
for (i in 1:nrow(top100)) {
	cand <- as.character(top100[i,1])
	df_train <- data.frame(promoter_meth = mmatrix_BRCA_PROMOTER[cand,c(G1_train,G2_train)], gb_meth = mmatrix_BRCA_BODY[cand,c(G1_train,G2_train)], expression = t(cpm_plusOne[cand,c(G1_train,G2_train)]), class=c(rep("AN",length(G1_train)),rep("T",length(G2_train))))
	df_predict <- data.frame(promoter_meth = mmatrix_BRCA_PROMOTER[cand,c(G1_pred,G2_pred)], gb_meth = mmatrix_BRCA_BODY[cand,c(G1_pred,G2_pred)], expression = t(cpm_plusOne[cand,c(G1_pred,G2_pred)]), class=c(rep("AN",length(G1_pred)),rep("T",length(G2_pred))))
	colnames(df_predict)[3] <- "expression"
	colnames(df_train)[3] <- "expression"
	
	# independent model
	models[[i]] <- glm(class ~ promoter_meth + gb_meth + expression, data=df_train,family="binomial")
	
	predictions_g2 <- predict(models[[i]],df_predict[G2_pred,],"response")
	predictions_g1 <- predict(models[[i]],df_predict[G1_pred,],"response")
	
	performances_binomial[i,1] <- length(which(predictions_g2 >= 0.5))/length(G2_pred)
	performances_binomial[i,2] <- length(which(predictions_g1 < 0.5))/length(G1_pred)
	performances_binomial[i,3] <- auc(predictor=c(predictions_g2,predictions_g1),response=c(rep("Pos",length(G2_pred)),rep("Neg",length(G1_pred))))
	
	# combined models
	colnames(df_train)[1:3] <- paste(colnames(df_train)[1:3],"_",i,sep="")
	colnames(df_predict)[1:3] <- paste(colnames(df_predict)[1:3],"_",i,sep="")
	dfs_train <- cbind(dfs_train,df_train[,1:3])
	dfs_predict <- cbind(dfs_predict,df_predict[,1:3])
	models_comb[[i]] <- glm(as.formula(paste(c("class ~ ",paste(colnames(dfs_predict[,-1]),collapse=" + ")),collapse="")), data=dfs_train,family="binomial")
	predictions_g2_comb <- predict(models_comb[[i]],dfs_predict[G2_pred,],"response")
	predictions_g1_comb <- predict(models_comb[[i]],dfs_predict[G1_pred,],"response")
	performances_binomial[i,4] <- length(which(predictions_g2_comb >= 0.5))/length(G2_pred)
	performances_binomial[i,5] <- length(which(predictions_g1_comb < 0.5))/length(G1_pred)
	performances_binomial[i,6] <- auc(predictor=c(predictions_g2_comb,predictions_g1_comb),response=c(rep("Pos",length(G2_pred)),rep("Neg",length(G1_pred))))
}


# find and build models using epxression data only

samples <- c(ANs,Ts)
top100_cpm <- matrix(ncol=101,nrow=length(samples))
top100_cpm[,1] <- c(rep("AN",length(ANs)),rep("T",length(Ts)))
top100_cpm[,2:101] <- t(cpm_BRCA_plusOne[workingList_BRCA[top100_e],c(ANs,Ts)])
colnames(top100_cpm) <- c("class",paste("G_",1:100,sep=""))
rownames(top100_cpm) <- c(ANs,Ts)
top100_cpm <- as.data.frame(top100_cpm)
as.numeric.factor <- function(x) {as.numeric(levels(x))[x]}
for (i in 1:100) top100_cpm[,1+i] <- as.numeric.factor(top100_cpm[,1+i])

samples <- c(ANs,Ts)
top100_cpm <- matrix(ncol=101,nrow=length(samples))
top100_cpm[,1] <- c(rep("AN",length(ANs)),rep("T",length(Ts)))
top100_cpm[,2:101] <- t(cpm_BRCA_plusOne[workingList_BRCA[top100_e],c(ANs,Ts)])
colnames(top100_cpm) <- c("class",paste("G_",1:100,sep=""))
rownames(top100_cpm) <- c(ANs,Ts)
top100_cpm <- as.data.frame(top100_cpm)
as.numeric.factor <- function(x) {as.numeric(levels(x))[x]}
for (i in 1:100) top100_cpm[,1+i] <- as.numeric.factor(top100_cpm[,1+i])

samples <- c(ANs_toPredict,Ts_toPredict)
top100_cpm_validation <- matrix(ncol=101,nrow=length(samples))
top100_cpm_validation[,1] <- c(rep("AN",length(ANs_toPredict)),rep("T",length(Ts_toPredict)))
top100_cpm_validation[,2:101] <- t(cpm_BRCA_top100e[workingList_BRCA[top100_e],c(ANs_toPredict,Ts_toPredict)])
colnames(top100_cpm_validation) <- c("class",paste("G_",1:100,sep=""))
rownames(top100_cpm_validation) <- c(ANs_toPredict,Ts_toPredict)
top100_cpm_validation <- as.data.frame(top100_cpm_validation)
for (i in 1:100) top100_cpm_validation[,1+i] <- as.numeric.factor(top100_cpm_validation[,1+i])

# search for model using training data
search_e <- step(glm(class ~ G_1 + G_2 + G_3 + G_4 + G_5 + G_6 + G_7 + G_8 + G_9 + G_10 + G_11 + G_12 + G_13 + G_14 + G_15 + G_16 + G_17 + G_18 + G_19 + G_20 + G_21 + G_22 + G_23 + G_24 + G_25 + G_26 + G_27 + G_28 + G_29 + G_30 + G_31 + G_32 + G_33 + G_34 + G_35 + G_36 + G_37 + G_38 + G_39 + G_40 + G_41 + G_42 + G_43 + G_44 + G_45 + G_46 + G_47 + G_48 + G_49 + G_50 + G_51 + G_52 + G_53 + G_54 + G_55 + G_56 + G_57 + G_58 + G_59 + G_60 + G_61 + G_62 + G_63 + G_64 + G_65 + G_66 + G_67 + G_68 + G_69 + G_70 + G_71 + G_72 + G_73 + G_74 + G_75 + G_76 + G_77 + G_78 + G_79 + G_80 + G_81 + G_82 + G_83 + G_84 + G_85 + G_86 + G_87 + G_88 + G_89 + G_90 + G_91 + G_92 + G_93 + G_94 + G_95 + G_96 + G_97 + G_98 + G_99 + G_100,data=top100_cpm,family="binomial"), ~.,trace=1000)

# build the top10 running models using train data and predict for validation data
maxModel <- G_92 + G_76 + G_31 + G_70 + G_83 + G_25 + G_63 + G_39 + G_15 + G_52
top10 <- c(92,76,31,70,83,25,63,39,15,52)
models_e <- NULL
models_e_s <- NULL
models_e[[1]] <- glm(class ~ G_92,data=top100_cpm,family="binomial")
models_e_s[[1]] <- glm(class ~ G_92,data=top100_cpm,family="binomial")
models_e[[2]] <- glm(class ~ G_92 + G_76,data=top100_cpm,family="binomial")
models_e_s[[2]] <- glm(class ~ G_76,data=top100_cpm,family="binomial")
models_e[[3]] <- glm(class ~ G_92 + G_76 + G_31,data=top100_cpm,family="binomial")
models_e_s[[3]] <- glm(class ~ G_31,data=top100_cpm,family="binomial")
models_e[[4]] <- glm(class ~ G_92 + G_76 + G_31 + G_70,data=top100_cpm,family="binomial")
models_e_s[[4]] <- glm(class ~ G_70,data=top100_cpm,family="binomial")
models_e[[5]] <- glm(class ~ G_92 + G_76 + G_31 + G_70 + G_83,data=top100_cpm,family="binomial")
models_e_s[[5]] <- glm(class ~ G_83,data=top100_cpm,family="binomial")
models_e[[6]] <- glm(class ~ G_92 + G_76 + G_31 + G_70 + G_83 + G_25,data=top100_cpm,family="binomial")
models_e_s[[6]] <- glm(class ~ G_25,data=top100_cpm,family="binomial")
models_e[[7]] <- glm(class ~ G_92 + G_76 + G_31 + G_70 + G_83 + G_25 + G_63,data=top100_cpm,family="binomial")
models_e_s[[7]] <- glm(class ~ G_63,data=top100_cpm,family="binomial")
models_e[[8]] <- glm(class ~ G_92 + G_76 + G_31 + G_70 + G_83 + G_25 + G_63 + G_39,data=top100_cpm,family="binomial")
models_e_s[[8]] <- glm(class ~ G_39,data=top100_cpm,family="binomial")
models_e[[9]] <- glm(class ~ G_92 + G_76 + G_31 + G_70 + G_83 + G_25 + G_63 + G_39 + G_15,data=top100_cpm,family="binomial")
models_e_s[[9]] <- glm(class ~ G_15,data=top100_cpm,family="binomial")
models_e[[10]] <- glm(class ~ G_92 + G_76 + G_31 + G_70 + G_83 + G_25 + G_63 + G_39 + G_15 + G_52,data=top100_cpm,family="binomial")
models_e_s[[10]] <- glm(class ~ G_52,data=top100_cpm,family="binomial")

# predict for validation data
performances_binomial_topE <- matrix(ncol=6,nrow=10)
colnames(performances_binomial_topE) <- c("sensitivity","specificity","AUC","running_sensitivity","running_specificity","running_AUC")
for (i in 1:10) {
	# running models
	predictions_t <- predict(models_e[[i]],top100_cpm_validation[Ts_toPredict,],"response")
	predictions_an <- predict(models_e[[i]],top100_cpm_validation[ANs_toPredict,],"response")
	
	performances_binomial_topE[i,4] <- length(which(predictions_t >= 0.5))/length(Ts_toPredict)
	performances_binomial_topE[i,5] <- length(which(predictions_an < 0.5))/length(ANs_toPredict)
	performances_binomial_topE[i,6] <- auc(predictor=c(predictions_t,predictions_an),response=c(rep("Pos",length(Ts_toPredict)),rep("Neg",length(ANs_toPredict))))
	
	# single models
	predictions_t <- predict(models_e_s[[i]],top100_cpm_validation[Ts_toPredict,],"response")
	predictions_an <- predict(models_e_s[[i]],top100_cpm_validation[ANs_toPredict,],"response")
	
	performances_binomial_topE[i,1] <- length(which(predictions_t >= 0.5))/length(Ts_toPredict)
	performances_binomial_topE[i,2] <- length(which(predictions_an < 0.5))/length(ANs_toPredict)
	performances_binomial_topE[i,3] <- auc(predictor=c(predictions_t,predictions_an),response=c(rep("Pos",length(Ts_toPredict)),rep("Neg",length(ANs_toPredict))))
}
