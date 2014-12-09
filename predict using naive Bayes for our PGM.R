#ids <- results_all[sort(rank(-results_all$Zs_2way),index.return=TRUE)$ix[1:20],10]
ids <- c(1657,9059,9227,11795,16267,9007,16764,1596,8298,7642,6508,12925,9878,11222,14399,14831,10755,12523,12127,9791)

# read predictions from files produced by 2way_predict.R
predictions <- NULL
for (i in 1:length(ids)) {
	predictions[[i]] <- read.table(paste(ids[i],".predicted",sep="",collapse=""),header=TRUE)
	rownames(predictions[[i]]) <- predictions[[i]][,1]
}
save(predictions,file="predictions_BRCA.RData")

library(Brobdingnag)
library(pROC)
G1_pred <- G1[56:82]
G2_pred <- G2[487:730]

performances_PGM <- matrix(ncol=6,nrow=100)
colnames(performances_PGM) <- c("sensitivity","specificity","AUC","runningNaiveBayes_sen","runningNaiveBayes_spe","AUC")
prod_loglik_g1_g1 <- rep(1,length(G1_pred))
prod_loglik_g1_g2 <- rep(1,length(G1_pred))
prod_loglik_g2_g1 <- rep(1,length(G2_pred))
prod_loglik_g2_g2 <- rep(1,length(G2_pred))

for (i in 1:length(mlogliks)) {
	rownames(mlogliks[[i]]) <- mlogliks[[i]][,1]
	performances_PGM[i,2] <- length(which(as.brob(exp(1))^-mlogliks[[i]][G1_pred,3] / (as.brob(exp(1))^-mlogliks[[i]][G1_pred,3] + as.brob(exp(1))^-mlogliks[[i]][G1_pred,4]) < 0.5))/length(G1_pred)
	performances_PGM[i,1] <- length(which(as.brob(exp(1))^-mlogliks[[i]][G2_pred,3] / (as.brob(exp(1))^-mlogliks[[i]][G2_pred,3] + as.brob(exp(1))^-mlogliks[[i]][G2_pred,4]) >= 0.5))/length(G2_pred)
	
	performances_PGM[i,3] <- auc(predictor=as.double(as.brob(exp(1))^-mlogliks[[i]][,3] / (as.brob(exp(1))^-mlogliks[[i]][,3] + as.brob(exp(1))^-mlogliks[[i]][,4])),response=c(rep("Pos",length(G1_pred)),rep("Neg",length(G2_pred))))
	
	# naive Bayes combination of results
	prod_loglik_g2_g2 <- prod_loglik_g2_g2 * as.brob(exp(-mlogliks[[i]][G2_pred,3]))
	prod_loglik_g2_g1 <- prod_loglik_g2_g1 * as.brob(exp(-mlogliks[[i]][G2_pred,4]))
	
	prod_loglik_g1_g1 <- prod_loglik_g1_g1 * as.brob(exp(-mlogliks[[i]][G1_pred,4]))
	prod_loglik_g1_g2 <- prod_loglik_g1_g2 * as.brob(exp(-mlogliks[[i]][G1_pred,3]))
	
	performances_PGM[i,5] <- length(which(as.double(prod_loglik_g1_g2 / (prod_loglik_g1_g2+prod_loglik_g1_g1)) < 0.5))/length(G1_pred)
	performances_PGM[i,4] <- length(which(as.double(prod_loglik_g2_g2 / (prod_loglik_g2_g2+prod_loglik_g2_g1)) >= 0.5))/length(G2_pred)
	
	performances_PGM[i,6] <- auc(predictor=c(as.double(prod_loglik_g2_g2 / (prod_loglik_g2_g2+prod_loglik_g2_g1)),as.double(prod_loglik_g1_g2 / (prod_loglik_g1_g2+prod_loglik_g1_g1))),response=c(rep("Pos",length(G2_pred)),rep("Neg",length(G1_pred))))
}


# selection of most predictive genes among the top100 using the GLM's logistic regression
top100_posteriors_train <- matrix(ncol=101,nrow=length(c(ANs,Ts)))
for (i in 1:length(predictions)) {
	top100_posteriors_train[,1+i] <- predictions[[i]][,2]
}
top100_posteriors_train[,1] <- c(rep("AN",length(ANs)),rep("T",length(Ts)))
rownames(top100_posteriors_train) <- c(ANs,Ts)
colnames(top100_posteriors_train) <- c("class",paste("G_",1:100,sep=""))
top100_posteriors_train <- as.data.frame(top100_posteriors_train)
as.numeric.factor <- function(x) {as.numeric(levels(x))[x]}
for (i in 1:100) top100_posteriors_train[,1+i] <- as.numeric.factor(top100_posteriors_train[,1+i])

search <- step(glm(class ~ G_1 + G_2 + G_3 + G_4 + G_5 + G_6 + G_7 + G_8 + G_9 + G_10 + G_11 + G_12 + G_13 + G_14 + G_15 + G_16 + G_17 + G_18 + G_19 + G_20 + G_21 + G_22 + G_23 + G_24 + G_25 + G_26 + G_27 + G_28 + G_29 + G_30 + G_31 + G_32 + G_33 + G_34 + G_35 + G_36 + G_37 + G_38 + G_39 + G_40 + G_41 + G_42 + G_43 + G_44 + G_45 + G_46 + G_47 + G_48 + G_49 + G_50 + G_51 + G_52 + G_53 + G_54 + G_55 + G_56 + G_57 + G_58 + G_59 + G_60 + G_61 + G_62 + G_63 + G_64 + G_65 + G_66 + G_67 + G_68 + G_69 + G_70 + G_71 + G_72 + G_73 + G_74 + G_75 + G_76 + G_77 + G_78 + G_79 + G_80 + G_81 + G_82 + G_83 + G_84 + G_85 + G_86 + G_87 + G_88 + G_89 + G_90 + G_91 + G_92 + G_93 + G_94 + G_95 + G_96 + G_97 + G_98 + G_99 + G_100,data=top100_posteriors_train,family="binomial"), ~.,trace=1000)
selected_validate_aic12 <- c(22,36,50,67,84)
selected_train_aic22 <- c(8,43,1,24,29,35,75,80,37,62) # top10
selected_train_aic12 <- c(1,8,24,43,62)
selected_train_aic6 <- c(8,43)


# predict using selected genes
library(Brobdingnag)
library(pROC)
selected <- selected_train_aic12
prod_loglik_an_an <- rep(1,length(ANs_toPredict))
prod_loglik_an_t <- rep(1,length(ANs_toPredict))
prod_loglik_t_an <- rep(1,length(Ts_toPredict))
prod_loglik_t_t <- rep(1,length(Ts_toPredict))
performances_PGM_selected <- matrix(ncol=6,nrow=length(selected))
colnames(performances_PGM_selected) <- c("sensitivity","specificity","AUC","runningNaiveBayes_sen","runningNaiveBayes_spe","AUC")
selected <- 1:100
for (i in 1:length(selected)) {
	performances_PGM_selected[i,1] <- length(which(predictions[[selected[i]]][Ts_toPredict,2] >= 0.5))/length(Ts_toPredict)
	performances_PGM_selected[i,2] <- length(which(predictions[[selected[i]]][ANs_toPredict,2] < 0.5))/length(ANs_toPredict)
	
	performances_PGM_selected[i,3] <- auc(predictor=c(predictions[[selected[i]]][Ts_toPredict,2],predictions[[selected[i]]][ANs_toPredict,2]),response=c(rep("Pos",length(Ts_toPredict)),rep("Neg",length(ANs_toPredict))))
	
	# naive Bayes combination of results
	prod_loglik_t_t <- prod_loglik_t_t * as.brob(exp(-predictions[[selected[i]]][Ts_toPredict,3]))
	prod_loglik_t_an <- prod_loglik_t_an * as.brob(exp(-predictions[[selected[i]]][Ts_toPredict,4]))
	
	prod_loglik_an_an <- prod_loglik_an_an * as.brob(exp(-predictions[[selected[i]]][ANs_toPredict,4]))
	prod_loglik_an_t <- prod_loglik_an_t * as.brob(exp(-predictions[[selected[i]]][ANs_toPredict,3]))
	
	performances_PGM_selected[i,5] <- length(which(as.double(prod_loglik_an_t / (prod_loglik_an_t+prod_loglik_an_an)) < 0.5))/length(ANs_toPredict)
	performances_PGM_selected[i,4] <- length(which(as.double(prod_loglik_t_t / (prod_loglik_t_t+prod_loglik_t_an)) >= 0.5))/length(Ts_toPredict)
	
	performances_PGM_selected[i,6] <- auc(predictor=c(as.double(prod_loglik_t_t / (prod_loglik_t_t+prod_loglik_t_an)),as.double(prod_loglik_an_t / (prod_loglik_an_t+prod_loglik_an_an))),response=c(rep("Pos",length(Ts_toPredict)),rep("Neg",length(ANs_toPredict))))
}



# predict using all significant genes in the progression dataset
library(Brobdingnag)
library(pROC)
performances_PGM <- matrix(ncol=6,nrow=length(predictions))
colnames(performances_PGM) <- c("sensitivity","specificity","AUC","runningNaiveBayes_sen","runningNaiveBayes_spe","AUC")

prod_loglik_G1_G2 <- rep(1,length(G1))
prod_loglik_G1_G1 <- rep(1,length(G1))

prod_loglik_G2_G2 <- rep(1,length(G2))
prod_loglik_G2_G1 <- rep(1,length(G2))

for (i in 1:length(predictions)) {
	performances_PGM[i,1] <- length(which(as.brob(exp(1))^-predictions[[i]][G1,3] / (as.brob(exp(1))^-predictions[[i]][G1,3] + as.brob(exp(1))^-predictions[[i]][G1,4]) >= 0.5))/length(G1)
	performances_PGM[i,2] <- length(which(as.brob(exp(1))^-predictions[[i]][G2,3] / (as.brob(exp(1))^-predictions[[i]][G2,3] + as.brob(exp(1))^-predictions[[i]][G2,4]) < 0.5))/length(G2)
	
	performances_PGM[i,3] <- auc(predictor=as.double(as.brob(exp(1))^-predictions[[i]][,3] / (as.brob(exp(1))^-predictions[[i]][,3] + as.brob(exp(1))^-predictions[[i]][,4])),response=c(rep("Pos",length(G1)),rep("Neg",length(G2))))
	
	# naive Bayes combination of results
	prod_loglik_G1_G1 <- prod_loglik_G1_G1 * as.brob(exp(1))^-predictions[[i]][G1,3]
	prod_loglik_G1_G2 <- prod_loglik_G1_G2 * as.brob(exp(1))^-predictions[[i]][G1,4]
	
	prod_loglik_G2_G2 <- prod_loglik_G2_G2 * as.brob(exp(1))^-predictions[[i]][G2,4]
	prod_loglik_G2_G1 <- prod_loglik_G2_G1 * as.brob(exp(1))^-predictions[[i]][G2,3]
	
	performances_PGM[i,4] <- length(which(as.double(prod_loglik_G1_G1 / (prod_loglik_G1_G1+prod_loglik_G1_G2)) >= 0.5))/length(G1)
	performances_PGM[i,5] <- length(which(as.double(prod_loglik_G2_G1 / (prod_loglik_G2_G1+prod_loglik_G2_G2)) < 0.5))/length(G2)
	
	performances_PGM[i,6] <- auc(predictor=c(as.double(prod_loglik_G1_G1 / (prod_loglik_G1_G1+prod_loglik_G1_G2)),as.double(prod_loglik_G2_G1 / (prod_loglik_G2_G1+prod_loglik_G2_G2))),response=c(rep("Pos",length(G1)),rep("Neg",length(G2))))
}

