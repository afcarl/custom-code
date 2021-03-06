# dead within first 1065 days
c("TCGA-A7-A3RF","TCGA-A7-A425","TCGA-LL-A5YM","TCGA-E9-A243","TCGA-A7-A13G","TCGA-A7-A26H","TCGA-LQ-A4E4","TCGA-A7-A13H","TCGA-A8-A08O","TCGA-E9-A226","TCGA-A2-A3XY","TCGA-E9-A2JS","TCGA-A2-A3XU","TCGA-AR-A5QQ")
# alive within after first 1065 days
c("TCGA-A7-A0CE","TCGA-A7-A0CH","TCGA-E9-A1RI","TCGA-E9-A1NE","TCGA-OL-A5RW","TCGA-E9-A1NA","TCGA-E9-A1N5","TCGA-A7-A0D9","TCGA-AR-A1AS","TCGA-AR-A2LN","TCGA-GM-A3NY","TCGA-A2-A3Y0","TCGA-E9-A22A","TCGA-AR-A2LO","TCGA-E9-A1NC","TCGA-A2-A3KD","TCGA-AR-A2LQ","TCGA-AC-A2FB","TCGA-GM-A3XG","TCGA-A2-A0YL","TCGA-A2-A3XW","TCGA-BH-A0HY","TCGA-EW-A2FS","TCGA-EW-A1P3","TCGA-BH-A0HA","TCGA-EW-A2FR","TCGA-AR-A255","TCGA-AR-A1AV","TCGA-AR-A2LJ","TCGA-AR-A1AX","TCGA-AR-A1AM","TCGA-AR-A2LJ","TCGA-AR-A1AW","TCGA-OL-A66J","TCGA-GM-A3XN","TCGA-GM-A3XL","TCGA-GM-A4E0","TCGA-AR-A254","TCGA-AR-A252","TCGA-AR-A24T","TCGA-AR-A1AU","TCGA-AR-A251","TCGA-AR-A24N","TCGA-AR-A24Z","TCGA-A2-A3XT","TCGA-AR-A24X","TCGA-AR-A24V","TCGA-B6-A401","TCGA-AR-A0U4","TCGA-AR-A0TT","TCGA-AR-A24R","TCGA-AR-A24M","TCGA-AR-A0TW","TCGA-AR-A24Q","TCGA-B6-A40B","TCGA-A2-A0EP","TCGA-A2-A0CR","TCGA-GM-A3NW","TCGA-AR-A0TP","TCGA-AR-A0U3","TCGA-AQ-A04L")

pval_adj[sort(pval_adj,index.return=TRUE)$ix[1:100]]
ids <- signif_progression <- c(10461,13744,162,372,297,8004,16709,7437,1947,16178,5379,5063,5076,3951,6008,15307,3504,2882,9783,10903,12935,16995,811,10911,12001,14625,9181,3125,8753,1155,2240,5174,9012,14959,16905,5756,5517,3370,7065,4496,1402,7808,9860,9305,2496,16613,2379,3237,16153,13145,3671,653,14897,1047,7901,11720,5925,3079,2401,14582,16948,8175,15702,13547,17195,13700,12661,304,5955,1470,468,17264,17265,5652,165,6440,8552,1447,2546,9890,16557,1376,17433,17499,10132,928,10118,10099,2318,6788,15235,4238,15226,169,17214,8171,14171,697,7640,14861,13836,1428,6500,1522,4591,8572,14373,13854,5856,2782,360,14558,12518,6844,16054,17188,14219,4788,13119,3401,3934,4423,5260,14029,13543,17408,3356,8065,2020,16354,8735,14620,9888,12179,3678,5290,1559,11460,15168,6567,3363,17427,13192,6749,7756,10090,4682,12106,8045,5834,15819,1658,6267,7466,2215,7453,2936,10667,11345,969,2373,13855,15031,7992,3774,11420,3389,13028,7253,1141,14743,15590,10868,11122,16277,11931,178,15650,5946,17622,6164,3590,12742,2640,6322,7306,11577,9653,16093,9944,2081,3326,11614,2929,14334,7444,5429,2422,3677,16194,3346,9557,3198,11733,4313)

# read predictions from files produced by 2way_predict.R
predictions <- NULL
for (i in 1:length(ids)) {
	predictions[[i]] <- read.table(paste(ids[i],".predicted",sep="",collapse=""),header=TRUE)
	rownames(predictions[[i]]) <- predictions[[i]][,1]
}

# recalculate failed posterior calculations on own machine using brobs

library(Brobdingnag)
for (i in 1:length(predictions_ML)) {
	predictions_ML[[i]][,2] <- as.double(as.brob(exp(1))^-predictions_ML[[i]][,3] / (as.brob(exp(1))^-predictions_ML[[i]][,3] + as.brob(exp(1))^-predictions_ML[[i]][,4]))
}

# selection of most predictive genes among the all 205 stat. signif. genes using the GLM's logistic regression
top205_posteriors_train <- matrix(ncol=206,nrow=length(c(G1,G2)))
for (i in 1:length(predictions_ML)) {
	top205_posteriors_train[,1+i] <- predictions_ML[[i]][,2]
}
top205_posteriors_train[,1] <- c(rep("G1",length(G1)),rep("G2",length(G2)))
rownames(top205_posteriors_train) <- c(G1,G2)
colnames(top205_posteriors_train) <- c("class",paste("G_",1:205,sep=""))
top205_posteriors_train <- as.data.frame(top205_posteriors_train)
as.numeric.factor <- function(x) {as.numeric(levels(x))[x]}
for (i in 1:205) top205_posteriors_train[,1+i] <- as.numeric.factor(top205_posteriors_train[,1+i])

search <- step(glm(class ~ G_1 + G_2 + G_3 + G_4 + G_5 + G_6 + G_7 + G_8 + G_9 + G_10 + G_11 + G_12 + G_13 + G_14 + G_15 + G_16 + G_17 + G_18 + G_19 + G_20 + G_21 + G_22 + G_23 + G_24 + G_25 + G_26 + G_27 + G_28 + G_29 + G_30 + G_31 + G_32 + G_33 + G_34 + G_35 + G_36 + G_37 + G_38 + G_39 + G_40 + G_41 + G_42 + G_43 + G_44 + G_45 + G_46 + G_47 + G_48 + G_49 + G_50 + G_51 + G_52 + G_53 + G_54 + G_55 + G_56 + G_57 + G_58 + G_59 + G_60 + G_61 + G_62 + G_63 + G_64 + G_65 + G_66 + G_67 + G_68 + G_69 + G_70 + G_71 + G_72 + G_73 + G_74 + G_75 + G_76 + G_77 + G_78 + G_79 + G_80 + G_81 + G_82 + G_83 + G_84 + G_85 + G_86 + G_87 + G_88 + G_89 + G_90 + G_91 + G_92 + G_93 + G_94 + G_95 + G_96 + G_97 + G_98 + G_99 + G_100 + G_101 + G_102 + G_103 + G_104 + G_105 + G_106 + G_107 + G_108 + G_109 + G_110 + G_111 + G_112 + G_113 + G_114 + G_115 + G_116 + G_117 + G_118 + G_119 + G_120 + G_121 + G_122 + G_123 + G_124 + G_125 + G_126 + G_127 + G_128 + G_129 + G_130 + G_131 + G_132 + G_133 + G_134 + G_135 + G_136 + G_137 + G_138 + G_139 + G_140 + G_141 + G_142 + G_143 + G_144 + G_145 + G_146 + G_147 + G_148 + G_149 + G_150 + G_151 + G_152 + G_153 + G_154 + G_155 + G_156 + G_157 + G_158 + G_159 + G_160 + G_161 + G_162 + G_163 + G_164 + G_165 + G_166 + G_167 + G_168 + G_169 + G_170 + G_171 + G_172 + G_173 + G_174 + G_175 + G_176 + G_177 + G_178 + G_179 + G_180 + G_181 + G_182 + G_183 + G_184 + G_185 + G_186 + G_187 + G_188 + G_189 + G_190 + G_191 + G_192 + G_193 + G_194 + G_195 + G_196 + G_197 + G_198 + G_199 + G_200 + G_201 + G_202 + G_203 + G_204 + G_205 ,data=top205_posteriors_train,family="binomial"), ~.,trace=1000)

selected_ML <- c(60,7,13,57,14,201,19,47,104,18,45,168,172,110,15,27,129,53,109,67)
selected <- selected_ML

# predict using all selected genes
library(Brobdingnag)
library(pROC)
performances_MLfound <- matrix(ncol=6,nrow=length(selected_ML))
colnames(performances_MLfound) <- c("sensitivity","specificity","AUC","runningNaiveBayes_sen","runningNaiveBayes_spe","AUC")

prod_loglik_G1_G2 <- rep(1,length(G1))
prod_loglik_G1_G1 <- rep(1,length(G1))

prod_loglik_G2_G2 <- rep(1,length(G2))
prod_loglik_G2_G1 <- rep(1,length(G2))

for (i in 1:length(selected_ML)) {
	performances_MLfound[i,1] <- length(which(as.brob(exp(1))^-predictions_LOU[[selected[i]]][G1,3] / (as.brob(exp(1))^-predictions_LOU[[selected[i]]][G1,3] + as.brob(exp(1))^-predictions_LOU[[selected[i]]][G1,4]) >= 0.5))/length(G1)
	performances_MLfound[i,2] <- length(which(as.brob(exp(1))^-predictions_LOU[[selected[i]]][G2,3] / (as.brob(exp(1))^-predictions_LOU[[selected[i]]][G2,3] + as.brob(exp(1))^-predictions_LOU[[selected[i]]][G2,4]) < 0.5))/length(G2)
	
	performances_MLfound[i,3] <- auc(predictor=as.double(as.brob(exp(1))^-predictions_LOU[[selected[i]]][,3] / (as.brob(exp(1))^-predictions_LOU[[selected[i]]][,3] + as.brob(exp(1))^-predictions_LOU[[selected[i]]][,4])),response=c(rep("Pos",length(G1)),rep("Neg",length(G2))))
	
	# naive Bayes combination of results
	prod_loglik_G1_G1 <- prod_loglik_G1_G1 * as.brob(exp(1))^-predictions_LOU[[selected[i]]][G1,3]
	prod_loglik_G1_G2 <- prod_loglik_G1_G2 * as.brob(exp(1))^-predictions_LOU[[selected[i]]][G1,4]
	
	prod_loglik_G2_G2 <- prod_loglik_G2_G2 * as.brob(exp(1))^-predictions_LOU[[selected[i]]][G2,4]
	prod_loglik_G2_G1 <- prod_loglik_G2_G1 * as.brob(exp(1))^-predictions_LOU[[selected[i]]][G2,3]
	
	performances_MLfound[i,4] <- length(which(as.double(prod_loglik_G1_G1 / (prod_loglik_G1_G1+prod_loglik_G1_G2)) >= 0.5))/length(G1)
	performances_MLfound[i,5] <- length(which(as.double(prod_loglik_G2_G1 / (prod_loglik_G2_G1+prod_loglik_G2_G2)) < 0.5))/length(G2)
	
	performances_MLfound[i,6] <- auc(predictor=c(as.double(prod_loglik_G1_G1 / (prod_loglik_G1_G1+prod_loglik_G1_G2)),as.double(prod_loglik_G2_G1 / (prod_loglik_G2_G1+prod_loglik_G2_G2))),response=c(rep("Pos",length(G1)),rep("Neg",length(G2))))
}



# read x-val data

Zs <- matrix(nrow=17728,ncol=7)
mlogliks <- NULL
for (i in 1:17728) {
	Zs[i,] <- t(read.table(paste("./",i,".result",sep=""),nrow=1,skip=0))
	mlogliks[[i]] <- read.table(paste("./",i,".result",sep=""),skip=1,header=TRUE)
}
rm(i)

# print top20 genes across folds
top20_genes_perFold_14 <- matrix(ncol=15,nrow=20)
for (i in 1:14) top20_genes_perFold_14[,i+1] <- workingList_BRCA[sort(rank(-Zs[,i]),index.return=TRUE)$ix[1:20]]

# subset and join the top20 x-val logliks
mlogliks_top20_sorted_14 <- mlogliks[1:20]
for (fold in 0:13) {
	index_G1_fold <- (fold+1):(fold+1)
	index_G2_fold <- (fold*4+1+14):((fold+1)*4+14)
	if (fold == 13) index_G2_fold <- c(index_G2_fold,71)
	
	top20i <- sort(rank(-Zs[,fold+1]),index.return=TRUE)$ix[1:20]
	for (top in 1:length(top20i)) {
		mlogliks_top20_sorted_14[[top]][c(index_G1_fold,index_G2_fold),] <- mlogliks[[top20i[top]]][c(index_G1_fold,index_G2_fold),]
	}
}
for (i in 1:20) rownames(mlogliks_top20_sorted_14[[i]]) <- mlogliks_top20_sorted_14[[i]][,1]

library(Brobdingnag)
library(pROC)
# calculate performances
mlogliks_top20_sorted <- mlogliks_top20_sorted_14

performances_XVal <- matrix(ncol=6,nrow=20)
colnames(performances_XVal) <- c("sensitivity","specificity","AUC","runningNaiveBayes_sen","runningNaiveBayes_spe","AUC")
prod_loglik_G1_G2 <- rep(1,length(G1))
prod_loglik_G1_G1 <- rep(1,length(G1))
prod_loglik_G2_G2 <- rep(1,length(G2))
prod_loglik_G2_G1 <- rep(1,length(G2))

for (i in 1:nrow(performances_XVal)) {
	performances_XVal[i,1] <- length(which(as.brob(exp(1))^-mlogliks_top20_sorted[[i]][G1,2] / (as.brob(exp(1))^-mlogliks_top20_sorted[[i]][G1,2] + as.brob(exp(1))^-mlogliks_top20_sorted[[i]][G1,3]) >= 0.5))/length(G1)
	performances_XVal[i,2] <- length(which(as.brob(exp(1))^-mlogliks_top20_sorted[[i]][G2,2] / (as.brob(exp(1))^-mlogliks_top20_sorted[[i]][G2,2] + as.brob(exp(1))^-mlogliks_top20_sorted[[i]][G2,3]) < 0.5))/length(G2)
	
	performances_XVal[i,3] <- auc(predictor=as.double(as.brob(exp(1))^-mlogliks_top20_sorted[[i]][,2] / (as.brob(exp(1))^-mlogliks_top20_sorted[[i]][,2] + as.brob(exp(1))^-mlogliks_top20_sorted[[i]][,3])),response=c(rep("Pos",length(G1)),rep("Neg",length(G2))))
	
	# naive Bayes combination of results
	prod_loglik_G1_G1 <- prod_loglik_G1_G1 * as.brob(exp(1))^-mlogliks_top20_sorted[[i]][G1,2]
	prod_loglik_G1_G2 <- prod_loglik_G1_G2 * as.brob(exp(1))^-mlogliks_top20_sorted[[i]][G1,3]
	
	prod_loglik_G2_G2 <- prod_loglik_G2_G2 * as.brob(exp(1))^-mlogliks_top20_sorted[[i]][G2,3]
	prod_loglik_G2_G1 <- prod_loglik_G2_G1 * as.brob(exp(1))^-mlogliks_top20_sorted[[i]][G2,2]
	
	performances_XVal[i,4] <- length(which(as.double(prod_loglik_G1_G1 / (prod_loglik_G1_G1+prod_loglik_G1_G2)) >= 0.5))/length(G1)
	performances_XVal[i,5] <- length(which(as.double(prod_loglik_G2_G1 / (prod_loglik_G2_G1+prod_loglik_G2_G2)) < 0.5))/length(G2)
	
	performances_XVal[i,6] <- auc(predictor=c(as.double(prod_loglik_G1_G1 / (prod_loglik_G1_G1+prod_loglik_G1_G2)),as.double(prod_loglik_G2_G1 / (prod_loglik_G2_G1+prod_loglik_G2_G2))),response=c(rep("Pos",length(G1)),rep("Neg",length(G2))))
}


# full 14-fold CV: the good way
library(edgeR)
library(pROC)

indiv_predictions_cv <- matrix(ncol=71,nrow=10)
colnames(indiv_predictions_cv) <- c(G1,G2)
indiv_predictions_cv_r <- matrix(ncol=71,nrow=10)
colnames(indiv_predictions_cv_r) <- c(G1,G2)
top10genes_xval_E <- matrix(ncol=14,nrow=10)
for (fold in 0:13) {
	index_G1_fold <- (fold+1):(fold+1)
	index_G2_fold <- (fold*4+1):((fold+1)*4)
	if (fold == 13) index_G2_fold <- c(index_G2_fold,57)
	
	y_cv <- DGEList(counts=counts[,c(G1[-index_G1_fold],G2[-index_G2_fold])],group=targets$Type[-c(index_G1_fold,14+index_G2_fold)])
	y_cv$samples$lib.size <- factors_ls[c(G1[-index_G1_fold],G2[-index_G2_fold])]*10^6
	y_cv <- estimateTagwiseDisp(calcNormFactors(y_cv))
	et_cv <- exactTest(y_cv,dispersion="auto")
	top10_cv <- rownames(et_cv[order(et_cv$table[,3])[1:10],])
	top10genes_xval_E[,fold+1] <- top10_cv
	
	temp_train <- as.data.frame(t(cpm[top10_cv,c(G1[-index_G1_fold],G2[-index_G2_fold])]))
	temp_train$class <- factor(c(rep("G1",length(G1[-index_G1_fold])),rep("G2",length(G2[-index_G2_fold]))))
	temp_eval <- as.data.frame(t(cpm[top10_cv,c(G1[index_G1_fold],G2[index_G2_fold])]))
	temp_eval$class <- factor(c(rep("G1",length(G1[index_G1_fold])),rep("G2",length(G2[index_G2_fold]))))
	for (i in 1:10) {
		model <- glm(as.formula(paste("class ~ ",top10_cv[i])), data=temp_train, family="binomial")
		indiv_predictions_cv[i,c(G1[index_G1_fold],G2[index_G2_fold])] <- predict(model,temp_eval,type="response")
		model <- glm(as.formula(paste("class ~ ", paste(top10_cv[1:i],collapse=" + "))), data=temp_train, family="binomial")
		indiv_predictions_cv_r[i,c(G1[index_G1_fold],G2[index_G2_fold])] <- predict(model,temp_eval,type="response")
	}
}

performances_binomial_fullCV_top10E <- matrix(ncol=6,nrow=10)
colnames(performances_binomial_fullCV_top10E) <- c("sensitivity","specificity","AUC","sensitivity_r","specificity_r","AUC_r")

for (i in 1:10) {
	performances_binomial_fullCV_top10E[i,1] <- length(which(indiv_predictions_cv[i,1:length(G1)] < 0.5))/length(G1)
	performances_binomial_fullCV_top10E[i,2] <- length(which(indiv_predictions_cv[i,(1+length(G1)):(length(G1)+length(G2))] >= 0.5))/length(G2)
	performances_binomial_fullCV_top10E[i,3] <- auc(predictor=indiv_predictions_cv[i,],response=c(rep("Pos",length(G1)),rep("Neg",length(G2))))
	performances_binomial_fullCV_top10E[i,4] <- length(which(indiv_predictions_cv_r[i,1:length(G1)] < 0.5))/length(G1)
	performances_binomial_fullCV_top10E[i,5] <- length(which(indiv_predictions_cv_r[i,(1+length(G1)):(length(G1)+length(G2))] >= 0.5))/length(G2)
	performances_binomial_fullCV_top10E[i,6] <- auc(predictor=indiv_predictions_cv_r[i,],response=c(rep("Pos",length(G1)),rep("Neg",length(G2))))
}

###############################################################
# methylation only classifiers full 14CV
library(IMA)
library(pROC)

indiv_predictions_cv <- matrix(ncol=71,nrow=10)
colnames(indiv_predictions_cv) <- c(G1,G2)
indiv_predictions_cv_r <- matrix(ncol=71,nrow=10)
colnames(indiv_predictions_cv_r) <- c(G1,G2)
top10genes_xval_M <- matrix(ncol=14,nrow=10)
for (fold in 0:13) {
	index_G1_fold <- (fold+1):(fold+1)
	index_G2_fold <- (fold*4+1):((fold+1)*4)
	if (fold == 13) index_G2_fold <- c(index_G2_fold,57)
	
	temp_res <- as.data.frame(matrix(ncol=3,nrow=2*length(workingList_BRCA)))
	
	pr_cv <- sitetest2(mmatrix_BRCA_PROMOTER_noncap[,c(G1[-index_G1_fold],G2[-index_G2_fold])],targets[-c(index_G1_fold,14+index_G2_fold),],gcase="G1",gcontrol = "G2",testmethod = "satterthwaite")
	temp_res[1:length(workingList_BRCA),] <- cbind(workingList_BRCA,"promoter",pr_cv[,1])
	gb_cv <- sitetest2(mmatrix_BRCA_BODY_noncap[,c(G1[-index_G1_fold],G2[-index_G2_fold])],targets[-c(index_G1_fold,14+index_G2_fold),],gcase="G1",gcontrol = "G2",testmethod = "satterthwaite")
	temp_res[(1+length(workingList_BRCA)):(2*length(workingList_BRCA)),] <- cbind(workingList_BRCA,"body",gb_cv[,1])
	temp_res[,3] <- as.numeric(temp_res[,3])
	top10_cv <- temp_res[order(temp_res[,3])[1:10],]
	top10genes_xval_M[,fold+1] <- top10_cv[,1]
	
	temp_train <- as.data.frame(matrix(nrow=length(c(G1[-index_G1_fold],G2[-index_G2_fold])),ncol=10))
	
	temp_eval <- as.data.frame(matrix(nrow=length(c(G1[index_G1_fold],G2[index_G2_fold])),ncol=10))
	colnames(temp_train) <- paste("G_",1:10,sep="")
	colnames(temp_eval) <- paste("G_",1:10,sep="")
	gene_symbols <- paste("G_",1:10,sep="")
	for (i in 1:10) {
		if (top10_cv[i,2] == "promoter") {
			temp_train[,i] <- mmatrix_BRCA_PROMOTER_noncap[top10_cv[i,1], c(G1[-index_G1_fold], G2[-index_G2_fold])]
			temp_eval[,i] <- mmatrix_BRCA_PROMOTER_noncap[top10_cv[i,1], c(G1[index_G1_fold], G2[index_G2_fold])]
		}
		if (top10_cv[i,2] == "body") {
			temp_train[,i] <- mmatrix_BRCA_BODY_noncap[top10_cv[i,1], c(G1[-index_G1_fold], G2[-index_G2_fold])]
			temp_eval[,i] <- mmatrix_BRCA_BODY_noncap[top10_cv[i,1], c(G1[index_G1_fold], G2[index_G2_fold])]
		}
	}
	temp_eval$class <- factor(c(rep("G1",length(G1[index_G1_fold])),rep("G2",length(G2[index_G2_fold]))))
	temp_train$class <- factor(c(rep("G1",length(G1[-index_G1_fold])),rep("G2",length(G2[-index_G2_fold]))))
	
	for (i in 1:10) {
		model <- glm(as.formula(paste("class ~ ",gene_symbols[i])), data=temp_train, family="binomial")
		indiv_predictions_cv[i,c(G1[index_G1_fold],G2[index_G2_fold])] <- predict(model,temp_eval,type="response")
		model <- glm(as.formula(paste("class ~ ", paste(gene_symbols[1:i],collapse=" + "))), data=temp_train, family="binomial")
		indiv_predictions_cv_r[i,c(G1[index_G1_fold],G2[index_G2_fold])] <- predict(model,temp_eval,type="response")
	}
}

performances_binomial_fullCV_top10M <- matrix(ncol=6,nrow=10)
colnames(performances_binomial_fullCV_top10M) <- c("sensitivity","specificity","AUC","sensitivity_r","specificity_r","AUC_r")

for (i in 1:10) {
	performances_binomial_fullCV_top10M[i,1] <- length(which(indiv_predictions_cv[i,1:length(G1)] < 0.5))/length(G1)
	performances_binomial_fullCV_top10M[i,2] <- length(which(indiv_predictions_cv[i,(1+length(G1)):(length(G1)+length(G2))] >= 0.5))/length(G2)
	performances_binomial_fullCV_top10M[i,3] <- auc(predictor=indiv_predictions_cv[i,],response=c(rep("Pos",length(G1)),rep("Neg",length(G2))))
	performances_binomial_fullCV_top10M[i,4] <- length(which(indiv_predictions_cv_r[i,1:length(G1)] < 0.5))/length(G1)
	performances_binomial_fullCV_top10M[i,5] <- length(which(indiv_predictions_cv_r[i,(1+length(G1)):(length(G1)+length(G2))] >= 0.5))/length(G2)
	performances_binomial_fullCV_top10M[i,6] <- auc(predictor=indiv_predictions_cv_r[i,],response=c(rep("Pos",length(G1)),rep("Neg",length(G2))))
}

###############################################################


# all 3 variables in the model
indiv_predictions_cv <- matrix(ncol=71,nrow=10)
colnames(indiv_predictions_cv) <- c(G1,G2)
indiv_predictions_cv_r <- matrix(ncol=71,nrow=10)
colnames(indiv_predictions_cv_r) <- c(G1,G2)
top10genes_xval_3var <- matrix(ncol=14,nrow=10)
for (fold in 0:13) {
	index_G1_fold <- (fold+1):(fold+1)
	index_G2_fold <- (fold*4+1):((fold+1)*4)
	if (fold == 13) index_G2_fold <- c(index_G2_fold,57)
	
	y_cv <- DGEList(counts=counts[,c(G1[-index_G1_fold],G2[-index_G2_fold])],group=targets$Type[-c(index_G1_fold,14+index_G2_fold)])
	y_cv$samples$lib.size <- factors_ls[c(G1[-index_G1_fold],G2[-index_G2_fold])]*10^6
	y_cv <- estimateTagwiseDisp(calcNormFactors(y_cv))
	et_cv <- exactTest(y_cv,dispersion="auto")
	
	temp_res <- matrix(ncol=4,nrow=length(workingList_BRCA))
	rownames(temp_res) <- workingList_BRCA
	temp_res[,1] <- et_cv$table[workingList_BRCA,3]
	
	pr_cv <- sitetest2(bmatrix_BRCA_PROMOTER[,c(G1[-index_G1_fold],G2[-index_G2_fold])],targets[-c(index_G1_fold,14+index_G2_fold),],gcase="G1",gcontrol = "G2",testmethod = "satterthwaite")
	temp_res[,2] <- pr_cv[,1]
	gb_cv <- sitetest2(bmatrix_BRCA_BODY[,c(G1[-index_G1_fold],G2[-index_G2_fold])],targets[-c(index_G1_fold,14+index_G2_fold),],gcase="G1",gcontrol = "G2",testmethod = "satterthwaite")
	temp_res[,3] <- gb_cv[,1]
	
	temp_res[,4] <- apply(temp_res[,1:3],1,fishersMethod)
	
	top10_cv <- rownames(temp_res[order(temp_res[,4])[1:10],])
	top10genes_xval_3var[,fold+1] <- top10_cv
	temp_train <- as.data.frame(t(cpm[top10_cv,c(G1[-index_G1_fold],G2[-index_G2_fold])]))
	temp_train <- cbind(temp_train,as.data.frame(t(mmatrix_BRCA_PROMOTER_noncap[top10_cv,c(G1[-index_G1_fold],G2[-index_G2_fold])])))
	temp_train <- cbind(temp_train,as.data.frame(t(mmatrix_BRCA_BODY_noncap[top10_cv,c(G1[-index_G1_fold],G2[-index_G2_fold])])))
	colnames(temp_train)[1:10] <- paste(colnames(temp_train)[1:10],"_E",sep="")
	colnames(temp_train)[11:20] <- paste(colnames(temp_train)[11:20],"_PR",sep="")
	colnames(temp_train)[21:30] <- paste(colnames(temp_train)[21:30],"_GB",sep="")
	temp_train$class <- factor(c(rep("G1",length(G1[-index_G1_fold])),rep("G2",length(G2[-index_G2_fold]))))
	
	temp_eval <- as.data.frame(t(cpm[top10_cv,c(G1[index_G1_fold],G2[index_G2_fold])]))
	temp_eval <- cbind(temp_eval,as.data.frame(t(mmatrix_BRCA_PROMOTER_noncap[top10_cv,c(G1[index_G1_fold],G2[index_G2_fold])])))
	temp_eval <- cbind(temp_eval,as.data.frame(t(mmatrix_BRCA_BODY_noncap[top10_cv,c(G1[index_G1_fold],G2[index_G2_fold])])))
	colnames(temp_eval)[1:10] <- paste(colnames(temp_eval)[1:10],"_E",sep="")
	colnames(temp_eval)[11:20] <- paste(colnames(temp_eval)[11:20],"_PR",sep="")
	colnames(temp_eval)[21:30] <- paste(colnames(temp_eval)[21:30],"_GB",sep="")
	temp_eval$class <- factor(c(rep("G1",length(G1[index_G1_fold])),rep("G2",length(G2[index_G2_fold]))))
	for (i in 1:10) {
		model <- glm(as.formula(paste("class ~ ",paste(paste(top10_cv[i],"_E",sep=""),paste(top10_cv[i],"_PR",sep=""),paste(top10_cv[i],"_GB",sep=""),sep=" + "))), data=temp_train, family="binomial")
		indiv_predictions_cv[i,c(G1[index_G1_fold],G2[index_G2_fold])] <- predict(model,temp_eval,type="response")
		
		model <- glm(as.formula(paste("class ~ ",paste(paste(top10_cv[1:i],"_E",sep="", collapse=" + "), paste(top10_cv[1:i],"_PR",sep="", collapse=" + "), paste(top10_cv[1:i],"_GB",sep="", collapse=" + "),sep=" + "))), data=temp_train, family="binomial")
		indiv_predictions_cv_r[i,c(G1[index_G1_fold],G2[index_G2_fold])] <- predict(model,temp_eval,type="response")
	}
}

performances_binomial_fullCV_top103var <- matrix(ncol=6,nrow=10)
colnames(performances_binomial_fullCV_top103var) <- c("sensitivity","specificity","AUC","sensitivity_r","specificity_r","AUC_r")

for (i in 1:10) {
	performances_binomial_fullCV_top103var[i,1] <- length(which(indiv_predictions_cv[i,1:length(G1)] < 0.5))/length(G1)
	performances_binomial_fullCV_top103var[i,2] <- length(which(indiv_predictions_cv[i,(1+length(G1)):(length(G1)+length(G2))] >= 0.5))/length(G2)
	performances_binomial_fullCV_top103var[i,3] <- auc(predictor=indiv_predictions_cv[i,],response=c(rep("Pos",length(G1)),rep("Neg",length(G2))))
	performances_binomial_fullCV_top103var[i,4] <- length(which(indiv_predictions_cv_r[i,1:length(G1)] < 0.5))/length(G1)
	performances_binomial_fullCV_top103var[i,5] <- length(which(indiv_predictions_cv_r[i,(1+length(G1)):(length(G1)+length(G2))] >= 0.5))/length(G2)
	performances_binomial_fullCV_top103var[i,6] <- auc(predictor=indiv_predictions_cv_r[i,],response=c(rep("Pos",length(G1)),rep("Neg",length(G2))))
}

performances_binomial_fullCV_top103var_i <- matrix(ncol=3,nrow=10)
colnames(performances_binomial_fullCV_top103var_i) <- c("sensitivity","specificity","AUC")

for (i in 1:10) {
	performances_binomial_fullCV_top103var_i[i,1] <- length(which(indiv_predictions_cv[i,1:length(G1)] < 0.5))/length(G1)
	performances_binomial_fullCV_top103var_i[i,2] <- length(which(indiv_predictions_cv[i,(1+length(G1)):(length(G1)+length(G2))] >= 0.5))/length(G2)
	performances_binomial_fullCV_top103var_i[i,3] <- auc(predictor=indiv_predictions_cv[i,],response=c(rep("Pos",length(G1)),rep("Neg",length(G2))))
}


top20_genes_perFold_14 <- top20_genes_perFold_14[,-1]
# build logistic regression models using PGM-found genes at each fold

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
	
	for (i in 1:10) {
		model <- glm(as.formula(paste("class ~ ",paste(paste(current[i],"_E",sep=""),paste(current[i],"_PR",sep=""),paste(current[i],"_GB",sep=""),sep=" + "))), data=temp_train, family="binomial")
		indiv_predictions_pgmBased_cv[i,c(G1[index_G1_fold],G2[index_G2_fold])] <- predict(model,temp_eval,type="response")
		
		model <- glm(as.formula(paste("class ~ ",paste(paste(current[1:i],"_E",sep="", collapse=" + "), paste(current[1:i],"_PR",sep="", collapse=" + "), paste(current[1:i],"_GB",sep="", collapse=" + "),sep=" + "))), data=temp_train, family="binomial")
		indiv_predictions_pgmBased_cv_r[i,c(G1[index_G1_fold],G2[index_G2_fold])] <- predict(model,temp_eval,type="response")
	}
}

performances_pgmBased_LR_cv <- matrix(ncol=6,nrow=10)
colnames(performances_pgmBased_LR_cv) <- c("sensitivity","specificity","AUC","sensitivity_r","specificity_r","AUC_r")

for (i in 1:10) {
	performances_pgmBased_LR_cv[i,1] <- length(which(indiv_predictions_pgmBased_cv[i,1:length(G1)] < 0.5))/length(G1)
	performances_pgmBased_LR_cv[i,2] <- length(which(indiv_predictions_pgmBased_cv[i,(1+length(G1)):(length(G1)+length(G2))] >= 0.5))/length(G2)
	performances_pgmBased_LR_cv[i,3] <- auc(predictor=indiv_predictions_pgmBased_cv[i,],response=c(rep("Pos",length(G1)),rep("Neg",length(G2))))
	performances_pgmBased_LR_cv[i,4] <- length(which(indiv_predictions_pgmBased_cv_r[i,1:length(G1)] < 0.5))/length(G1)
	performances_pgmBased_LR_cv[i,5] <- length(which(indiv_predictions_pgmBased_cv_r[i,(1+length(G1)):(length(G1)+length(G2))] >= 0.5))/length(G2)
	performances_pgmBased_LR_cv[i,6] <- auc(predictor=indiv_predictions_pgmBased_cv_r[i,],response=c(rep("Pos",length(G1)),rep("Neg",length(G2))))
}

require(ggplot2)
require(pROC)
require(Brobdingnag)
# make and plot ROC
prod_loglik_G1_G2 <- rep(1,length(G1))
prod_loglik_G1_G1 <- rep(1,length(G1))
prod_loglik_G2_G2 <- rep(1,length(G2))
prod_loglik_G2_G1 <- rep(1,length(G2))

for (i in 1:2) {
	# naive Bayes combination of results
	prod_loglik_G1_G1 <- prod_loglik_G1_G1 * as.brob(exp(1))^-mlogliks_top20_sorted[[i]][G1,2]
	prod_loglik_G1_G2 <- prod_loglik_G1_G2 * as.brob(exp(1))^-mlogliks_top20_sorted[[i]][G1,3]
	
	prod_loglik_G2_G2 <- prod_loglik_G2_G2 * as.brob(exp(1))^-mlogliks_top20_sorted[[i]][G2,3]
	prod_loglik_G2_G1 <- prod_loglik_G2_G1 * as.brob(exp(1))^-mlogliks_top20_sorted[[i]][G2,2]
}

roc_pgm <- roc(cases=as.double(prod_loglik_G2_G2 / (prod_loglik_G2_G2+prod_loglik_G2_G1)), controls=as.double(prod_loglik_G1_G2 / (prod_loglik_G1_G2+prod_loglik_G1_G1)),ci=T)
roc_logistic <- roc(cases=indiv_predictions_pgmBased_cv_r[2,G1], controls=indiv_predictions_pgmBased_cv_r[2,G2],ci=T)
roc_logistic_9 <- roc(cases=indiv_predictions_pgmBased_cv_r[9,G1], controls=indiv_predictions_pgmBased_cv_r[9,G2],ci=T)

plotter <- data.frame(
Sensitivity=c(
	roc_pgm$sensitivities,
	roc_logistic$sensitivities,
	roc_logistic_9$sensitivities),
Specificity=c(
	roc_pgm$specificities,
	roc_logistic$specificities,
	roc_logistic_9$specificities),
Model=factor(c(
	rep("top2 PGMs combined",length(roc_pgm$specificities)),
	rep("top2 logistic regression combined",length(roc_logistic$specificities)),
	rep("top9 logistic regression combined",length(roc_logistic_9$specificities))),levels=c("top2 PGMs combined","top2 logistic regression combined","top9 logistic regression combined")))

plot2 <- ggplot(plotter,aes(x=Specificity,y=Sensitivity,colour=Model)) + geom_polygon(alpha=0) + scale_x_reverse() + theme_bw() + theme(legend.position="bottom") + geom_abline(intercept=1,slope=1,size=0.75) + scale_colour_brewer(palette="Set1", labels=paste(c("top2 PGM-combined: AUC=","top2 LR-combined: AUC=","top9 LR-combined: AUC="),c(signif(roc_pgm$auc,4), signif(roc_logistic$auc,4), signif(roc_logistic_9$auc,4))," (", c(capture.output(roc_pgm$ci), capture.output(roc_logistic$ci), capture.output(roc_logistic_9$ci)),")", sep = "")) + ggtitle("ROC plot for progression predictive models") + theme(plot.title = element_text(face="bold", size=14),  axis.title.x = element_text(face="bold", size=12), axis.title.y = element_text(face="bold", size=12, angle=90), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.justification=c(1,0), legend.position=c(1,0), legend.title=element_blank(), legend.key = element_blank(), legend.text = element_text(size=7))

pdf(file="progression ROCs comparison v3.pdf")
print(plot2)
dev.off()