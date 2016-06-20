G1 <- 82
G2 <- 730
nfolds <- length(folds_g1)
length_g1 <- length(unlist(folds_g1))
length_g2 <- length(unlist(folds_g2))
xfold_LR_probabilities <- matrix(nrow=G1+G2,ncol=10, data=0) # a matrix to save sample class posteriors from respective folds, single gene models
xfold_LR_comb_probabilities <- matrix(nrow=G1+G2,ncol=10, data=0) # a matrix to save sample class probabilities from respective folds, running combinations of gene models
df <- data.frame(class = factor(c(rep("G1",G1), rep("G2",G2))))
df_all_data <- list() # a list with data frames to hold running data for each separate fold
for (j in 1:nfolds) {df_all_data[[j]] <- df}
for (j in 1:10) { # top 10 ranks
	for (fold in 1:nfolds) { # iterate through nfolds folds
		g1_xfold_training <- setdiff((1:length_g1),unlist(folds_g1[[fold]]))
		g2_xfold_training <- setdiff(length_g1+(1:length_g2),unlist(folds_g2[[fold]]))
		i <- top[j,fold]
		df <- data_BRCA[[i]] %>%
			as.data.frame() %>%
			mutate(EXPR = (read_count+1) / lib_size) %>%
			mutate(GB_overall = select(., starts_with("gb")) %>% rowMeans(., na.rm=T)) %>%
			mutate(PR_overall = select(., starts_with("pr")) %>% rowMeans(., na.rm=T)) %>%
			select(-starts_with("pr", ignore.case = F), -starts_with("gb", ignore.case = F), -lib_size, -read_count) %>%
			mutate(class = factor(c(rep("G1",G1), rep("G2",G2))))
			
		# single rank/fold LR gene model
		model <- glm(as.formula("class ~ EXPR + GB_overall + PR_overall"), data=df[c(g1_xfold_training,g2_xfold_training),], family="binomial")
		xfold_LR_probabilities[unlist(c(folds_g1[fold],folds_g2[fold])), j] <- predict(model,df[-c(g1_xfold_training,g2_xfold_training),],type="response")
		colnames(df) <- paste(colnames(df),j,sep="_")
		df_all_data[[fold]] <- cbind(df_all_data[[fold]], df[,1:3])
		
		# running combination of ranks/fold LR gene model
		model <- glm(as.formula(paste("class ~ ",paste(c(paste("EXPR",1:j,sep="_"),paste("GB_overall",1:j,sep="_"),paste("PR_overall",1:j,sep="_")),collapse =" + "))), data=df_all_data[[fold]][c(g1_xfold_training,g2_xfold_training),], family="binomial")
		xfold_LR_comb_probabilities[unlist(c(folds_g1[fold],folds_g2[fold])), j] <- predict(model,df_all_data[[fold]][-c(g1_xfold_training,g2_xfold_training),],type="response")
	}
}
# calculate AUCs
library(pROC)
aucs_single_lr <- vector() # single rank AUCs
aucs_running_combinations_lr <- vector() # running combination of posteriors across ranks AUCs
aucs_running_combinations_lr_comb <- vector()  # runningcombination of multi-gene features under one LR model AUCs
for (j in 1:10) {
	aucs_single_lr[j] <- auc(predictor=xfold_LR_probabilities[,j],response=c(rep("G1",G1),rep("G2",G2)))
	if (j > 1) aucs_running_combinations_lr[j] <- auc(predictor=apply(xfold_LR_probabilities[,1:j],1,prod), response=c(rep("G1",G1),rep("G2",G2)))
	aucs_running_combinations_lr_comb[j] <- auc(predictor=xfold_LR_comb_probabilities[,j], response=c(rep("G1",G1),rep("G2",G2)))
}
View(cbind(aucs_single_lr, aucs_running_combinations_lr, aucs_running_combinations_lr_comb))