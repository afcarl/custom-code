top20 <- rownames(counts_BRCA_top20_train)
samples <- c(ANs,Ts)

mmatrix_BRCA_BODY_top20_train <- matrix(ncol=length(samples),nrow=20,data=0)
colnames(mmatrix_BRCA_BODY_top20_train) <- samples
row.names(mmatrix_BRCA_BODY_top20_train) <- top20
for (i in 1:20) {
	IDs_body <- unique(c(eval(parse(text = paste('GENEBODYInd$SID$','"',top20[i],'"',sep=""))), eval(parse(text = paste('UTR3Ind$SID$','"',top20[i],'"',sep="")))))
	temp = mmatrix_pc[IDs_body,samples]
	if (length(IDs_body)==1) {mmatrix_BRCA_BODY_top20_train[i, ] = temp} else {mmatrix_BRCA_BODY_top20_train[i, ] = apply(temp, 2, eval(mean), na.rm = TRUE)}
}

mmatrix_BRCA_PROMOTER_top20_train <- matrix(ncol=length(samples),nrow=20,data=0)
colnames(mmatrix_BRCA_PROMOTER_top20_train) <- samples
row.names(mmatrix_BRCA_PROMOTER_top20_train) <- top20
for (i in 1:20) {
	IDs_promoter <- unique(c(eval(parse(text = paste('TSS1500Ind$SID$','"',top20[i],'"',sep=""))),eval(parse(text = paste('TSS200Ind$SID$','"',top20[i],'"',sep=""))), eval(parse(text = paste('UTR5Ind$SID$','"',top20[i],'"',sep="")))),eval(parse(text = paste('EXON1Ind$SID$','"',top20[i],'"',sep=""))))
	temp = mmatrix_pc[IDs_promoter,samples]
	if (length(IDs_promoter)==1) {mmatrix_BRCA_PROMOTER_top20_train[i, ] = temp} else {mmatrix_BRCA_PROMOTER_top20_train[i, ] = apply(temp, 2, eval(mean), na.rm = TRUE)}
}
mmatrix_BRCA_PROMOTER_top20_train[which(mmatrix_BRCA_PROMOTER_top20_train < -7)] <- -7
mmatrix_BRCA_PROMOTER_top20_train[which(mmatrix_BRCA_PROMOTER_top20_train > 7)] <- 7
mmatrix_BRCA_BODY_top20_train[which(mmatrix_BRCA_BODY_top20_train < -7)] <- -7
mmatrix_BRCA_BODY_top20_train[which(mmatrix_BRCA_BODY_top20_train > 7)] <- 7

samples <- c(ANs_toPredict,Ts_toPredict)


mmatrix_BRCA_BODY_top20_predict <- matrix(ncol=length(samples),nrow=20,data=0)
colnames(mmatrix_BRCA_BODY_top20_predict) <- samples
row.names(mmatrix_BRCA_BODY_top20_predict) <- top20
for (i in 1:20) {
	IDs_body <- unique(c(eval(parse(text = paste('GENEBODYInd$SID$','"',top20[i],'"',sep=""))), eval(parse(text = paste('UTR3Ind$SID$','"',top20[i],'"',sep="")))))
	temp = mmatrix_pc[IDs_body,samples]
	if (length(IDs_body)==1) {mmatrix_BRCA_BODY_top20_predict[i, ] = temp} else {mmatrix_BRCA_BODY_top20_predict[i, ] = apply(temp, 2, eval(mean), na.rm = TRUE)}
}

mmatrix_BRCA_PROMOTER_top20_predict <- matrix(ncol=length(samples),nrow=20,data=0)
colnames(mmatrix_BRCA_PROMOTER_top20_predict) <- samples
row.names(mmatrix_BRCA_PROMOTER_top20_predict) <- top20
for (i in 1:20) {
	IDs_promoter <- unique(c(eval(parse(text = paste('TSS1500Ind$SID$','"',top20[i],'"',sep=""))),eval(parse(text = paste('TSS200Ind$SID$','"',top20[i],'"',sep=""))), eval(parse(text = paste('UTR5Ind$SID$','"',top20[i],'"',sep="")))),eval(parse(text = paste('EXON1Ind$SID$','"',top20[i],'"',sep=""))))
	temp = mmatrix_pc[IDs_promoter,samples]
	if (length(IDs_promoter)==1) {mmatrix_BRCA_PROMOTER_top20_predict[i, ] = temp} else {mmatrix_BRCA_PROMOTER_top20_predict[i, ] = apply(temp, 2, eval(mean), na.rm = TRUE)}
}
mmatrix_BRCA_PROMOTER_top20_predict[which(mmatrix_BRCA_PROMOTER_top20_predict < -7)] <- -7
mmatrix_BRCA_PROMOTER_top20_predict[which(mmatrix_BRCA_PROMOTER_top20_predict > 7)] <- 7
mmatrix_BRCA_BODY_top20_predict[which(mmatrix_BRCA_BODY_top20_predict < -7)] <- -7
mmatrix_BRCA_BODY_top20_predict[which(mmatrix_BRCA_BODY_top20_predict > 7)] <- 7

cpm_BRCA_top20_predict <- counts_BRCA_plusOne_all[top20,samples]
for (i in 1:ncol(cpm_BRCA_top20_predict)) cpm_BRCA_top20_predict[,i] <- cpm_BRCA_top20_predict[,i]/factors_ls[i]

# data for predictions using GLMs
save(mmatrix_BRCA_PROMOTER_top20_predict,mmatrix_BRCA_BODY_top20_predict,cpm_BRCA_top20_predict,ANs_toPredict,Ts_toPredict,samples,file="predictData.RData")

# data for predictions using our PGMs
# find constituent probes
IDs <- NULL
for (i in 1:100) {
	IDs_body <- unique(c(eval(parse(text = paste('GENEBODYInd$SID$','"',top100[i],'"',sep=""))), eval(parse(text = paste('UTR3Ind$SID$','"',top100[i],'"',sep="")))))
	IDs_promoter <- unique(c(eval(parse(text = paste('TSS1500Ind$SID$','"',top100[i],'"',sep=""))),eval(parse(text = paste('TSS200Ind$SID$','"',top100[i],'"',sep=""))), eval(parse(text = paste('UTR5Ind$SID$','"',top100[i],'"',sep="")))),eval(parse(text = paste('EXON1Ind$SID$','"',top100[i],'"',sep=""))))
	IDs <- c(IDs,IDs_body,IDs_promoter)
}
IDs <- unique(IDs)
mmatrix_pc_top20 <- mmatrix_pc[IDs,samples]
mmatrix_pc_top20[which(mmatrix_pc_top20 < -7)] <- -7
mmatrix_pc_top20[which(mmatrix_pc_top20 > 7)] <- 7

counts_BRCA_top20 <- counts_BRCA_plusOne_all[top20,samples]

save(mmatrix_pc_top20,counts_BRCA_top20,ANs_toPredict,Ts_toPredict,factors_ls,file="BRCA_predictionData.RData")
