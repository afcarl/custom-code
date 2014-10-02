bmatrix_BRCA_BODY2 <- matrix(ncol=150,nrow=17728,data=0)
colnames(bmatrix_BRCA_BODY2) <- c(ANs2,Ts2)
row.names(bmatrix_BRCA_BODY2) <- workingList_BRCA4
for (i in 1:17728) {
	IDs_body <- unique(c(eval(parse(text = paste('GENEBODYInd$SID$','"',workingList_BRCA4[i],'"',sep=""))), eval(parse(text = paste('UTR3Ind$SID$','"',workingList_BRCA4[i],'"',sep="")))))
	temp = bmatrix[IDs_body,c(ANs2,Ts2)]
	if (length(IDs_body)==1) {bmatrix_BRCA_BODY2[i, ] = temp} else {bmatrix_BRCA_BODY2[i, ] = apply(temp, 2, eval(mean), na.rm = TRUE)}
}

bmatrix_BRCA_PROMOTER2 <- matrix(ncol=150,nrow=17728,data=0)
colnames(bmatrix_BRCA_PROMOTER2) <- c(ANs2,Ts2)
row.names(bmatrix_BRCA_PROMOTER2) <- workingList_BRCA4
for (i in 1:17728) {
	IDs_promoter <- unique(c(eval(parse(text = paste('TSS1500Ind$SID$','"',workingList_BRCA4[i],'"',sep=""))),eval(parse(text = paste('TSS200Ind$SID$','"',workingList_BRCA4[i],'"',sep=""))), eval(parse(text = paste('UTR5Ind$SID$','"',workingList_BRCA4[i],'"',sep="")))),eval(parse(text = paste('EXON1Ind$SID$','"',workingList_BRCA4[i],'"',sep=""))))
	temp = bmatrix[IDs_promoter,c(ANs2,Ts2)]
	if (length(IDs_promoter)==1) {bmatrix_BRCA_PROMOTER2[i, ] = temp} else {bmatrix_BRCA_PROMOTER2[i, ] = apply(temp, 2, eval(mean), na.rm = TRUE)}
}

top20 <- rownames(counts_BRCA_top20)

mmatrix_BRCA_BODY_top20 <- matrix(ncol=length(samples),nrow=20,data=0)
colnames(mmatrix_BRCA_BODY_top20) <- samples
row.names(mmatrix_BRCA_BODY_top20) <- top20
for (i in 1:20) {
	IDs_body <- unique(c(eval(parse(text = paste('GENEBODYInd$SID$','"',top20[i],'"',sep=""))), eval(parse(text = paste('UTR3Ind$SID$','"',top20[i],'"',sep="")))))
	temp = mmatrix_pc_top20[IDs_body,]
	if (length(IDs_body)==1) {mmatrix_BRCA_BODY_top20[i, ] = temp} else {mmatrix_BRCA_BODY_top20[i, ] = apply(temp, 2, eval(mean), na.rm = TRUE)}
}

mmatrix_BRCA_PROMOTER_top20 <- matrix(ncol=length(samples),nrow=20,data=0)
colnames(mmatrix_BRCA_PROMOTER_top20) <- samples
row.names(mmatrix_BRCA_PROMOTER_top20) <- top20
for (i in 1:20) {
	IDs_promoter <- unique(c(eval(parse(text = paste('TSS1500Ind$SID$','"',top20[i],'"',sep=""))),eval(parse(text = paste('TSS200Ind$SID$','"',top20[i],'"',sep=""))), eval(parse(text = paste('UTR5Ind$SID$','"',top20[i],'"',sep="")))),eval(parse(text = paste('EXON1Ind$SID$','"',top20[i],'"',sep=""))))
	temp = mmatrix_pc_top20[IDs_promoter,]
	if (length(IDs_promoter)==1) {mmatrix_BRCA_PROMOTER_top20[i, ] = temp} else {mmatrix_BRCA_PROMOTER_top20[i, ] = apply(temp, 2, eval(mean), na.rm = TRUE)}
}
mmatrix_BRCA_PROMOTER_top20[which(mmatrix_BRCA_PROMOTER_top20 == -Inf)] <- -7
mmatrix_BRCA_PROMOTER_top20[which(mmatrix_BRCA_PROMOTER_top20 == Inf)] <- 7
mmatrix_BRCA_BODY_top20[which(mmatrix_BRCA_BODY_top20 == -Inf)] <- -7
mmatrix_BRCA_BODY_top20[which(mmatrix_BRCA_BODY_top20 == Inf)] <- 7