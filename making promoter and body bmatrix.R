mmatrix_BRCA_BODY <- matrix(ncol=length(c(ANs,Ts)),nrow=17728,data=0)
colnames(mmatrix_BRCA_BODY) <- c(ANs,Ts)
row.names(mmatrix_BRCA_BODY) <- workingList_BRCA
mmatrix_BRCA_PROMOTER <- matrix(ncol=length(c(ANs,Ts)),nrow=17728,data=0)
colnames(mmatrix_BRCA_PROMOTER) <- c(ANs,Ts)
row.names(mmatrix_BRCA_PROMOTER) <- workingList_BRCA
for (i in 1:17728) {
	IDs_body <- unique(c(eval(parse(text = paste('GENEBODYInd$SID$','"',workingList_BRCA[i],'"',sep=""))), eval(parse(text = paste('UTR3Ind$SID$','"',workingList_BRCA[i],'"',sep="")))))
	temp = mmatrix_pc[IDs_body,c(ANs,Ts)]
	if (length(IDs_body)==1) {mmatrix_BRCA_BODY[i, ] = temp} else {mmatrix_BRCA_BODY[i, ] = apply(temp, 2, eval(mean), na.rm = TRUE)}
	IDs_promoter <- unique(c(eval(parse(text = paste('TSS1500Ind$SID$','"',workingList_BRCA[i],'"',sep=""))),eval(parse(text = paste('TSS200Ind$SID$','"',workingList_BRCA[i],'"',sep=""))), eval(parse(text = paste('UTR5Ind$SID$','"',workingList_BRCA[i],'"',sep="")))),eval(parse(text = paste('EXON1Ind$SID$','"',workingList_BRCA[i],'"',sep=""))))
	temp = mmatrix_pc[IDs_promoter,c(ANs,Ts)]
	if (length(IDs_promoter)==1) {mmatrix_BRCA_PROMOTER[i, ] = temp} else {mmatrix_BRCA_PROMOTER[i, ] = apply(temp, 2, eval(mean), na.rm = TRUE)}
}

mmatrix_BRCA_PROMOTER[which(mmatrix_BRCA_PROMOTER < -7)] <- -7
mmatrix_BRCA_PROMOTER[which(mmatrix_BRCA_PROMOTER > 7)] <- 7
mmatrix_BRCA_BODY[which(mmatrix_BRCA_BODY < -7)] <- -7
mmatrix_BRCA_BODY[which(mmatrix_BRCA_BODY > 7)] <- 7
