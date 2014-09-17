data_BRCA <- list()
for (i in 1:10) {
	IDs_promoter <- unique(c(eval(parse(text = paste('TSS1500Ind$SID$','"',workingList_BRCA[i],'"',sep=""))),eval(parse(text = paste('TSS200Ind$SID$','"',workingList_BRCA[i],'"',sep=""))), eval(parse(text = paste('UTR5Ind$SID$','"',workingList_BRCA[i],'"',sep="")))),eval(parse(text = paste('EXON1Ind$SID$','"',workingList_BRCA[i],'"',sep=""))))
	IDs_body <- unique(c(eval(parse(text = paste('GENEBODYInd$SID$','"',workingList_BRCA[i],'"',sep=""))), eval(parse(text = paste('UTR3Ind$SID$','"',workingList_BRCA[i],'"',sep="")))))
	ncol = length(IDs_promoter) + length(IDs_body) + 2
	temp <- matrix(nrow=150,ncol=ncol)
	temp[,1] <- factors_ls[c(ANs,Ts)]*10^6
	temp[,2] <- counts_BRCA_plusOne[workingList_BRCA[i],c(ANs,Ts)]
	temp[,3:(2+length(IDs_promoter))] <- t(mmatrix[IDs_promoter,c(ANs,Ts)])
	temp[,(3+length(IDs_promoter)):(2+length(IDs_promoter)+length(IDs_body))] <- t(mmatrix[IDs_body,c(ANs,Ts)])
	colnames(temp) <- c("lib_size","read_count",paste("pr_",1:length(IDs_promoter),sep=""),paste("gb_",1:length(IDs_body),sep=""))
	rownames(temp) <- c(ANs,Ts)
	data_BRCA[[workingList_BRCA[i]]] <- temp
}