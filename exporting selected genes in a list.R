data_BRCA <- list()
for (i in 1:17728) {
	IDs_promoter <- unique(c(eval(parse(text = paste('TSS1500Ind$SID$','"',workingList[i],'"',sep=""))),eval(parse(text = paste('TSS200Ind$SID$','"',workingList[i],'"',sep=""))), eval(parse(text = paste('UTR5Ind$SID$','"',workingList[i],'"',sep="")))),eval(parse(text = paste('EXON1Ind$SID$','"',workingList[i],'"',sep=""))))
	IDs_body <- unique(c(eval(parse(text = paste('GENEBODYInd$SID$','"',workingList[i],'"',sep=""))), eval(parse(text = paste('UTR3Ind$SID$','"',workingList[i],'"',sep="")))))
	ncol = length(IDs_promoter) + length(IDs_body) + 2
	nrow = length(c(G1,G2))
	temp <- matrix(nrow=nrow,ncol=ncol)
	temp[,1] <- unlist(factors_ls[c(G1,G2)])*10^6
	temp[,2] <- unlist(counts[workingList[i],c(G1,G2)])
	temp[,3:(2+length(IDs_promoter))] <- t(mmatrix_pc[IDs_promoter,c(G1,G2)])
	temp[,(3+length(IDs_promoter)):(2+length(IDs_promoter)+length(IDs_body))] <- t(mmatrix_pc[IDs_body,c(G1,G2)])
	colnames(temp) <- c("lib_size","read_count",paste("pr_",1:length(IDs_promoter),sep=""),paste("gb_",1:length(IDs_body),sep=""))
	rownames(temp) <- c(G1,G2)
	data_BRCA[[workingList[i]]] <- temp
}

data_BRCA_progressing <- list()
for (i in 1:17728) {
	IDs_promoter <- unique(c(eval(parse(text = paste('TSS1500Ind$SID$','"',workingList[i],'"',sep=""))),eval(parse(text = paste('TSS200Ind$SID$','"',workingList[i],'"',sep=""))), eval(parse(text = paste('UTR5Ind$SID$','"',workingList[i],'"',sep="")))),eval(parse(text = paste('EXON1Ind$SID$','"',workingList[i],'"',sep=""))))
	IDs_body <- unique(c(eval(parse(text = paste('GENEBODYInd$SID$','"',workingList[i],'"',sep=""))), eval(parse(text = paste('UTR3Ind$SID$','"',workingList[i],'"',sep="")))))
	ncol = length(IDs_promoter) + length(IDs_body) + 2
	nrow = length(c(progressed,nonProgressed))
	temp <- matrix(nrow=nrow,ncol=ncol)
	temp[,1] <- true_ls[c(progressed,nonProgressed)]
	temp[,2] <- unlist(counts[workingList[i],c(progressed,nonProgressed)])
	temp[,3:(2+length(IDs_promoter))] <- t(mmatrix_pc[IDs_promoter,c(progressed,nonProgressed)])
	temp[,(3+length(IDs_promoter)):(2+length(IDs_promoter)+length(IDs_body))] <- t(mmatrix_pc[IDs_body,c(progressed,nonProgressed)])
	colnames(temp) <- c("lib_size","read_count",paste("pr_",1:length(IDs_promoter),sep=""),paste("gb_",1:length(IDs_body),sep=""))
	rownames(temp) <- c(progressed,nonProgressed)
	data_BRCA_progressing[[workingList[i]]] <- temp
}