data_BRCA_all <- list()
for (i in 1:17728) {
	IDs_promoter <- unique(c(eval(parse(text = paste('TSS1500Ind$SID$','"',workingList[i],'"',sep=""))),eval(parse(text = paste('TSS200Ind$SID$','"',workingList[i],'"',sep=""))), eval(parse(text = paste('UTR5Ind$SID$','"',workingList[i],'"',sep="")))),eval(parse(text = paste('EXON1Ind$SID$','"',workingList[i],'"',sep=""))))
	IDs_body <- unique(c(eval(parse(text = paste('GENEBODYInd$SID$','"',workingList[i],'"',sep=""))), eval(parse(text = paste('UTR3Ind$SID$','"',workingList[i],'"',sep="")))))
	ncol = length(IDs_promoter) + length(IDs_body) + 2
	nrow = length(c(ANs,Ts))
	temp <- matrix(nrow=nrow,ncol=ncol)
	temp[,1] <- unlist(true_ls[c(ANs,Ts)])
	temp[,2] <- unlist(counts[workingList[i],c(ANs,Ts)])
	temp[,3:(2+length(IDs_promoter))] <- t(mmatrix_pc[IDs_promoter,c(ANs,Ts)])
	temp[,(3+length(IDs_promoter)):(2+length(IDs_promoter)+length(IDs_body))] <- t(mmatrix_pc[IDs_body,c(ANs,Ts)])
	colnames(temp) <- c("lib_size","read_count",paste("pr_",1:length(IDs_promoter),sep=""),paste("gb_",1:length(IDs_body),sep=""))
	rownames(temp) <- c(ANs,Ts)
	data_BRCA_all[[workingList[i]]] <- temp
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


# Find and order all CpG sites according to genomic location
IDs_promoter <- NULL
IDs_body <- NULL
for (i in 1:length(GB_IDs)) {
  IDs_promoter[[i]] <- unique(c(TSS1500Ind$SID[[Pr_IDs[i], exact=TRUE]], TSS200Ind$SID[[Pr_IDs[i], exact=TRUE]], UTR5Ind$SID[[Pr_IDs[i], exact=TRUE]], EXON1Ind$SID[[Pr_IDs[i], exact=TRUE]]))
	IDs_promoter[[i]] <- IDs_promoter[[i]][order(fullannot[IDs_promoter[[i]], "pos"])]
	
	IDs_body[[i]] <- unique(c(GENEBODYInd$SID[[GB_IDs[i], exact=TRUE]], UTR3Ind$SID[[GB_IDs[i], exact=TRUE]]))
	IDs_body[[i]] <- IDs_body[[i]][order(fullannot[IDs_body[[i]], "pos"])]
}

for (i in length(GB_IDs):length(Pr_IDs)) {
	IDs_promoter[[i]] <- unique(c(TSS1500Ind$SID[[Pr_IDs[i], exact=TRUE]], TSS200Ind$SID[[Pr_IDs[i], exact=TRUE]], UTR5Ind$SID[[Pr_IDs[i], exact=TRUE]], EXON1Ind$SID[[Pr_IDs[i], exact=TRUE]]))
	IDs_promoter[[i]] <- IDs_promoter[[i]][order(fullannot[IDs_promoter[[i]], "pos"])]
}
names(IDs_promoter) <- Pr_IDs
names(IDs_body) <- GB_IDs

workingList <- intersect(intersect(rownames(counts), names(IDs_body)), names(IDs_promoter))
data_BRCA_all <- list()
for (i in 1:length(workingList)) {
    current_ID <- workingList[i]
    ncol = length(IDs_promoter[[current_ID]]) + length(IDs_body[[current_ID]]) + 2
    nrow = length(c(ANs,Ts))
    temp <- matrix(nrow=nrow,ncol=ncol)
    temp[,1] <- unlist(true_ls[c(ANs,Ts)])
    temp[,2] <- unlist(counts[current_ID,c(ANs,Ts)])
    temp[,3:(2+length(IDs_promoter[[current_ID]]))] <- t(mmatrix_pc[IDs_promoter[[current_ID]],c(ANs,Ts)])
    temp[,(3+length(IDs_promoter[[current_ID]])):(2+length(IDs_promoter[[current_ID]])+length(IDs_body[[current_ID]]))] <- t(mmatrix_pc[IDs_body[[current_ID]],c(ANs,Ts)])
    colnames(temp) <- c("lib_size","read_count",paste("pr_",1:length(IDs_promoter[[current_ID]]),sep=""),paste("gb_",1:length(IDs_body[[current_ID]]),sep=""))
    rownames(temp) <- c(ANs,Ts)
    data_BRCA_all[[current_ID]] <- temp
}
save(data_BRCA_all, ANs, Ts, file="data_BRCA_all.RData")

data_BRCA_progressing <- list()
for (i in 1:length(workingList)) {
    current_ID <- workingList[i]
    ncol = length(IDs_promoter[[current_ID]]) + length(IDs_body[[current_ID]]) + 2
    nrow = length(c(progressed,nonProgressed))
    temp <- matrix(nrow=nrow,ncol=ncol)
    temp[,1] <- unlist(true_ls[c(progressed,nonProgressed)])
    temp[,2] <- unlist(counts[current_ID,c(progressed,nonProgressed)])
    temp[,3:(2+length(IDs_promoter[[current_ID]]))] <- t(mmatrix_pc[IDs_promoter[[current_ID]],c(progressed,nonProgressed)])
    temp[,(3+length(IDs_promoter[[current_ID]])):(2+length(IDs_promoter[[current_ID]])+length(IDs_body[[current_ID]]))] <- t(mmatrix_pc[IDs_body[[current_ID]],c(progressed,nonProgressed)])
    colnames(temp) <- c("lib_size","read_count",paste("pr_",1:length(IDs_promoter[[current_ID]]),sep=""),paste("gb_",1:length(IDs_body[[current_ID]]),sep=""))
    rownames(temp) <- c(progressed,nonProgressed)
    data_BRCA_all[[current_ID]] <- temp
}
save(data_BRCA_progressing, progressed, nonProgressed, file="data_BRCA_progressing.RData")