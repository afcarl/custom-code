for (i in 1:length(workingList_BRCA4)){
	IDs_promoter <- unique(c(eval(parse(text = paste('TSS1500Ind$SID$','"',workingList_BRCA4[i],'"',sep=""))),eval(parse(text = paste('TSS200Ind$SID$','"',workingList_BRCA4[i],'"',sep=""))), eval(parse(text = paste('UTR5Ind$SID$','"',workingList_BRCA4[i],'"',sep="")))),eval(parse(text = paste('EXON1Ind$SID$','"',workingList_BRCA4[i],'"',sep=""))))
	IDs_body <- unique(c(eval(parse(text = paste('GENEBODYInd$SID$','"',workingList_BRCA4[i],'"',sep=""))), eval(parse(text = paste('UTR3Ind$SID$','"',workingList_BRCA4[i],'"',sep="")))))
	for (current_cample in 1:150){
		bw = sd(bmatrix[IDs_body,c(Ts2,ANs2)[current_sample]])
		if (!is.finite(bw)) print(paste("non-finite body 'bw'",i,current_cample,sep=" "))
		bw = sd(bmatrix[IDs_promoter,c(Ts2,ANs2)[current_sample]])
		if (!is.finite(bw)) print(paste("non-finite promoter 'bw'",i,current_cample,sep=" "))
	}
}