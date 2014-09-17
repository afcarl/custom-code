promoter_ncpgs <- vector(length(workingList_BRCA),mode="numeric")
body_ncpgs <- vector(length(workingList_BRCA),mode="numeric")
for (i in 1:length(workingList_BRCA)) {
	promoter_ncpgs[i] <- length(unique(c(eval(parse(text = paste('TSS1500Ind$SID$','"',workingList_BRCA[i],'"',sep=""))),eval(parse(text = paste('TSS200Ind$SID$','"',workingList_BRCA[i],'"',sep=""))), eval(parse(text = paste('UTR5Ind$SID$','"',workingList_BRCA[i],'"',sep="")))),eval(parse(text = paste('EXON1Ind$SID$','"',workingList_BRCA[i],'"',sep="")))))
	body_ncpgs[i] <- length(unique(c(eval(parse(text = paste('GENEBODYInd$SID$','"',workingList_BRCA[i],'"',sep=""))), eval(parse(text = paste('UTR3Ind$SID$','"',workingList_BRCA[i],'"',sep=""))))))
}