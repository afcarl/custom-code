mmatrix_BRCA_BODY <- matrix(ncol=length(c(G1,G2)),nrow=17728,data=0)
colnames(mmatrix_BRCA_BODY) <- c(G1,G2)
row.names(mmatrix_BRCA_BODY) <- workingList_BRCA
mmatrix_BRCA_PROMOTER <- matrix(ncol=length(c(G1,G2)),nrow=17728,data=0)
colnames(mmatrix_BRCA_PROMOTER) <- c(G1,G2)
row.names(mmatrix_BRCA_PROMOTER) <- workingList_BRCA
for (i in 1:17728) {
	IDs_body <- unique(c(eval(parse(text = paste('GENEBODYInd$SID$','"',workingList_BRCA[i],'"',sep=""))), eval(parse(text = paste('UTR3Ind$SID$','"',workingList_BRCA[i],'"',sep="")))))
	temp = mmatrix_pc[IDs_body,c(G1,G2)]
	if (length(IDs_body)==1) {mmatrix_BRCA_BODY[i, ] = temp} else {mmatrix_BRCA_BODY[i, ] = apply(temp, 2, eval(mean), na.rm = TRUE)}
	IDs_promoter <- unique(c(eval(parse(text = paste('TSS1500Ind$SID$','"',workingList_BRCA[i],'"',sep=""))),eval(parse(text = paste('TSS200Ind$SID$','"',workingList_BRCA[i],'"',sep=""))), eval(parse(text = paste('UTR5Ind$SID$','"',workingList_BRCA[i],'"',sep="")))),eval(parse(text = paste('EXON1Ind$SID$','"',workingList_BRCA[i],'"',sep=""))))
	temp = mmatrix_pc[IDs_promoter,c(G1,G2)]
	if (length(IDs_promoter)==1) {mmatrix_BRCA_PROMOTER[i, ] = temp} else {mmatrix_BRCA_PROMOTER[i, ] = apply(temp, 2, eval(mean), na.rm = TRUE)}
}

mmatrix_BRCA_PROMOTER[which(mmatrix_BRCA_PROMOTER < -7)] <- -7
mmatrix_BRCA_PROMOTER[which(mmatrix_BRCA_PROMOTER > 7)] <- 7
mmatrix_BRCA_BODY[which(mmatrix_BRCA_BODY < -7)] <- -7
mmatrix_BRCA_BODY[which(mmatrix_BRCA_BODY > 7)] <- 7



mmatrix_BRCA_BODY <- matrix(ncol=length(c(G1,G2)),nrow=17728,data=0)
colnames(mmatrix_BRCA_BODY) <- c(G1,G2)
row.names(mmatrix_BRCA_BODY) <- workingList_BRCA
mmatrix_BRCA_PROMOTER <- matrix(ncol=length(c(G1,G2)),nrow=17728,data=0)
colnames(mmatrix_BRCA_PROMOTER) <- c(G1,G2)
row.names(mmatrix_BRCA_PROMOTER) <- workingList_BRCA
for (i in 1:17728) {
	IDs_body <- unique(c(eval(parse(text = paste('GENEBODYInd$SID$','"',workingList_BRCA[i],'"',sep=""))), eval(parse(text = paste('UTR3Ind$SID$','"',workingList_BRCA[i],'"',sep="")))))
	temp = mmatrix_pc[IDs_body,c(G1,G2)]
	if (length(IDs_body)==1) {mmatrix_BRCA_BODY[i, ] = temp} else {mmatrix_BRCA_BODY[i, ] = apply(temp, 2, eval(mean), na.rm = TRUE)}
	IDs_promoter <- unique(c(eval(parse(text = paste('TSS1500Ind$SID$','"',workingList_BRCA[i],'"',sep=""))),eval(parse(text = paste('TSS200Ind$SID$','"',workingList_BRCA[i],'"',sep=""))), eval(parse(text = paste('UTR5Ind$SID$','"',workingList_BRCA[i],'"',sep="")))),eval(parse(text = paste('EXON1Ind$SID$','"',workingList_BRCA[i],'"',sep=""))))
	temp = mmatrix_pc[IDs_promoter,c(G1,G2)]
	if (length(IDs_promoter)==1) {mmatrix_BRCA_PROMOTER[i, ] = temp} else {mmatrix_BRCA_PROMOTER[i, ] = apply(temp, 2, eval(mean), na.rm = TRUE)}
}
mmatrix_BRCA_BODY_noncap <- mmatrix_BRCA_BODY
mmatrix_BRCA_PROMOTER_noncap <- mmatrix_BRCA_PROMOTER

mmatrix_BRCA_PROMOTER[which(mmatrix_BRCA_PROMOTER < -7)] <- -7
mmatrix_BRCA_PROMOTER[which(mmatrix_BRCA_PROMOTER > 7)] <- 7
mmatrix_BRCA_BODY[which(mmatrix_BRCA_BODY < -7)] <- -7
mmatrix_BRCA_BODY[which(mmatrix_BRCA_BODY > 7)] <- 7




bmatrix_BRCA_BODY <- matrix(ncol=length(c(G1,G2)),nrow=17728,data=0)
colnames(bmatrix_BRCA_BODY) <- c(G1,G2)
row.names(bmatrix_BRCA_BODY) <- workingList_BRCA
bmatrix_BRCA_PROMOTER <- matrix(ncol=length(c(G1,G2)),nrow=17728,data=0)
colnames(bmatrix_BRCA_PROMOTER) <- c(G1,G2)
row.names(bmatrix_BRCA_PROMOTER) <- workingList_BRCA
for (i in 1:17728) {
	IDs_body <- unique(c(eval(parse(text = paste('GENEBODYInd$SID$','"',workingList_BRCA[i],'"',sep=""))), eval(parse(text = paste('UTR3Ind$SID$','"',workingList_BRCA[i],'"',sep="")))))
	temp = bmatrix_pc[IDs_body,c(G1,G2)]
	if (length(IDs_body)==1) {bmatrix_BRCA_BODY[i, ] = temp} else {bmatrix_BRCA_BODY[i, ] = apply(temp, 2, eval(mean), na.rm = TRUE)}
	IDs_promoter <- unique(c(eval(parse(text = paste('TSS1500Ind$SID$','"',workingList_BRCA[i],'"',sep=""))),eval(parse(text = paste('TSS200Ind$SID$','"',workingList_BRCA[i],'"',sep=""))), eval(parse(text = paste('UTR5Ind$SID$','"',workingList_BRCA[i],'"',sep="")))),eval(parse(text = paste('EXON1Ind$SID$','"',workingList_BRCA[i],'"',sep=""))))
	temp = bmatrix_pc[IDs_promoter,c(G1,G2)]
	if (length(IDs_promoter)==1) {bmatrix_BRCA_PROMOTER[i, ] = temp} else {bmatrix_BRCA_PROMOTER[i, ] = apply(temp, 2, eval(mean), na.rm = TRUE)}
}