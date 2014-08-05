function (dataf, list11Rdata)
{
    beta = dataf@bmatrix
    groupinfo = dataf@groupinfo
    list11 = c("TSS1500Ind", "TSS200Ind", "UTR5Ind", "EXON1Ind", 
               "GENEBODYInd", "UTR3Ind", "ISLANDInd", "NSHELFInd", "NSHOREInd", 
               "SSHELFInd", "SSHOREInd")
    gcase="T"
    gcontrol=c("N","AN")
    paired=FALSE
    testmethod="limma"
    Padj="BH"
    concov = "OFF"
    annot.df = as.data.frame(dataf@annot)
    
    index = row.names(annot.df[grepl("TSS1500",annot.df[,"UCSC_REFGENE_GROUP"]),])
    betasubset = dataf@bmatrix[index,]
    TSS1500Indtest = as.data.frame(testfunc2(eset=betasubset,concov=concov,testmethod=testmethod,Padj=Padj,groupinfo = groupinfo,paired = paired,gcase = gcase,gcontrol = gcontrol))
    
    index = row.names(annot.df[grepl("TSS200",annot.df[,"UCSC_REFGENE_GROUP"]),])
    betasubset = dataf@bmatrix[index,]
    TSS200Indtest = as.data.frame(testfunc2(eset=betasubset,concov=concov,testmethod=testmethod,Padj=Padj,groupinfo = groupinfo,paired = paired,gcase = gcase,gcontrol = gcontrol))
    
    index = row.names(annot.df[grepl("3'UTR",annot.df[,"UCSC_REFGENE_GROUP"]),])
    betasubset = dataf@bmatrix[index,]
    UTR3Indtest = as.data.frame(testfunc2(eset=betasubset,concov=concov,testmethod=testmethod,Padj=Padj,groupinfo = groupinfo,paired = paired,gcase = gcase,gcontrol = gcontrol))
    
    index = row.names(annot.df[grepl("1stExon",annot.df[,"UCSC_REFGENE_GROUP"]),])
    betasubset = dataf@bmatrix[index,]
    EXON1Indtest = as.data.frame(testfunc2(eset=betasubset,concov=concov,testmethod=testmethod,Padj=Padj,groupinfo = groupinfo,paired = paired,gcase = gcase,gcontrol = gcontrol))
    
    index = row.names(annot.df[grepl("Body",annot.df[,"UCSC_REFGENE_GROUP"]),])
    betasubset = dataf@bmatrix[index,]
    GENEBODYIndtest = as.data.frame(testfunc2(eset=betasubset,concov=concov,testmethod=testmethod,Padj=Padj,groupinfo = groupinfo,paired = paired,gcase = gcase,gcontrol = gcontrol))
    
    index = row.names(annot.df[grepl("5'UTR",annot.df[,"UCSC_REFGENE_GROUP"]),])
    betasubset = dataf@bmatrix[index,]
    UTR5Indtest = as.data.frame(testfunc2(eset=betasubset,concov=concov,testmethod=testmethod,Padj=Padj,groupinfo = groupinfo,paired = paired,gcase = gcase,gcontrol = gcontrol))
    
    index = dataf@annot[,"RELATION_TO_UCSC_CPG_ISLAND"]%in%"N_Shelf"
    betasubset = dataf@bmatrix[index,]
    NSHELFIndtest = as.data.frame(testfunc2(eset=betasubset,concov=concov,testmethod=testmethod,Padj=Padj,groupinfo = groupinfo,paired = paired,gcase = gcase,gcontrol = gcontrol))
    
    index = dataf@annot[,"RELATION_TO_UCSC_CPG_ISLAND"]%in%"N_Shore"
    betasubset = dataf@bmatrix[index,]
    NSHOREIndtest = as.data.frame(testfunc2(eset=betasubset,concov=concov,testmethod=testmethod,Padj=Padj,groupinfo = groupinfo,paired = paired,gcase = gcase,gcontrol = gcontrol))
    
    index = dataf@annot[,"RELATION_TO_UCSC_CPG_ISLAND"]%in%"Island"
    betasubset = dataf@bmatrix[index,]
    ISLANDIndtest = as.data.frame(testfunc2(eset=betasubset,concov=concov,testmethod=testmethod,Padj=Padj,groupinfo = groupinfo,paired = paired,gcase = gcase,gcontrol = gcontrol))
    
    index = dataf@annot[,"RELATION_TO_UCSC_CPG_ISLAND"]%in%"S_Shore"
    betasubset = dataf@bmatrix[index,]
    SSHOREIndtest = as.data.frame(testfunc2(eset=betasubset,concov=concov,testmethod=testmethod,Padj=Padj,groupinfo = groupinfo,paired = paired,gcase = gcase,gcontrol = gcontrol))
    
    index = dataf@annot[,"RELATION_TO_UCSC_CPG_ISLAND"]%in%"S_Shelf"
    betasubset = dataf@bmatrix[index,]
    SSHELFIndtest = as.data.frame(testfunc2(eset=betasubset,concov=concov,testmethod=testmethod,Padj=Padj,groupinfo = groupinfo,paired = paired,gcase = gcase,gcontrol = gcontrol))
    
    
    
    eval(parse(text = paste("save(", paste(paste(list11, "test", 
                                                 sep = ""), collapse = ","), ",file = list11Rdata)", sep = "")))
	
    
}
