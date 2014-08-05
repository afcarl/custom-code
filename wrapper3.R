function (dataf, list11Rdata)
{
    beta = dataf@bmatrix
    groupinfo = dataf@groupinfo
    list11 = c("TSS1500Ind", "GENEBODYInd", "ISLANDInd")
    
    gcase="T"
    gcontrol=c("N","AN")
    paired=FALSE
    testmethod="limma"
    Padj="BH"
    concov = "OFF"
	rawpcut = 0.01
	adjustpcut = 0.01
    annot.df = as.data.frame(dataf@annot)
    
    indexTSS1500 = which(grepl("TSS1500",annot.df[,"UCSC_REFGENE_GROUP"]))
    indexGENEBODY = which(grepl("Body",annot.df[,"UCSC_REFGENE_GROUP"]))
    indexISLAND = which(grepl("Island",annot.df[,"RELATION_TO_UCSC_CPG_ISLAND"]))
    
    indexISLANDu = indexISLAND[!indexISLAND[] %in% indexGENEBODY[]]
    indexISLANDu = indexISLANDu[!indexISLANDu[] %in% indexTSS1500[]]
    
    indexGENEBODYu = indexGENEBODY[!indexGENEBODY[] %in% indexTSS1500[]]
    indexGENEBODYu = indexGENEBODYu[!indexGENEBODYu[] %in% indexISLAND[]]
    
    indexTSS1500u = indexTSS1500[!indexTSS1500[] %in% indexGENEBODY[]]
    indexTSS1500u = indexTSS1500u[!indexTSS1500u[] %in% indexISLAND[]]
    
    
    indexlistTSS1500 = dataf@TSS1500Ind
    
    eset = indexregionfunc2(indexlistTSS1500, indexTSS1500u, beta, "mean")
    TSS1500Indtest = as.data.frame(testfunc2(eset,concov=concov,testmethod=testmethod,Padj=Padj,groupinfo = groupinfo,paired = paired,gcase = gcase,gcontrol = gcontrol))
    TSS1500Indtest=TSS1500Indtest[!is.na(TSS1500Indtest[,3]),]
	TSS1500Indtest=outputDMfunc(TSS1500Indtest,rawpcut = rawpcut,adjustpcut = adjustpcut,betadiffcut = betadiffcut)
	
	
    indexlistGENEBODY = dataf@GENEBODYInd
    
    eset = indexregionfunc2(indexlistGENEBODY,indexGENEBODYu, beta, "mean")
    GENEBODYIndtest = as.data.frame(testfunc2(eset,concov=concov,testmethod=testmethod,Padj=Padj,groupinfo = groupinfo,paired = paired,gcase = gcase,gcontrol = gcontrol))
    GENEBODYIndtest=GENEBODYIndtest[!is.na(GENEBODYIndtest[,3]),]
	GENEBODYIndtest=outputDMfunc(GENEBODYIndtest,rawpcut = rawpcut,adjustpcut = adjustpcut,betadiffcut = betadiffcut)
    
    
    indexlistISLAND = dataf@ISLANDInd
    
    eset = indexregionfunc2(indexlistISLAND,indexISLANDu, beta, "mean")
    ISLANDIndtest = as.data.frame(testfunc2(eset,concov=concov,testmethod=testmethod,Padj=Padj,groupinfo = groupinfo,paired = paired,gcase = gcase,gcontrol = gcontrol))
    ISLANDIndtest=ISLANDIndtest[!is.na(ISLANDIndtest[,3]),]    
    ISLANDIndtest=outputDMfunc(ISLANDIndtest,rawpcut = rawpcut,adjustpcut = adjustpcut,betadiffcut = betadiffcut)
	
    
    eval(parse(text = paste("save(", paste(paste(list11, "test", 
                                                 sep = ""), collapse = ","), ",file = list11Rdata)", sep = "")))
    
    
}
