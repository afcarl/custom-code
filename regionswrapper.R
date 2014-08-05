function (dataf, indexmethod = c("mean", "median", "tbrm"), gcase = "g2", 
    gcontrol = "g1", testmethod = c("wilcox", "limma", "pooled", 
        "satterthwaite"), Padj = "BH", concov = c("OFF", "ON"), 
    paired = FALSE, list11excel, list11Rdata, rawpcut = NULL, 
    adjustpcut = NULL, betadiffcut = NULL) 
{
    beta = dataf@bmatrix
    groupinfo = dataf@groupinfo
    indexlist = NULL
    list11 = c("TSS1500Ind", "TSS200Ind", "UTR5Ind", "EXON1Ind", 
        "GENEBODYInd", "UTR3Ind", "ISLANDInd", "NSHELFInd", "NSHOREInd", 
        "SSHELFInd", "SSHOREInd")
    for (i in 1:11) {
        cat("calculating", list11[i], "\n")
        eval(parse(text = paste("indexlist=dataf@", list11[i], 
            sep = "")))
        eset = indexregionfunc(indexlist, beta, indexmethod)
        cat(list11[i], ": Collected the loci for each region and derived an index of overall region methylation value\nStart testing\n")
        eval(parse(text = paste(list11[i], "test=as.data.frame(testfunc(eset,concov = concov,testmethod=testmethod,Padj=Padj,groupinfo = groupinfo,paired = paired,gcase = gcase,gcontrol = gcontrol))", 
            sep = "")))
        eval(parse(text = paste(list11[i], "test=outputDMfunc(", 
            list11[i], "test,rawpcut = rawpcut,adjustpcut = adjustpcut,betadiffcut = betadiffcut)", 
            sep = "")))
    }
    eval(parse(text = paste("save(", paste(paste(list11, "test", 
        sep = ""), collapse = ","), ",file = list11Rdata)", sep = "")))
    require(WriteXLS)
    WriteXLS(paste(list11, "test", sep = ""), ExcelFileName = list11excel, 
        SheetNames = list11, row.names = TRUE)
}
