IMA.methy450PP2 <- function (data, na.omit = FALSE, peakcorrection = FALSE, normalization = FALSE, transfm = c(FALSE, "arcsinsqr", "logit"), samplefilterdetectP = c(FALSE, 1e-05), samplefilterperc = 0.75, sitefilterdetectP = c(FALSE, 0.05), sitefilterperc = 0.75, locidiff = c(FALSE, 0.01), locidiffgroup = list("g1", "g2"), XYchrom = c(FALSE, "X", "Y", c("X", "Y")), snpfilter = c(FALSE, "snpsites.txt")) 
{
    bmatrix = data@bmatrix
    detect_p = data@detectP
    annotation = data@annot
    groupinfo = data@groupinfo
    orignalrownm = rownames(bmatrix)
    if (samplefilterdetectP) {
        goodsample = colSums(detect_p <= samplefilterdetectP) >= 
            samplefilterperc * nrow(detect_p)
        bmatrix = bmatrix[, goodsample]
        detect_p = detect_p[, goodsample]
        cat(ncol(bmatrix) - length(goodsample), "samples removed with at least", 
            samplefilterperc * 100, "percentage sites having pvalue greater than", 
            samplefilterdetectP, "\n")
        groupinfo = groupinfo[goodsample, ]
    }
    if (snpfilter != FALSE) {
        snpsites = read.delim(snpfilter, sep = "\t", stringsAsFactors = FALSE)[, 
            "TargetID"]
        index = rownames(bmatrix) %in% snpsites
        bmatrix = bmatrix[!index, ]
        detect_p = detect_p[!index, ]
        annotation = annotation[!index, ]
        cat(sum(index), "sites contain snps and removed,\n")
    }
    if (XYchrom[1] != FALSE) {
        chr = annotation[, "CHR"]
        index = which(chr %in% XYchrom)
        good_chrom = rownames(annotation)[-index]
        cat(length(index), "sites on chr", XYchrom, "are removed\n")
    }
    else {
        good_chrom = rownames(bmatrix)
    }
    if (sitefilterdetectP) {
        good_loci = rownames(detect_p)[rowSums(detect_p <= sitefilterdetectP) >= 
            sitefilterperc * ncol(detect_p)]
        cat(nrow(detect_p) - length(good_loci), "sites had at least", 
            sitefilterperc * 100, "% samples with pvalue great than", 
            sitefilterdetectP, "and are removed\n")
    }
    else {
        good_loci = rownames(bmatrix)
    }
    if (locidiff) {
        c1 = groupinfo[, 2] %in% locidiffgroup[[1]]
        c2 = groupinfo[, 2] %in% locidiffgroup[[2]]
        con_mean = apply(bmatrix[, c1], 1, mean)
        trt_mean = apply(bmatrix[, c2], 1, mean)
        good_diff = rownames(bmatrix)[abs(trt_mean - con_mean) >= 
            locidiff]
        cat(length(good_diff), "sites had the beta difference between group great than", 
            locidiff, "and are kept for the downstream analysis \n")
    }
    else {
        good_diff = rownames(bmatrix)
    }
    all_good = intersect(intersect(good_chrom, good_loci), good_diff)
    cat(length(all_good), "sites were retained from the original", 
        length(orignalrownm), "sites\n")
    bmatrix = bmatrix[all_good, ]
    annotation = annotation[all_good, ]
    detect_p = detect_p[all_good, ]
    if (peakcorrection) {
        cat("Peak correction...\nThis part of code was provided by Matthieu Defrance <defrance@bigre.ulb.ac.be>\n")
        cat("Thanks for sharing the code with us.\n")
        cat("Dimension of beta matrix", dim(bmatrix), "\n")
        cat("Dimension of annotation", dim(annotation), "\n")
        bmatrix = peak.correction(bmatrix, annotation)
        bmatrix = bmatrix[rownames(annotation), ]
    }
    if (normalization) {
        require(preprocessCore)
        bmatrix = normalize.quantiles(as.matrix(bmatrix))
        colnames(bmatrix) = colnames(detect_p)
        rownames(bmatrix) = rownames(detect_p)
        cat("Quantile normalization Performed\n")
    }
    if (transfm == "arcsinsqr") {
        if (na.omit) {
            bmatrix = asin(sqrt(bmatrix))
            cat("Transfer beta matrix by the arcsin square root\n")
        }
        else {
            cat("\tMissing value exist in the orignial data,\nPlease remove the missing value before transformation,use na.omit = TRUE\n")
            stop
        }
    }
    if (transfm == "logit") {
		bmatrix[bmatrix == 0] <- min(bmatrix[bmatrix > 0 & !is.na(bmatrix_pc)], 0.001)/10
		bmatrix[bmatrix == 1] <- max(bmatrix[bmatrix < 1 & !is.na(bmatrix_pc)], 0.999) + (1 - max(bmatrix[bmatrix < 1], 0.999))/100
		bmatrix = log2(bmatrix/(1 - bmatrix))
		cat("Transfer beta matrix by the logit transformation \n")
    }
    cat(".......\nSplit the annotation file to 11 annotated region categories\n.......\n\n")
    annot = annotation
    name = "UCSC_REFGENE_NAME"
    cpGsite = as.character(annot[, 1])
    genelist = strsplit(as.character(annot[, name]), ";")
    genelist[which(genelist == "character(0)")] = "NA"
    name = "UCSC_REFGENE_GROUP"
    refgene = strsplit(as.character(annot[, name]), ";")
    refgene[which(refgene == "character(0)")] = "NA"
    listlength = lapply(refgene, length)
    listlength[listlength == 0] = 1
    col0 = rep(1:nrow(annot), listlength)
    col1 = rep(cpGsite, listlength)
    col2 = unlist(genelist)
    col3 = unlist(refgene)
    col4 = rep(as.character(annotation[, "RELATION_TO_UCSC_CPG_ISLAND"]), 
        listlength)
    col5 = rep(as.character(annotation[, "UCSC_CPG_ISLANDS_NAME"]), 
        listlength)
    splitToRegionlist = function(grepname = c("TSS1500", "TSS200", 
        "5'UTR", "1stExon", "Gene Body", "3'UTR")) {
        index = col3 == grepname
        col1sub = col1[index]
        col2sub = col2[index]
        temp = split(col1sub, col2sub)
        returnSID = lapply(temp, unique)
        col0sub = col0[index]
        temp = split(col0sub, col2sub)
        returnPID = lapply(temp, unique)
        return(Ind = list(SID = returnSID, PID = returnPID))
    }
    TSS1500Ind = splitToRegionlist(grepname = "TSS1500")
    TSS200Ind = splitToRegionlist(grepname = "TSS200")
    UTR5Ind = splitToRegionlist(grepname = "5'UTR")
    EXON1Ind = splitToRegionlist(grepname = "1stExon")
    GENEBODYInd = splitToRegionlist(grepname = "Body")
    UTR3Ind = splitToRegionlist(grepname = "3'UTR")
    cat("TSS1500 region contains:", length(TSS1500Ind$SID), 
        "UCSC REFGENE region \nTSS200 region contains:", length(TSS200Ind$SID), 
        "UCSC REFGENE region\n5'UTR region contains:", length(UTR5Ind$SID), 
        "UCSC REFGENE region\n1st Exon region contains:", length(EXON1Ind$SID), 
        "UCSC REFGENE region\nGene body region contains:", length(GENEBODYInd$SID), 
        "UCSC REFGENE region\n3'UTR region contains:", length(UTR3Ind$SID), 
        "UCSC REFGENE region\n")
    splitToRegionlist2 = function(grepname = c("Island", "N_Shore", 
        "S_Shore", "N_Shelf", "S_Shelf")) {
        index = col4 == grepname
        col1sub = col1[index]
        col5sub = col5[index]
        temp = split(col1sub, col5sub)
        returnSID = lapply(temp, unique)
        col0sub = col0[index]
        temp = split(col0sub, col5sub)
        returnPID = lapply(temp, unique)
        return(Ind = list(SID = returnSID, PID = returnPID))
    }
    ISLANDInd = splitToRegionlist2(grepname = "Island")
    NSHOREInd = splitToRegionlist2(grepname = "N_Shore")
    SSHOREInd = splitToRegionlist2(grepname = "S_Shore")
    NSHELFInd = splitToRegionlist2(grepname = "N_Shelf")
    SSHELFInd = splitToRegionlist2(grepname = "S_Shelf")
    cat("Island region contains:", length(ISLANDInd$SID), "UCSC CPG ISLAND region\nN_Shore region contains", 
        length(NSHOREInd$SID), "UCSC CPG ISLAND region\nS_Shore region contains", 
        length(SSHOREInd$SID), "UCSC CPG ISLAND region\nN_Shelf region contains", 
        length(NSHELFInd$SID), "UCSC CPG ISLAND region\nS_Shelf region contains", 
        length(SSHELFInd$SID), "UCSC CPG ISLAND region\n")
    setClass("methy450batch", representation(bmatrix = "matrix", 
        annot = "matrix", detectP = "matrix", groupinfo = "data.frame", 
        TSS1500Ind = "list", TSS200Ind = "list", UTR5Ind = "list", 
        EXON1Ind = "list", GENEBODYInd = "list", UTR3Ind = "list", 
        ISLANDInd = "list", NSHOREInd = "list", SSHOREInd = "list", 
        NSHELFInd = "list", SSHELFInd = "list"), where = topenv(parent.frame()))
    x.methy450 = new("methy450batch", bmatrix = as.matrix(bmatrix), 
        annot = as.matrix(annotation), detectP = as.matrix(detect_p), 
        groupinfo = groupinfo, TSS1500Ind = TSS1500Ind, TSS200Ind = TSS200Ind, 
        UTR5Ind = UTR5Ind, EXON1Ind = EXON1Ind, GENEBODYInd = GENEBODYInd, 
        UTR3Ind = UTR3Ind, ISLANDInd = ISLANDInd, NSHOREInd = NSHOREInd, 
        SSHOREInd = SSHOREInd, NSHELFInd = NSHELFInd, SSHELFInd = SSHELFInd)
    cat("\nA methy450batch class is created and the slotNames are:\n", 
        slotNames(x.methy450), "\n")
    return(x.methy450)
}
