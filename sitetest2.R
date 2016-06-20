sitetest2 <- function (bmatrix, group, groupCol=2, gcase = "g2", gcontrol = "g1", testmethod = c("wilcox", "limma", "pooled", "satterthwaite"), Padj = "BH", concov = "OFF", rawpcut = NULL, adjustpcut = NULL, betadiffcut = NULL, paired = FALSE) 
{
    beta = bmatrix
    grouplev = group[, groupCol]
    if (concov == "ON") {
        cat("Performing linear regression....\n")
        require(MASS)
        testout = apply(beta, 1, function(x) {
            temp = summary(lm(x ~ as.numeric(as.character(grouplev))))
            pvalue = temp$coefficients[2, c(1, 4)]
            return(pvalue)
        })
        adjustP = p.adjust(testout[2, ], method = Padj)
        out = cbind(testout[2, ], adjustP, testout[1, ])
        rownames(out) = rownames(beta)
        colnames(out) = c("P-Value", "Adjust Pval", "Coefficient")
    }
    else {
        caseind = which(grouplev %in% gcase)
        controlind = which(grouplev %in% gcontrol)
        if (paired == TRUE) {
            lev1 = caseind[order(group[caseind, 3])]
            lev2 = controlind[order(group[controlind, 3])]
        }
        else {
            lev1 = caseind
            lev2 = controlind
        }
        eset = beta[, c(lev1, lev2)]
        if (testmethod == "wilcox") {
            cat("Performing Wilcox testing...\n")
            testout = apply(eset, 1, function(x) {
                wilcox.test(x[1:length(lev1)], x[(length(lev1) + 1):(length(lev1) + length(lev2))], paired = paired)$p.value
            })
        }
        if (testmethod == "limma") {
            require(limma)
            cat("Performing limma...\n")
            TS = as.factor(c(rep("T", length(lev1)), rep("C", length(lev2))))
            SS = as.factor(rep(1:length(lev1), 2))
            if (paired == FALSE) {
                design = model.matrix(~0 + TS)
                rownames(design) = colnames(eset)
                colnames(design) = c("C", "T")
                fit = lmFit(eset, design)
                cont.matrix = makeContrasts(comp = T - C, levels = design)
                fit2 = contrasts.fit(fit, cont.matrix)
                fit2 = eBayes(fit2)
                result1 = topTable(fit2, coef = 1, adjust.method = Padj, number = nrow(fit2))
            }
            else {
                design = model.matrix(~SS + TS)
                cat("Here is your design matrix\n")
                print(design)
                fit = lmFit(eset, design)
                fit2 = eBayes(fit)
                result1 = topTable(fit2, coef = "TST", adjust.method = Padj, number = nrow(fit2))
            }
            testout = result1[rownames(eset), "P.Value"]
        }
        if (testmethod == "pooled") {
            cat("Performing pooled t.test...\n")
            testout = apply(eset, 1, function(x) {
                t.test(x[1:length(lev1)], x[(length(lev1) + 1):(length(lev1) + length(lev2))], var.equal = TRUE, paired = paired)$p.value
            })
        }
        if (testmethod == "satterthwaite") {
            cat("Performing satterthwaite t.test...\n")
            testout = apply(eset, 1, function(x) {
                t.test(x[1:length(lev1)], x[(length(lev1) + 1):(length(lev1) + length(lev2))], paired = paired)$p.value
            })
        }
        adjustP = p.adjust(testout, method = Padj)
        difb = apply(eset, 1, function(x) {
            mean(x[1:length(lev1)], na.rm=T) - mean(x[(length(lev1) + 1):ncol(eset)], na.rm=T)
        })
        out = cbind(testout, adjustP, difb, rowMeans(eset[, 1:length(lev1)], na.rm=T), rowMeans(eset[, (length(lev1) + 1):ncol(eset)], na.rm=T))
        rownames(out) = rownames(eset)
        colnames(out) = c("P-Value", "Adjust Pval", "Beta-Difference", paste("Mean", paste(gcase, collapse = "_"), sep = "_"), paste("Mean", paste(gcontrol, collapse = "_"), sep = "_"))
    }
    out = outputDMfunc(out = out, rawpcut = rawpcut, adjustpcut = adjustpcut, betadiffcut = betadiffcut)
    return(out)
}
outputDMfunc <-
function(out,rawpcut = 0.05,adjustpcut =0.05,betadiffcut = 0.14){
     
          if(is.null(rawpcut)&is.null(adjustpcut)& is.null(betadiffcut)){
               cat("Kept the full comparison result. Please specify the significance criteria if you are only interested in the differentially methylated regions/sites only\n");
                out = out;
          }else{
                 if(!is.null(rawpcut)){
                        rawpcutout = out[,1]<=rawpcut
                  }else{
                        rawpcutout = out[,1]<=1
                  }
                  if(!is.null(adjustpcut)){
                        adjustpcutout = out[,2]<= adjustpcut
                  }else{
                       adjustpcutout = out[,2]<=1
                  }
                  if(!is.null(betadiffcut)){
                       betadiffcutout = abs(out[,3])>=betadiffcut;
                  }else{
                       betadiffcutout = abs(out[,3])>=0;
                  }
                  out = out[rawpcutout&adjustpcutout&betadiffcutout,]
           }
  }