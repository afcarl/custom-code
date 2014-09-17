load("essentials_SIM.RData")
library(edgeR)
library(IMA)

fishersMethod <- function(x) pchisq(-2 * sum(log(x)),df=2*length(x),lower=FALSE)

ids <- c("001","002","003","004","005","006","007","008","009","010","011","012","013","014","015","016","017","018","019","020","021","022","023","024","025","026","027","028","029","030","031","032","033","034","035","036","037","038","039","040","041","042","043","044","045","046","047","048","049","050","051","052","053","054","055","056","057","058","059","060","061","062","063","064","065","066","067","068","069","070","071","072","073","074","075","076","077","078","079","080","081","082","083","084","085","086","087","088","089","090","091","092","093","094","095","096","097","098","099","100")
ANs <- paste("AN_",ids,sep="")
Ts <- paste("TU_",ids,sep="")

# expression
expr_data <- matrix(ncol=200,nrow=length(sim.data))
colnames(expr_data) <- c(ANs,Ts)
rownames(expr_data) <- paste("id_",seq(1,nrow(expr_data)),sep="")
for (i in 1:nrow(expr_data)) expr_data[i,] <- c(sim.data[[i]][[1]][,1],sim.data[[i]][[2]][,1])

group <- factor(c(rep("AN",100),rep("T",100)))
y <- DGEList(counts=expr_data,group=group)
y$samples[,2] <- 10^7
y$samples[,3] <- 1
y <- estimateCommonDisp(y)
y <- estimateTagwiseDisp(y)
et <- exactTest(y,dispersion="tagwise")

Ts <- paste("TU_",ids,sep="")
ANs <- paste("AN_",ids,sep="")
groupinfo <- as.data.frame(rbind(cbind(ANs,"AN"),cbind(Ts,"T")))
rownames(groupinfo) <- groupinfo[,1]
colnames(groupinfo) <- c("Sample.name","Group")

gb_meth <- matrix(ncol=200,nrow=length(sim.data))
for (i in 1:nrow(gb_meth)) gb_meth[i,] <- c(apply(sim.data[[i]][[1]][,4:14],1,mean),apply(sim.data[[i]][[2]][,4:14],1,mean))
rownames(gb_meth) <- paste("id_",seq(1,nrow(gb_meth)),sep="")
colnames(gb_meth) <- c(ANs,Ts)
gb_res <- sitetest2(gb_meth,groupinfo,gcase="T",gcontrol = "AN",testmethod = "satterthwaite")

pr_meth <- matrix(ncol=200,nrow=length(sim.data))
for (i in 1:nrow(pr_meth)) pr_meth[i,] <- c(apply(sim.data[[i]][[1]][,15:25],1,mean),apply(sim.data[[i]][[2]][,15:25],1,mean))
rownames(pr_meth) <- paste("id_",seq(1,nrow(pr_meth)),sep="")
colnames(pr_meth) <- c(ANs,Ts)
pr_res <- sitetest2(pr_meth,groupinfo,gcase="T",gcontrol = "AN",testmethod = "satterthwaite")

results <- cbind(gb_res[,1],pr_res[,1],et$table[,3])
results <- as.data.frame(results)
results$fishers <- apply(results,1,fishersMethod)

library(pROC)
calculateAUC <- function(x,vector) {auc(predictor=c(vector[(1+x*1000):(1000+x*1000)],vector[(21001+x*1000):(22000+x*1000)]),response=c(rep("Pos",1000),rep("Neg",1000)))}
calculateAUC <- function(x,vector) {auc(predictor=c(vector[(1+x*1000):(1000+x*1000)],vector[(11001+x*1000):(12000+x*1000)]),response=c(rep("Pos",1000),rep("Neg",1000)))}
sapply(0:10,FUN=calculateAUC,vector=results[,4])

sapply(0:20,FUN=calculateAUC,vector=results[,4])

for (i in 0:10) print(summary(results[(1+i*1000):(1000+i*1000),1]))
for (i in 0:10) print(summary(results[(1+i*1000):(1000+i*1000),2]))
for (i in 0:10) print(summary(results[(1+i*1000):(1000+i*1000),3]))
for (i in 0:10) print(summary(results[(1+i*1000):(1000+i*1000),4]))

for (i in 0:10) print(summary(results[(1+i*1000):(1000+i*1000),1]<=0.05))
for (i in 0:10) print(summary(results[(1+i*1000):(1000+i*1000),2]<=0.05))
for (i in 0:10) print(summary(results[(1+i*1000):(1000+i*1000),3]<=0.05))
for (i in 0:10) print(summary(results[(1+i*1000):(1000+i*1000),4]<=0.05))


for (i in 0:20) print(summary(results[(1+i*1000):(1000+i*1000),1]))
for (i in 0:20) print(summary(results[(1+i*1000):(1000+i*1000),2]))
for (i in 0:20) print(summary(results[(1+i*1000):(1000+i*1000),3]))
for (i in 0:20) print(summary(results[(1+i*1000):(1000+i*1000),4]))

for (i in 0:20) print(summary(results[(1+i*1000):(1000+i*1000),1]<=0.05))
for (i in 0:20) print(summary(results[(1+i*1000):(1000+i*1000),2]<=0.05))
for (i in 0:20) print(summary(results[(1+i*1000):(1000+i*1000),3]<=0.05))
for (i in 0:20) print(summary(results[(1+i*1000):(1000+i*1000),4]<=0.05))


sim_results_22k <- matrix(ncol=13,nrow=nrow(expr_data))
colnames(sim_results_22k) <- c("edgeR_pval","IMA_pr_pval","IMA_gb_pval","SoTA_fishers_pval","Ds_exprOnly","exprOnly_pval","Ds_prOnly","prOnly_pval","Ds_gbOnly","gbOnly_pval","Ds_2way_ind","pval_2way","indPGM_fishers_pval")
sim_results_22k <- as.data.frame(sim_results_22k)
sim_results_22k[,1] <- et$table[,3]


sitetest2 <- function (beta,group, gcase = "g2", gcontrol = "g1", testmethod = c("wilcox", "limma", "pooled", "satterthwaite"), Padj = "BH", concov = "OFF",rawpcut = NULL, adjustpcut = NULL, betadiffcut = NULL, paired = FALSE)
{
    
    grouplev = group[, 2]
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
                wilcox.test(x[1:length(lev1)], x[(length(lev1) + 
                                                      1):(length(lev1) + length(lev2))], paired = paired)$p.value
            })
        }
        if (testmethod == "limma") {
            require(limma)
            cat("Performing limma...\n")
            TS = as.factor(c(rep("T", length(lev1)), rep("C", 
                                                         length(lev2))))
            SS = rep(1:length(lev1), 2)
            if (paired == FALSE) {
                design = model.matrix(~0 + TS)
                rownames(design) = colnames(eset)
                colnames(design) = c("C", "T")
                fit = lmFit(eset, design)
                cont.matrix = makeContrasts(comp = T - C, levels = design)
                fit2 = contrasts.fit(fit, cont.matrix)
                fit2 = eBayes(fit2)
                result1 = topTable(fit2, coef = 1, adjust.method = Padj, 
                                   number = nrow(fit2))
            }
            else {
                design = model.matrix(~SS + TS)
                colnames(design) = c("Intercept", "Pairorder", 
                                     "TST")
                rownames(design) = colnames(eset)
                cat("Here is your design matrix\n")
                print(design)
                fit = lmFit(eset, design)
                fit2 = eBayes(fit)
                result1 = topTable(fit2, coef = "TST", adjust.method = Padj, 
                                   number = nrow(fit2))
            }
            testout = result1[match(rownames(eset), result1[, 
                                                            1]), "P.Value"]
        }
        if (testmethod == "pooled") {
            cat("Performing pooled t.test...\n")
            testout = apply(eset, 1, function(x) {
                t.test(x[1:length(lev1)], x[(length(lev1) + 
                                                 1):(length(lev1) + length(lev2))], var.equal = TRUE, 
                       paired = paired)$p.value
            })
        }
        if (testmethod == "satterthwaite") {
            cat("Performing satterthwaite t.test...\n")
            testout = apply(eset, 1, function(x) {
                t.test(x[1:length(lev1)], x[(length(lev1) + 
                                                 1):(length(lev1) + length(lev2))], paired = paired)$p.value
            })
        }
        adjustP = p.adjust(testout, method = Padj)
        difb = apply(eset, 1, function(x) {
            mean(x[1:length(lev1)]) - mean(x[(length(lev1) + 
                                                  1):ncol(eset)])
        })
        out = cbind(testout, adjustP, difb, rowMeans(eset[, 
                                                          1:length(lev1)]), rowMeans(eset[, (length(lev1) + 
                                                                                                 1):ncol(eset)]))
        rownames(out) = rownames(eset)
        colnames(out) = c("P-Value", "Adjust Pval", "Beta-Difference", 
                          paste("Mean", paste(gcase, collapse = "_"), sep = "_"), 
                          paste("Mean", paste(gcontrol, collapse = "_"), sep = "_"))
    }
    out = outputDMfunc(out = out, rawpcut = rawpcut, adjustpcut = adjustpcut, 
                       betadiffcut = betadiffcut)
    return(out)
}

