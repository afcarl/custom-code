library(pROC)

calculateAUC <- function(x) {auc(predictor=c(pval_null_zscore[(1+x*1000):(1000+x*1000)],pval_null_zscore[(11001+x*1000):(12000+x*1000)]),response=c(rep("Pos",1000),rep("Neg",1000)))}

calculateSEN <- function(x) {length(which((pval_null_zscore[(1+x*1000):(1000+x*1000)]<=0.05)==TRUE))/1000}

calculateSPE <- function(x) {length(which((pval_null_zscore[(11001+x*1000):(11000+x*1000)]>0.05)==TRUE))/1000}

pval_null_zscore <- NULL
means <- NULL
sds <- NULL
Ds <- NULL
Zs <- NULL
for (i in 1:11000) {
	temp <- t(read.table(paste("./",i,".result",sep=""),nrow=1,skip=0))
	pval_null_zscore[[i]] <- temp[1,]
	Ds[[i]] <- temp[2,]
	means[[i]] <- temp[3,]
	sds[[i]] <- temp[4,]
	Zs[[i]] <- temp[5,]
	
	temp <- t(read.table(paste("./",i+11000,".result",sep=""),nrow=1,skip=0))
	pval_null_zscore[[i+11000]] <- temp[1,]
	Ds[[i+11000]] <- temp[2,]
	means[[i+11000]] <- temp[3,]
	sds[[i+11000]] <- temp[4,]
	Zs[[i+11000]] <- temp[5,]
}
pval_null_zscore <- unlist(pval_null_zscore)
means <- unlist(means)
sds <- unlist(sds)
Ds <- unlist(Ds)
Zs <- unlist(Zs)
rm(i)
rm(temp)

cat("sen: ")
sapply(0:10,FUN=calculateSEN)
cat("spe: ")
sapply(0:10,FUN=calculateSPE)
cat("AUC: ")
sapply(0:10,FUN=calculateAUC)