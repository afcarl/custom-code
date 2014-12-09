pval_null_zscore <- NULL
means <- NULL
sds <- NULL
Ds <- NULL
Zs <- NULL
for (i in 1:21000) {
	temp <- t(read.table(paste("./",i,".result",sep=""),nrow=1,skip=0))
	pval_null_zscore[[i]] <- temp[1,]
	Ds[[i]] <- temp[2,]
	means[[i]] <- temp[3,]
	sds[[i]] <- temp[4,]
	Zs[[i]] <- temp[5,]
	temp <- t(read.table(paste("./",i+21000,".result",sep=""),nrow=1,skip=0))
	pval_null_zscore[[i+21000]] <- temp[1,]
	Ds[[i+21000]] <- temp[2,]
	means[[i+21000]] <- temp[3,]
	sds[[i+21000]] <- temp[4,]
	Zs[[i+21000]] <- temp[5,]
}
pval_null_zscore <- unlist(pval_null_zscore)
means <- unlist(means)
sds <- unlist(sds)
Ds <- unlist(Ds)
Zs <- unlist(Zs)
rm(i)
rm(temp)

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

pval_null_zscore <- NULL
means <- NULL
sds <- NULL
Ds <- NULL
Zs <- NULL
for (i in 1:17728) {
	temp <- t(read.table(paste("./",i,".result",sep=""),nrow=1,skip=0))
	pval_null_zscore[[i]] <- temp[1,]
	Ds[[i]] <- temp[2,]
	means[[i]] <- temp[3,]
	sds[[i]] <- temp[4,]
	Zs[[i]] <- temp[5,]
}
pval_null_zscore <- unlist(pval_null_zscore)
means <- unlist(means)
sds <- unlist(sds)
Ds <- unlist(Ds)
Zs <- unlist(Zs)
rm(i)
rm(temp)


# read x-val data

Zs <- matrix(nrow=17728,ncol=14)
mlogliks <- NULL
for (i in 1:17728) {
	Zs[i,] <- t(read.table(paste("./",i,".result",sep=""),nrow=1,skip=0))
	mlogliks[[i]] <- read.table(paste("./",i,".result",sep=""),skip=1,header=TRUE)
}
rm(i)


targets <- read.delim("commandfile_predictions.txt",header=F)
mlogliks <- NULL
for (i in 1:100) {
	mlogliks[[i]] <- read.table(paste("./",targets[i,2],".predicted",sep=""),header=TRUE)
}
rm(i)
