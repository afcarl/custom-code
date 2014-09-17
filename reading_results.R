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

pval_null_zscore <- NULL
means <- NULL
sds <- NULL
Ds <- NULL
Zs <- NULL
for (i in 11001:12000) {
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
summary(pval_null_zscore <= 0.05)

pval_null_zscore <- NULL
means <- NULL
sds <- NULL
Ds <- NULL
Zs <- NULL
for (i in 3341:4000) {
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
summary(pval_null_zscore <= 0.05)

pval_null_zscore <- NULL
means <- NULL
sds <- NULL
Ds <- NULL
Zs <- NULL
for (i in 2001:3000) {
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
summary(pval_null_zscore <= 0.05)

pval_null_zscore <- NULL
means <- NULL
sds <- NULL
Ds <- NULL
Zs <- NULL
for (i in 1:1000) {
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
summary(pval_null_zscore <= 0.05)