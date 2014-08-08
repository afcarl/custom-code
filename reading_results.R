Ds <- NULL
for (i in 1:42000) Ds[[i]] <- t(read.table(paste("./",i,".result",sep=""),nrow=1,skip=2))[,1]
rm(i)
Ds <- unlist(Ds)
summary(Ds[1:21000])
summary(Ds[21001:42000])

Ds <- NULL
for (i in 1:42000) Ds[[i]] <- t(read.table(paste("./",i,".result",sep=""),nrow=1,skip=3))[,1]
rm(i)
Ds <- unlist(Ds)
summary(Ds[1:21000])
summary(Ds[21001:42000])


pval_null_zscore <- NULL
means <- NULL
sds <- NULL
Ds <- NULL
for (i in 1:42000) {
	#Ds[[i]] <- t(read.table(paste("./",i,".result",sep=""),nrow=1,skip=2))[,1]
	temp <- t(read.table(paste("./",i,".result",sep=""),nrow=1,skip=0))
	pval_null_zscore[[i]] <- temp[1,]
	Ds[[i]] <- temp[2,]
	means[[i]] <- temp[3,]
	sds[[i]] <- temp[4,]
}
pval_null_zscore <- unlist(pval_null_zscore)
means <- unlist(means)
sds <- unlist(sds)
Ds <- unlist(Ds)

length(which(pval_null_zscore[1:21000]<=0.05))
length(which(pval_null_zscore[21001:42000]>0.05))

length(which(pval_null_chisq[1:21000]<=0.05))
length(which(pval_null_chisq[21001:42000]>0.05))


pval_null_chisq_neg <- NULL
pval_null_zscore_neg <- NULL
means_neg <- NULL
sds_neg <- NULL
Ds_neg <- NULL
for (i in 1:21000) {
	Ds_neg[[i]] <- t(read.table(paste("./",i,".result",sep=""),nrow=1,skip=2))[,1]
	temp <- t(read.table(paste("./",i,".result",sep=""),nrow=1,skip=3))
	pval_null_chisq_neg[[i]] <- temp[2,]
	pval_null_zscore_neg[[i]] <- temp[1,]
	means_neg[[i]] <- temp[3,]
	sds_neg[[i]] <- temp[4,]
}
pval_null_chisq_neg <- unlist(pval_null_chisq_neg)
pval_null_zscore_neg <- unlist(pval_null_zscore_neg)
means_neg <- unlist(means_neg)
sds_neg <- unlist(sds_neg)
Ds_neg <- unlist(Ds_neg)
pval_null_chisq[1:21000] <- pval_null_chisq_neg
pval_null_zscore[1:21000] <- pval_null_zscore_neg
means[1:21000] <- means_neg
sds[1:21000] <- sds_neg
Ds[1:21000] <- Ds_neg


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
for (i in 1:22000) {
	#Ds[[i]] <- t(read.table(paste("./",i,".result",sep=""),nrow=1,skip=2))[,1]
	temp <- t(read.table(paste("./",i,".result",sep=""),nrow=1,skip=0))
	pval_null_zscore[[i]] <- temp[1,]
	Ds[[i]] <- temp[2,]
	means[[i]] <- temp[3,]
	sds[[i]] <- temp[4,]
}
pval_null_zscore <- unlist(pval_null_zscore)
means <- unlist(means)
sds <- unlist(sds)
Ds <- unlist(Ds)

