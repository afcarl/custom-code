# ls *.result > ready.txt
ready <- read.table(file="ready.txt")
pvals <- vector()
for (i in 1:nrow(ready)) pvals[i] <- read.table(file=paste("./",ready[i,1], sep=""),nrow=1,skip=1)[1]
pvals <- unlist(pvals)
summary(pvals)