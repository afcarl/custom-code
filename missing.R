found <- read.table("findings.txt",header=FALSE)
found <- found[,1]
all <- c(3001:4000,11001:12000
all[-which(all %in% found)]
