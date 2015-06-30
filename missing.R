found <- read.table("findings.txt",header=FALSE)
found <- found[,1]
all <- 1:17728
all[-which(all %in% found)]
