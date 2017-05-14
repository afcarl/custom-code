load("data/RData/data_LUSC_progressing.RData")

expr <- matrix(nrow = length(data_LUSC_progressing), ncol = nrow(data_LUSC_progressing[[1]]))
colnames(expr) <- rownames(data_LUSC_progressing[[1]])
rownames(expr) <- names(data_LUSC_progressing)
counts <- expr_plusOne <- expr
for (i in 1:length(data_LUSC_progressing)) { # extract expression data
	counts[i,] <- data_LUSC_progressing[[i]][,2]
	expr[i,] <- data_LUSC_progressing[[i]][,2]/data_LUSC_progressing[[i]][,1]
	expr_plusOne[i,] <- (data_LUSC_progressing[[i]][,2] + 1)/data_LUSC_progressing[[i]][,1]
}

for (i in 1:length(data_LUSC_progressing)) { # extract methylation data
	data_LUSC_progressing[[i]] <- cbind(data_LUSC_progressing[[i]], data_LUSC_progressing[[i]][,2]/data_LUSC_progressing[[i]][,1], (data_LUSC_progressing[[i]][,2] + 1)/data_LUSC_progressing[[i]][,1])
	colnames(data_LUSC_progressing[[i]])[ncol(data_LUSC_progressing[[i]])-1] <- "expr"
	colnames(data_LUSC_progressing[[i]])[ncol(data_LUSC_progressing[[i]])] <- "expr_plusOne"
	assign(names(data_LUSC_progressing)[i], data_LUSC_progressing[[i]])
}
lib_sizes <- data_LUSC_progressing[[1]][,1]
rm(i)
rm(data_LUSC_progressing)
save.image("data/RData/data_LUSC_progressing_lazyReady.RData")
quit()
n

# May 2017 update with new progression sets definitions:
load("data_KIRC_all.RData")
load("progressionIDs052017_KIRC.RData")
progressed <- intersect(progressed_KIRC, rownames(data_KIRC_all[[1]]))
nonProgressed <- intersect(nonProgressed_KIRC, rownames(data_KIRC_all[[1]]))

data_KIRC_progressing <- data_KIRC_all
expr <- matrix(nrow = length(data_KIRC_progressing), ncol = length(c(progressed, nonProgressed)))
colnames(expr) <- c(progressed, nonProgressed)
rownames(expr) <- names(data_KIRC_progressing)
counts <- expr_plusOne <- expr

for (i in 1:length(data_KIRC_progressing)) {
	data_KIRC_progressing[[i]] <- data_KIRC_progressing[[i]][c(progressed, nonProgressed),]
	counts[i,] <- data_KIRC_progressing[[i]][,2]
	expr[i,] <- data_KIRC_progressing[[i]][,2]/data_KIRC_progressing[[i]][,1]
	expr_plusOne[i,] <- (data_KIRC_progressing[[i]][,2] + 1)/data_KIRC_progressing[[i]][,1]
}
for (i in 1:length(data_KIRC_progressing)) { # extract methylation data
	data_KIRC_progressing[[i]] <- cbind(data_KIRC_progressing[[i]], data_KIRC_progressing[[i]][,2]/data_KIRC_progressing[[i]][,1], (data_KIRC_progressing[[i]][,2] + 1)/data_KIRC_progressing[[i]][,1])
	colnames(data_KIRC_progressing[[i]])[ncol(data_KIRC_progressing[[i]])-1] <- "expr"
	colnames(data_KIRC_progressing[[i]])[ncol(data_KIRC_progressing[[i]])] <- "expr_plusOne"
	assign(names(data_KIRC_progressing)[i], data_KIRC_progressing[[i]])
}
lib_sizes <- data_KIRC_progressing[[1]][,1]
rm(i, ANs, Ts, data_BRCA_progressing, data_BRCA_all, progressed_BRCA, nonProgressed_BRCA)

save.image("data_BRCA_progressing_lazyReady.RData")
quit()
n
# make lazyLoad-ready database - note - it doesn't work with too large objects (2.5G+)
# e = local({load("data/RData/data_KIRC_all_lazyReady.RData"); environment()})
# tools:::makeLazyLoadDB(e, "data/RData/data_KIRC_all_lazyReady")

# using SOAR package

library(SOAR)
load("data/RData/data_BRCA_progressing_lazyReady.RData")
Sys.setenv(R_LOCAL_CACHE="data/RData/data_BRCA_progressing_lazyReady.R_CACHE")
Store(ls())
quit()
n
# check that it works
library(SOAR)
Sys.setenv(R_LOCAL_CACHE="data/RData/data_BRCA_all_lazyReady.R_CACHE")
Objects()