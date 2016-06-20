sampNames <- t(Phenotype)[1,]
sampGroups <- as.factor(Phenotype[,2])

require(RColorBrewer)
pal = brewer.pal(8, "Dark2")
col <- pal[sampGroups]
pchs <- c(1,4,2)[sampGroups]

o_10000 <- order(-apply(bmatrix_pc,1,sd))[1:10000]
d <- dist(t(bmatrix_pc[o_10000,]))
fit <- cmdscale(d)
xlim <- range(fit[, 1]) * 1.2
ylim <- range(fit[, 2]) * 1.2

main <- sprintf("Beta MDS\n10,000 most variable probes")
plot(fit, type = "p", xlim = xlim, ylim = ylim, xlab = "",ylab = "", main = main, pch=pchs)
legend("topleft", legend = levels(sampGroups), ncol = 3, pch=c(1,4,2))

