plotter1 <- data.frame(rank(-results_all$Zs_2way),rank(results_all$Fishers_combinedP),"all")
colnames(plotter1) <- c("PGM","SotA_fish","set")

plotter2 <- data.frame(rank(-results_all$Zs_2way)[which(workingList_BRCA %in% genes_fromPaper)],rank(results_all$Fishers_combinedP)[which(workingList_BRCA %in% genes_fromPaper)],"TCGA_BRCA_originalPaper")
colnames(plotter2) <- c("PGM","SotA_fish","set")

plotter3 <- data.frame(rank(-results_all$Zs_2way)[which(workingList_BRCA %in% genes_fromNature)],rank(results_all$Fishers_combinedP)[which(workingList_BRCA %in% genes_fromNature)],"Vogelstein_genes")
colnames(plotter3) <- c("PGM","SotA_fish","set")

plotter4 <- data.frame(rank(-results_all$Zs_2way)[which(workingList_BRCA %in% genes_fromCosmic)],rank(results_all$Fishers_combinedP)[which(workingList_BRCA %in% genes_fromCosmic)],"Cosmic_genes")
colnames(plotter4) <- c("PGM","SotA_fish","set")

ggplot(plotter1,aes(x=PGM,y=SotA_fish)) +theme_bw() + theme_bw() + geom_point(alpha=0.03) + geom_point(data=plotter2,aes(x=PGM,y=SotA_fish),colour="blue",alpha=0.5) + geom_point(data=plotter3,aes(x=PGM,y=SotA_fish),colour="green",alpha=0.5) + geom_point(data=plotter4,aes(x=PGM,y=SotA_fish),colour="red",alpha=0.5) + ggtitle("cor=0.725") + xlab("our full PGM") + ylab("edgeR + IMA results combined with Fisher's method")

ggplot(plotter1,aes(x=PGM,y=SotA_fish)) +theme_bw() + theme_bw() + geom_point(alpha=0.03) + ggtitle("cor=0.725") + xlab("rank according to our full PGM") + ylab("rank according to edgeR + IMA\n results combined with Fisher's method")
ggplot(plotter1,aes(x=PGM,y=SotA_fish)) +theme_bw() + theme_bw() + geom_point(alpha=0.03) + geom_point(data=plotter2,aes(x=PGM,y=SotA_fish),colour="blue",alpha=0.5) + ggtitle("cor=0.725") + xlab("rank according to our full PGM") + ylab("rank according to edgeR + IMA\n results combined with Fisher's method")
ggplot(plotter1,aes(x=PGM,y=SotA_fish)) +theme_bw() + theme_bw() + geom_point(alpha=0.03) + geom_point(data=plotter3,aes(x=PGM,y=SotA_fish),colour="green",alpha=0.5) + ggtitle("cor=0.725") + xlab("rank according to our full PGM") + ylab("rank according to edgeR + IMA\n results combined with Fisher's method")
ggplot(plotter1,aes(x=PGM,y=SotA_fish)) +theme_bw() + theme_bw() + geom_point(alpha=0.03) + geom_point(data=plotter4,aes(x=PGM,y=SotA_fish),colour="red",alpha=0.5) + ggtitle("cor=0.725") + xlab("rank according to our full PGM") + ylab("rank according to edgeR + IMA\n results combined with Fisher's method")