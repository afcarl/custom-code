SotA_Fish <- c(0.994639,0.989985,0.974896,0.951959,0.927059,0.860719,0.772325,0.6761,0.578703,0.506866,0.494861)
PGM_Fish <- c(0.999957,0.999406,0.997602,0.990097,0.979989,0.950004,0.877052,0.785053,0.626624,0.540399,0.500966)
PGM_full <- c(0.99995,0.999624,0.998532,0.995332,0.984446,0.94928,0.879564,0.771318,0.628272,0.539664,0.508325)

plotter_auc <- as.data.frame(rbind(cbind(seq(0,1,0.1),PGM_full,"full_PGM"),cbind(seq(0,1,0.1),PGM_Fish,"fishersCombined_PGM"),cbind(seq(0,1,0.1),SotA_Fish,"fishersCombined_SotA")))
colnames(plotter_auc) <- c("fraction_tumours","AUC","method")
plotter_auc[,1] <- as.numeric.factor(plotter_auc[,1])
plotter_auc[,2] <- as.numeric.factor(plotter_auc[,2])
ggplot(plotter_auc,aes(x=fraction_tumours,y=AUC)) +theme_bw() + geom_line(aes(colour=method),lwd=0.5) + geom_point(aes(colour=method)) + xlab("fraction of true tumours") + scale_x_continuous(labels=c("1","0.75","0.5","0.25","0")) + theme(legend.position="bottom")

SotA_Fish <- c(0.638291,0.648743,0.634719,0.661737,0.644815,0.627265,0.633608,0.624091,0.649164,0.632478,0.649038)
PGM_Fish <- c(0.690236,0.682235,0.657925,0.673176,0.664372,0.695690,0.694812,0.673271,0.666555,0.650082,0.684331)
PGM_full <-c(0.787947,0.795443,0.823789,0.886148,0.944583,0.979411,0.992761,0.998194,0.999578,0.999770,1.000000)

plotter_auc <- as.data.frame(rbind(cbind(seq(0,1,0.1),PGM_full,"full_PGM"),cbind(seq(0,1,0.1),PGM_Fish,"fishersCombined_PGM"),cbind(seq(0,1,0.1),SotA_Fish,"fishersCombined_SotA")))
colnames(plotter_auc) <- c("delta_cor","AUC","method")
plotter_auc[,1] <- as.numeric.factor(plotter_auc[,1])
plotter_auc[,2] <- as.numeric.factor(plotter_auc[,2])
ggplot(plotter_auc,aes(x=delta_cor,y=AUC)) +theme_bw() + geom_line(aes(colour=method),lwd=0.5) + geom_point(aes(colour=method)) + xlab("delta(correlation)") + ylim(c(0.5,1)) + theme(legend.position="bottom")


# using facets

SotA_Fish <- c(0.994639,0.989985,0.974896,0.951959,0.927059,0.860719,0.772325,0.6761,0.578703,0.506866,0.494861)
PGM_Fish <- c(0.999957,0.999406,0.997602,0.990097,0.979989,0.950004,0.877052,0.785053,0.626624,0.540399,0.500966)
PGM_full <- c(0.99995,0.999624,0.998532,0.995332,0.984446,0.94928,0.879564,0.771318,0.628272,0.539664,0.508325)
plotter_auc <- as.data.frame(cbind(rbind(cbind(seq(0,1,0.1),PGM_full,"full_PGM"),cbind(seq(0,1,0.1),PGM_Fish,"fishersCombined_PGM"),cbind(seq(0,1,0.1),SotA_Fish,"fishersCombined_SotA")),rep("Dilution analysis",33)))

SotA_Fish <- c(0.638291,0.648743,0.634719,0.661737,0.644815,0.627265,0.633608,0.624091,0.649164,0.632478,0.649038)
PGM_Fish <- c(0.690236,0.682235,0.657925,0.673176,0.664372,0.695690,0.694812,0.673271,0.666555,0.650082,0.684331)
PGM_full <-c(0.787947,0.795443,0.823789,0.886148,0.944583,0.979411,0.992761,0.998194,0.999578,0.999770,1.000000)

plotter_auc <- rbind(plotter_auc,as.data.frame(cbind(rbind(cbind(seq(0,1,0.1),PGM_full,"full_PGM"),cbind(seq(0,1,0.1),PGM_Fish,"fishersCombined_PGM"),cbind(seq(0,1,0.1),SotA_Fish,"fishersCombined_SotA")),rep("Delta correlation analysis",33))))

colnames(plotter_auc) <- c("variable","AUC","method","analysis")
plotter_auc[,1] <- as.numeric.factor(plotter_auc[,1])
plotter_auc[,2] <- as.numeric.factor(plotter_auc[,2])

ggplot(plotter_auc,aes(x=variable,y=AUC)) + facet_wrap(facets=~analysis) +theme_bw() + geom_line(aes(colour=method),lwd=0.5) + geom_point(aes(colour=method)) + theme(legend.position="right", panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.background = element_rect(fill="white"))