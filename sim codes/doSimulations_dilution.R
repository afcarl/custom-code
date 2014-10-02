library(edgeR)
library(IMA)
library(pROC)
source("sitetest2.R")
fishersMethod <- function(x) pchisq(-2 * sum(log(x)),df=2*length(x),lower=FALSE)
calculateAUC <- function(x,vector) {auc(predictor=c(vector[(1+x*1000):(1000+x*1000)],vector[(11001+x*1000):(12000+x*1000)]),response=c(rep("Pos",1000),rep("Neg",1000)))}
#calculateAUC <- function(x,vector) {auc(predictor=c(vector[(1+x*100):(100+x*100)],vector[(1101+x*100):(1200+x*100)]),response=c(rep("Pos",100),rep("Neg",100)))}

cat("sim pos\n")
source("simulating_data_mc_mvnorm_dilution.R")
sim.data_pos <- sim.data
cat("sim neg\n")
source("simulating_data_mc_mvnorm_dilution_neg.R")
sim.data <- c(sim.data_pos,sim.data)

ids <- c("001","002","003","004","005","006","007","008","009","010","011","012","013","014","015","016","017","018","019","020","021","022","023","024","025","026","027","028","029","030","031","032","033","034","035","036","037","038","039","040","041","042","043","044","045","046","047","048","049","050","051","052","053","054","055","056","057","058","059","060","061","062","063","064","065","066","067","068","069","070","071","072","073","074","075","076","077","078","079","080","081","082","083","084","085","086","087","088","089","090","091","092","093","094","095","096","097","098","099","100")
ANs <- paste("AN_",ids,sep="")
Ts <- paste("TU_",ids,sep="")

# expression
cat("do expr\n")
expr_data <- matrix(ncol=200,nrow=length(sim.data))
colnames(expr_data) <- c(ANs,Ts)
rownames(expr_data) <- paste("id_",seq(1,nrow(expr_data)),sep="")
for (i in 1:nrow(expr_data)) expr_data[i,] <- c(sim.data[[i]][[1]][,1],sim.data[[i]][[2]][,1])

group <- factor(c(rep("AN",100),rep("T",100)))
y <- DGEList(counts=expr_data,group=group)
y$samples[,2] <- 10^7
y$samples[,3] <- 1
y <- estimateCommonDisp(y)
y <- estimateTagwiseDisp(y)
et <- exactTest(y,dispersion="tagwise")

Ts <- paste("TU_",ids,sep="")
ANs <- paste("AN_",ids,sep="")
groupinfo <- as.data.frame(rbind(cbind(ANs,"AN"),cbind(Ts,"T")))
rownames(groupinfo) <- groupinfo[,1]
colnames(groupinfo) <- c("Sample.name","Group")

cat("do gb meth\n")
gb_meth <- matrix(ncol=200,nrow=length(sim.data))
for (i in 1:nrow(gb_meth)) gb_meth[i,] <- c(apply(sim.data[[i]][[1]][,4:14],1,mean),apply(sim.data[[i]][[2]][,4:14],1,mean))
rownames(gb_meth) <- paste("id_",seq(1,nrow(gb_meth)),sep="")
colnames(gb_meth) <- c(ANs,Ts)
gb_res <- sitetest2(gb_meth,groupinfo,gcase="T",gcontrol = "AN",testmethod = "satterthwaite")

cat("do pr meth\n")
pr_meth <- matrix(ncol=200,nrow=length(sim.data))
for (i in 1:nrow(pr_meth)) pr_meth[i,] <- c(apply(sim.data[[i]][[1]][,15:25],1,mean),apply(sim.data[[i]][[2]][,15:25],1,mean))
rownames(pr_meth) <- paste("id_",seq(1,nrow(pr_meth)),sep="")
colnames(pr_meth) <- c(ANs,Ts)
pr_res <- sitetest2(pr_meth,groupinfo,gcase="T",gcontrol = "AN",testmethod = "satterthwaite")

results <- cbind(gb_res[,1],pr_res[,1],et$table[,3])
results <- as.data.frame(results)
results$fishers <- apply(results,1,fishersMethod)

cat("gb\n")
print(sapply(0:10,FUN=calculateAUC,vector=results[,1]))
cat("pr\n")
print(sapply(0:10,FUN=calculateAUC,vector=results[,2]))
cat("expr\n")
print(sapply(0:10,FUN=calculateAUC,vector=results[,3]))
cat("joint\n")
print(sapply(0:10,FUN=calculateAUC,vector=results[,4]))
