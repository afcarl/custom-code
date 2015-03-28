require(ggplot2)
as.numeric.factor <- function(x) {as.numeric(levels(x))[x]}
cbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

min <- floor(min(c(dane1,dane2)))
max <- ceiling(max(c(dane1,dane2)))
hist1 <- hist(dane1,breaks=seq(min,max,1),plot=FALSE)
hist2 <- hist(dane2,breaks=seq(min,max,1),plot=FALSE)

df <- as.data.frame(rbind(cbind(hist1$mids,hist1$density,rep("MAS64",length(hist1$mids))),cbind(hist2$mids,hist2$density,rep("MAS75",length(hist1$mids)))))
colnames(df) <- c("size","height","polymer")
df$height <- as.numeric.factor(df$height)*100
df$size <- factor(df$size,levels=c("0.5","1.5","2.5","3.5","4.5","5.5","6.5","7.5","8.5","9.5","10.5","11.5","12.5"))
df$set <- as.character.factor(df$set)
ggplot(df,aes(x=size,y=height,fill=polymer)) +ggtitle("4")  + geom_bar(position="dodge",stat="identity") +ylab("Percentage") +xlab("particle size (um)") +scale_fill_manual(values=cbPalette)



# 3 populacje
min <- floor(min(c(dane1,dane2,dane3)))
max <- ceiling(max(c(dane1,dane2,dane3)))
hist1 <- hist(dane1,breaks=seq(min,max,1),plot=FALSE)
hist2 <- hist(dane2,breaks=seq(min,max,1),plot=FALSE)
hist3 <- hist(dane3,breaks=seq(min,max,1),plot=FALSE)
df <- as.data.frame(rbind(cbind(hist1$mids,hist1$density,rep("24h",length(hist1$mids))),cbind(hist2$mids,hist2$density,rep("48h",length(hist1$mids))),cbind(hist3$mids,hist3$density,rep("72h",length(hist1$mids)))))
colnames(df) <- c("size","height","time")
df$height <- as.numeric.factor(df$height)*100
df$size <- factor(df$size,levels=c("0.5","1.5","2.5","3.5","4.5","5.5","6.5","7.5","8.5","9.5","10.5","11.5","12.5"))
df$set <- as.character.factor(df$set)
ggplot(df,aes(x=size,y=height,fill=time)) + geom_bar(position="dodge",stat="identity") +ylab("Percentage") +xlab("particle size (um)") +scale_fill_manual(values=cbPalette)
