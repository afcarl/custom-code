norm_cpm_BRCA <- cpm_BRCA_plusOne
s2 <- 1
m2 <- 0

for (i in 1:nrow(norm_cpm_BRCA)) {
	xi <- norm_cpm_BRCA[i,ANs2]
	s1 <- sd(xi)
	m1 <- mean(xi)
	norm_cpm_BRCA[i,ANs2] <- m2 + (xi-m1) * (s2/s1)
	
	xi <- norm_cpm_BRCA[i,Ts2]
	s1 <- sd(xi)
	m1 <- mean(xi)
	norm_cpm_BRCA[i,Ts2] <- m2 + (xi-m1) * (s2/s1)
}


norm_bmatrix_BRCA_BODY <- bmatrix_BRCA_BODY2
for (i in 1:nrow(norm_bmatrix_BRCA_BODY)) {
	xi <- norm_bmatrix_BRCA_BODY[i,ANs2]
	s1 <- sd(xi)
	m1 <- mean(xi)
	norm_bmatrix_BRCA_BODY[i,ANs2] <- m2 + (xi-m1) * (s2/s1)
	
	xi <- norm_bmatrix_BRCA_BODY[i,Ts2]
	s1 <- sd(xi)
	m1 <- mean(xi)
	norm_bmatrix_BRCA_BODY[i,Ts2] <- m2 + (xi-m1) * (s2/s1)
}

norm_bmatrix_BRCA_PROMOTER <- bmatrix_BRCA_PROMOTER2
for (i in 1:nrow(norm_bmatrix_BRCA_PROMOTER)) {
	xi <- norm_bmatrix_BRCA_PROMOTER[i,ANs2]
	s1 <- sd(xi)
	m1 <- mean(xi)
	norm_bmatrix_BRCA_PROMOTER[i,ANs2] <- m2 + (xi-m1) * (s2/s1)
	
	xi <- norm_bmatrix_BRCA_PROMOTER[i,Ts2]
	s1 <- sd(xi)
	m1 <- mean(xi)
	norm_bmatrix_BRCA_PROMOTER[i,Ts2] <- m2 + (xi-m1) * (s2/s1)
}

temp <- cbind(norm_bmatrix_BRCA_PROMOTER[which,ANs2[1]],norm_bmatrix_BRCA_BODY[which,ANs2[1]],norm_cpm_BRCA[which,ANs2[1]])
for (i in 2:75){
	temp <- rbind(temp,cbind(norm_bmatrix_BRCA_PROMOTER[which,ANs2[i]],norm_bmatrix_BRCA_BODY[which,ANs2[i]],norm_cpm_BRCA[which,ANs2[i]]))
}
colnames(temp) <- c("PROMOTER","BODY","EXPRESSION")
temp <- as.data.frame(temp)

temp <- cbind(norm_bmatrix_BRCA_PROMOTER[which,Ts2[1]],norm_bmatrix_BRCA_BODY[which,Ts2[1]],norm_cpm_BRCA[which,Ts2[1]])
for (i in 2:75){
	temp <- rbind(temp,cbind(norm_bmatrix_BRCA_PROMOTER[which,Ts2[i]],norm_bmatrix_BRCA_BODY[which,Ts2[i]],norm_cpm_BRCA[which,Ts2[i]]))
}
colnames(temp) <- c("PROMOTER","BODY","EXPRESSION")
temp <- as.data.frame(temp)

tempAN <- cbind(norm_bmatrix_BRCA_PROMOTER[,ANs2[i]],norm_bmatrix_BRCA_BODY[,ANs2[i]],norm_cpm_BRCA[,ANs2[i]])
colnames(tempAN) <- c("PROMOTER","BODY","EXPRESSION")
tempAN <- as.data.frame(tempAN)
tempT <- cbind(norm_bmatrix_BRCA_PROMOTER[,Ts2[i]],norm_bmatrix_BRCA_BODY[,Ts2[i]],norm_cpm_BRCA[,Ts2[i]])
colnames(tempT) <- c("PROMOTER","BODY","EXPRESSION")
tempT <- as.data.frame(tempT)