fac_2d <- function (matrix2,matrix1) {
	result <- matrix(ncol=ncol(matrix1),nrow=ncol(matrix1),data=rep(0,ncol(matrix1)*ncol(matrix2)))
	for (i in 1:nrow(matrix1)){ # iterate through each observation
		for (k in 1:ncol(matrix1)) { # iterate through each bin
			result[k,] <- result[k,] + matrix1[i,] * matrix2[i,k]
			result[,k] <- result[,k] + matrix1[i,k] * matrix2[i,]
		}
	}
	for (i in 1:20) result[,i] <- result[,i]/sum(result[,i])
	for (i in 1:20) result[i,] <- result[i,]/sum(result[i,])
	return(result)
}

tensor_product <- function(matrix1,matrix2) {
	if (is.matrix(matrix1) && is.matrix(matrix2)) result <- matrix(ncol=ncol(matrix1),nrow=ncol(matrix1),data=rep(0,ncol(matrix1)*ncol(matrix2))) else result <- matrix(ncol=length(matrix1),nrow=length(matrix1),data=rep(0,length(matrix1)*length(matrix2)))
	
	if (is.matrix(matrix1) && is.matrix(matrix2)) for (i in 1:nrow(matrix1)) {
		result <- result + matrix(nrow=20,ncol=20,byrow=TRUE,data=apply(expand.grid(matrix1[i,],matrix2[i,]), 1, prod))
	} else for (i in 1:length(matrix1)) {
		result <- result + matrix(nrow=20,ncol=20,byrow=TRUE,data=apply(expand.grid(matrix1,matrix2), 1, prod))
	}
	if (is.matrix(matrix1) && is.matrix(matrix2)) {
		for (i in 1:nrow(result)) result[i,] <- result[i,]/sum(result[i,])
	}
	return(result)
}