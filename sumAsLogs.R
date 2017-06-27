sum_asLogs <- function (x) {
	if (x[1] >= x[2]) x[1] + log(1+ exp(x[2]-x[1])) else x[2] + log(1+ exp(x[1]-x[2]))
}