runmean <- function (x, k,endrule = c("mean", "NA", "trim", "keep", "constant", "func"), align = c("center","left", "right"))
{
	EndRule <- function (x, y, k, dimx, endrule = c("NA", "trim", "keep", "constant", "func"), align = c("center", "left", "right"), Func, ...) {
	align = match.arg(align)
	k = as.integer(k)
	k2 = k%/%2
	if (k2 < 1) 
		k2 = 1
	yIsVec = is.null(dimx)
	if (yIsVec) 
		dimx = c(length(y), 1)
	dim(x) <- dimx
	dim(y) <- dimx
	n = nrow(x)
	m = ncol(x)
	if (k > n) 
		k2 = (n - 1)%/%2
	k1 = k - k2 - 1
	if (align == "center" && k == 2) 
	align = "right"
	if (endrule == "trim") {
		y = y[(k1 + 1):(n - k2), ]
	}
	else if (align == "center") {
		idx1 = 1:k1
		idx2 = (n - k2 + 1):n
		if (endrule == "NA") {
			y[idx1, ] = NA
			y[idx2, ] = NA
		}
		else if (endrule == "keep") {
			y[idx1, ] = x[idx1, ]
			y[idx2, ] = x[idx2, ]
		}
		else if (endrule == "constant") {
			y[idx1, ] = y[k1 + 1 + integer(m), ]
			y[idx2, ] = y[n - k2 + integer(m), ]
		}
		else if (endrule == "func" || !yIsVec) {
			for (j in 1:m) {
				for (i in idx1) y[i, j] = Func(x[1:(i + k2),j], ...)
				for (i in idx2) y[i, j] = Func(x[(i - k1):n,j], ...)
			}
		}
	}
	else if (align == "left") {
		y[1:(n - k1), ] = y[(k1 + 1):n, ]
		idx = (n - k + 2):n
		if (endrule == "NA") {
			y[idx, ] = NA
		}
		else if (endrule == "keep") {
			y[idx, ] = x[idx, ]
		}
		else if (endrule == "constant") {
			y[idx, ] = y[n - k + integer(m) + 1, ]
		}
		else {
			for (j in 1:m) for (i in idx) y[i, j] = Func(x[i:n,j], ...)
		}
	}
	else if (align == "right") {
		y[(k2 + 1):n, ] = y[1:(n - k2), ]
		idx = 1:(k - 1)
		if (endrule == "NA") {
			y[idx, ] = NA
		}
		else if (endrule == "keep") {
			y[idx, ] = x[idx, ]
		}
		else if (endrule == "constant") {
			y[idx, ] = y[k + integer(m), ]
		}
		else {
			for (j in 1:m) for (i in idx) y[i, j] = Func(x[1:i,j], ...)
		}
	}
	if (yIsVec) 
		y = as.vector(y)
	return(y)
}
	
	endrule = match.arg(endrule)
	align = match.arg(align)
	dimx = dim(x)
	x = as.vector(x)
	n = length(x)
	if (k <= 1) 
		return(x)
	if (k > n) 
		k = n
	k2 = k%/%2
	y = double(n)
	k1 = k - k2 - 1
	y = c(sum(x[1:k]), diff(x, k))
	y = cumsum(y)/k
	y = c(rep(0, k1), y, rep(0, k2))
	if (endrule == "mean") 
		endrule = "func"
	y = EndRule(x, y, k, dimx, endrule, align, mean, na.rm = TRUE)
	return(y)
}
