function (indexlist, subset, beta, indexmethod = c("mean", "median", 
    "tbrm")) 
{
    nr = length(indexlist$PID)
    temp2 = matrix(NA, nrow = nr, ncol = ncol(beta))
    rownames(temp2) = names(indexlist$SID)
    colnames(temp2) = colnames(beta)
    for (i in 1:nr) {
        var = indexlist$PID[[i]][indexlist$PID[[i]] %in% subset]
		temp = beta[var, ]
		if (length(var == 0)) {
            temp2[i, ] = NA
        }
        if (length(var == 1)) {
            temp2[i, ] = temp
        }
        else {
            if (indexmethod == "tbrm") {
                temp2[i, ] = apply(temp, 2, eval(indexmethod))
            }
            else {
                temp2[i, ] = apply(temp, 2, eval(indexmethod), 
                  na.rm = TRUE)
            }
        }
    }
	temp2=temp2[!is.na(temp2[,1])]
    return(temp2)
}
