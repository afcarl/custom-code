function (indexlist, beta, indexmethod = c("mean", "median", 
    "tbrm")) 
{
    nr = length(indexlist$PID)
    temp2 = matrix(NA, nrow = nr, ncol = ncol(beta))
    rownames(temp2) = names(indexlist$SID)
    colnames(temp2) = colnames(beta)
    for (i in 1:nr) {
        temp = beta[indexlist$PID[[i]], ]
        if (length(indexlist$PID[[i]]) == 1) {
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
    return(temp2)
}
