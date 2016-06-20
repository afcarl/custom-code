peak.correction <- function (data, anno) 
{
    anno <- anno[row.names(data), ]
    TI <- anno[, "INFINIUM_DESIGN_TYPE"] == "I"
    TII <- anno[, "INFINIUM_DESIGN_TYPE"] == "II"
    corrected.data <- apply(data, 2, function(B) {
        SI <- summits(B[TI])
        SII <- summits(B[TII])
        BI <- correctI(as.vector(B[TI]), SI, SII)
		BII <- correctII(as.vector(B[TII]), SI, SII)
        return(c(BI, BII))
    })
    row.names(corrected.data) <- c(row.names(data[TI, ]), row.names(data[TII, ]))
    return(corrected.data)
}
correctI <- function (BetaValues, SI, SII) 
{
    return(BetaValues)
}

correctII <- function (BetaValues, SI, SII) 
{
    M <- Beta2M(BetaValues)
    sigma_u <- SII[1]/SI[1]
    sigma_m <- SII[2]/SI[2]
    M <- sapply(M, function(x) {
        if (is.na(x)) return(NA)  ##LS##
        if (x < 0) 
            return(x/sigma_u)
        else return(x/sigma_m)
    })
    return(M2Beta(M))
}

summits <- function (BetaValues) 
{
    d <- density(BetaValues,na.rm=TRUE)
    yneg <- d$y[1:which(d$x > M2Beta(0))[1]]
    ypos <- d$y[which(d$x > M2Beta(0))[1]:length(d$y)]
    sa <- d$x[which(d$y == max(yneg))]
    sb <- d$x[which(d$y == max(ypos))]
    return(c(Beta2M(sa), Beta2M(sb)))
}

M2Beta <- function (M) 
{
	return((2^M)/(2^M + 1))
}

Beta2M <- function (B) 
{
	return(log2(B/(1 - B)))
}
