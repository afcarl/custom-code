eval(parse(text = paste(list11[i], "eset=dataf_PC@", list11[i], sep = "")))


for (i in 1:11) {
        cat("calculating", list11[i], "\n")
        eval(parse(text = paste("indexlist=dataf_PC@", list11[i], 
                                sep = "")))
	eval(parse(text = paste(list11[i], "eset=indexregionfunc(indexlist,dataf_PC@bmatrix,indexmethod)", sep = "")))
}