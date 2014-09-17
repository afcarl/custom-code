findDiff <- function(x,y) {
	if (x < y) out <- y - x
	if (x > y) out <- -(x - y)
	return(out)
}