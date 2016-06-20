meth_Optims <- function(data1,data2,noBreaks,maxit=5000,method="Nelder-Mead") {
	flanks <- c(-7.01,7.01)
	final_KLs <- vector(length=noBreaks,mode="numeric")
	return <- NULL
	
	KL <- function(newBreaks) {
		########## preps #######################
		epsilon <- 1/length(c(data1,data2))
		breaks <- sort(c(flanks,newBreaks))
		all_labels <- as.character(seq(1,length(breaks)-1,1))
		############## distr. 1 ###################
		dens1 <- density(data1,bw=0.14,from=-7,to=7,n=1401)
		dens1$y <- dens1$y/sum(dens1$y)
		dens1 <- as.data.frame(cbind(dens1$x,dens1$y))
		colnames(dens1) <- c("X","density")
		dens1$interval <- findInterval(dens1[,1],breaks)
		dens1 <- aggregate(cbind(X,density) ~ interval, data=dens1, sum)
		dens1$interval <- as.factor(dens1$interval)
		labels <- colnames(t(summary(dens1$interval)))
		frequencies1_t <- dens1$density
		frequencies1 <- rep(0,length(breaks)-1)
		for (label in 1:length(all_labels)) {
			if (!is.na(match(all_labels[label], labels)))
				frequencies1[label] <- frequencies1_t[match(all_labels[label], labels)]
			}
		frequencies1 <- frequencies1 + epsilon
		frequencies1 <- frequencies1/sum(frequencies1)
		############## distr. 2 ###################
		dens2 <- density(data2,bw=0.14,from=-7,to=7,n=1401)
		dens2$y <- dens2$y/sum(dens2$y)
		dens2 <- as.data.frame(cbind(dens2$x,dens2$y))
		colnames(dens2) <- c("X","density")
		dens2$interval <- findInterval(dens2[,1],breaks)
		dens2 <- aggregate(cbind(X,density) ~ interval, data=dens2, sum)
		dens2$interval <- as.factor(dens2$interval)
		labels <- colnames(t(summary(dens2$interval)))
		frequencies2_t <- dens2$density
		frequencies2 <- rep(0,length(breaks)-1)
		for (label in 1:length(all_labels)) {
			if (!is.na(match(all_labels[label], labels)))
				frequencies2[label] <- frequencies2_t[match(all_labels[label], labels)]
			}
		frequencies2 <- frequencies2 + epsilon
		frequencies2 <- frequencies2/sum(frequencies2)
		############## KL calculation ###################
		KL <- sum(log(frequencies1/frequencies2)*frequencies1)
		return(-KL)
	}
	
	control <- NULL
	control$maxit <- maxit
	breaks_sets <- NULL
	for (i in 1:noBreaks) {
		breaksProposal <- seq(from=-7.01,to=7.01,by=14.02/(i+1))[-c(1,i+2)]
		temp <- suppressWarnings(optim(breaksProposal,lower=-7,upper=7,KL,method=method,control=control))
		final_KLs[i] <- -temp$value
		breaks_sets[[i]] <- sort(c(flanks,temp$par))
	}
	return <- NULL
	return$breaks <- breaks_sets
	return$KLs <- final_KLs
	return(return)
}


meth_greedyOptims <- function(data1,data2,noBreaks,maxit=5000,method="Nelder-Mead") {
  breaks <- c(-7.01,7.01)
  final_KLs <- vector(length=noBreaks,mode="numeric")
  return <- NULL
  
	KL2 <- function(newBreak) {
		########## preps #######################
		epsilon <- 1/length(c(data1,data2))
		breaks <- unique(sort(c(breaks,newBreak)))
		all_labels <- as.character(seq(1,length(breaks)-1,1))
		############## distr. 1 ###################
		dens1 <- density(data1,bw=0.14,from=-7,to=7,n=1401)
		dens1$y <- dens1$y/sum(dens1$y)
		dens1 <- as.data.frame(cbind(dens1$x,dens1$y))
		colnames(dens1) <- c("X","density")
		dens1$interval <- findInterval(dens1[,1],breaks)
		dens1 <- aggregate(cbind(X,density) ~ interval, data=dens1, sum)
		dens1$interval <- as.factor(dens1$interval)
		labels <- colnames(t(summary(dens1$interval)))
		frequencies1_t <- dens1$density
		frequencies1 <- rep(0,length(breaks)-1)
		for (label in 1:length(all_labels)) {
			if (!is.na(match(all_labels[label], labels)))
				frequencies1[label] <- frequencies1_t[match(all_labels[label], labels)]
			}
		frequencies1 <- frequencies1 + epsilon
		frequencies1 <- frequencies1/sum(frequencies1)
		############## distr. 2 ###################
		dens2 <- density(data2,bw=0.14,from=-7,to=7,n=1401)
		dens2$y <- dens2$y/sum(dens2$y)
		dens2 <- as.data.frame(cbind(dens2$x,dens2$y))
		colnames(dens2) <- c("X","density")
		dens2$interval <- findInterval(dens2[,1],breaks)
		dens2 <- aggregate(cbind(X,density) ~ interval, data=dens2, sum)
		dens2$interval <- as.factor(dens2$interval)
		labels <- colnames(t(summary(dens2$interval)))
		frequencies2_t <- dens2$density
		frequencies2 <- rep(0,length(breaks)-1)
		for (label in 1:length(all_labels)) {
			if (!is.na(match(all_labels[label], labels)))
				frequencies2[label] <- frequencies2_t[match(all_labels[label], labels)]
			}
		frequencies2 <- frequencies2 + epsilon
		frequencies2 <- frequencies2/sum(frequencies2)
		############## KL calculation ###################
		KL <- sum(log(frequencies1/frequencies2)*frequencies1)
		return(-KL)
	}
  
  control <- NULL
  control$maxit <- maxit
  breaks_sets <- NULL
  for (i in 1:noBreaks) {
    breakProposals <- vector(length=100,mode="numeric")
    proposal_KLs <- vector(length=100,mode="numeric")
    for (j in 1:100) {
      breakProposals[j] <- runif(1,min=min(c(data1,data2)),max=max(c(data1,data2)))
      proposal_KLs[j] <- -KL2(breakProposals[j])
    }
    breakProposal <- breakProposals[which(proposal_KLs %in% max(proposal_KLs))[1]]
    temp <- suppressWarnings(optim(breakProposal,lower=min(c(data1,data2)),upper=max(c(data1,data2)),KL2,method=method,control=control))
    final_KLs[i] <- -temp$value
    breaks <- sort(c(breaks,temp$par))
    breaks_sets[[i]] <- breaks
  }
  return <- NULL
  return$breaks <- breaks_sets
  return$KLs <- final_KLs
  return(return)
}

integrand <- function(x,k) {dpois(k,x)}
expr_greedyOptims <- function(data1,data2,noBreaks,maxit=5000,method="Nelder-Mead") {
	breaks <- c(0,2*max(c(data1,data2))/10)
	final_KLs <- vector(length=noBreaks,mode="numeric")
	return <- NULL
  
	KL3 <- function(newBreak) {
		########## preps #######################
		epsilon <- 1/length(c(data1,data2))
		breaks <- unique(sort(c(breaks,newBreak)))
		all_labels <- as.character(seq(1,length(breaks)-1,1))
		############## distr. 1 ###################
		temp1_d <- temp1
		temp1_d$interval <- findInterval(temp1_d[,1],breaks)
		temp1_d <- aggregate(cbind(cpm,density) ~ interval, data=temp1_d, sum)
		temp1_d$interval <- as.factor(temp1_d$interval)
		labels <- colnames(t(summary(temp1_d$interval)))
		frequencies1_t <- temp1_d$density
		frequencies1 <- rep(0,length(breaks)-1)
		for (label in 1:length(all_labels)) {
			if (!is.na(match(all_labels[label], labels)))
				frequencies1[label] <- frequencies1_t[match(all_labels[label], labels)]
			}
		frequencies1 <- frequencies1 + epsilon
		frequencies1 <- frequencies1/sum(frequencies1)
		############## distr. 2 ###################
		temp2_d <- temp2
		temp2_d$interval <- findInterval(temp2_d[,1],breaks)
		temp2_d <- aggregate(cbind(cpm,density) ~ interval, data=temp2_d, sum)
		temp2_d$interval <- as.factor(temp2_d$interval)
		labels <- colnames(t(summary(temp2_d$interval)))
		frequencies2_t <- temp2_d$density
		frequencies2 <- rep(0,length(breaks)-1)
		for (label in 1:length(all_labels)) {
			if (!is.na(match(all_labels[label], labels)))
				frequencies2[label] <- frequencies2_t[match(all_labels[label], labels)]
			}
		frequencies2 <- frequencies2 + epsilon
		frequencies2 <- frequencies2/sum(frequencies2)
		############## KL calculation ###################
		KL <- sum(log(frequencies1/frequencies2)*frequencies1)
		return(-KL)
	}
	
	temp1 <- matrix(ncol=2)
	temp2 <- matrix(ncol=2)
	colnames(temp1) <- c("cpm","density")
	colnames(temp2) <- c("cpm","density")
	for (i in 1:length(data1)) {
		lambda <- data1[i]
		X <- seq(round(max(lambda-(4*lambda*lambda^(-1/2)),1)),round(lambda+(4*lambda*lambda^(-1/2))))
		current <- 10
		temp1 <- rbind(temp1,cbind(X/current,dpois(X,lambda=lambda)*current))
	}
	temp1 <- as.data.frame(temp1[-c(1,which(is.na(temp1[,1]))),])
	temp1 <- temp1[order(temp1$cpm),]
	temp1[,3] <- cumsum(temp1[,2])
	temp1[,3] <- temp1[,3]/max(temp1[,3])
	for (i in 1:length(data2)) {
		lambda <- data2[i]
		X <- seq(round(max(lambda-(4*lambda*lambda^(-1/2)),1)),round(lambda+(4*lambda*lambda^(-1/2))))
		current <- 10
		temp2 <- rbind(temp2,cbind(X/current,dpois(X,lambda=lambda)*current))
	}
	temp2 <- as.data.frame(temp2[-c(1,which(is.na(temp2[,1]))),])
	temp2 <- temp2[order(temp2$cpm),]
	temp2[,3] <- cumsum(temp2[,2])
	temp2[,3] <- temp2[,3]/max(temp2[,3])
  
  control <- NULL
  control$maxit <- maxit
  breaks_sets <- NULL
  for (i in 1:noBreaks) {
    breakProposals <- vector(length=100,mode="numeric")
    proposal_KLs <- vector(length=100,mode="numeric")
    for (j in 1:10) {
      breakProposals[j] <- runif(1,min=min(c(data1,data2))/10,max=max(c(data1,data2))/10)
      proposal_KLs[j] <- -KL3(breakProposals[j])
    }
    breakProposal <- breakProposals[which(proposal_KLs %in% max(proposal_KLs))[1]]
    temp <- suppressWarnings(optim(breakProposal,lower=min(c(data1,data2))/10,upper=max(c(data1,data2))/10,KL3,method=method,control=control))
    final_KLs[i] <- -temp$value
    breaks <- sort(c(breaks,temp$par))
    breaks_sets[[i]] <- breaks
  }
  return <- NULL
  return$breaks <- breaks_sets
  return$KLs <- final_KLs
  return(return)
}


########### Sudhakar's ################################

find_greedy_breaks <- function(data1,data2,noBreaks,resolution,maxBreak) {
	# data1 and data2 are vectors with the values to be discretized
	# noBreaks is the requested number of binning boundaries - i.e requesting 1 break will result in 2 bins.
	# resolution is the value at which the range of values will be probed at in the greedy evaluation
	# maxBreak is there to allow Sudhakar to specify the maximal value his data can take - in his case it was 1.
	breaks <- sort(c(min(c(data1,data2)),max(c(data1,data2)))) # initialize breaks to be the min and max of the data
	final_KLs <- vector(length=noBreaks,mode="numeric") # a variable holding the KL values of selected breaks at each iteration
	breaks_sets <- NULL # a list containing selected break sets at each iteration 
	proposed_breaks <- seq(min(c(data1,data2)),maxBreak,resolution)
	for (i in 1:noBreaks) {
		which <- which(proposed_breaks %in% breaks)
		if (length(which) > 0) proposed_breaks <- proposed_breaks[-which] # make sure the proposed breaks do not happen to be positioned at the data. also, remove the already selected break points.
		KLs <- vector(length=length(proposed_breaks),mode="numeric")
		for (j in 1:length(proposed_breaks)) {
			KLs[j] <- KL(data1,data2,sort(c(breaks,proposed_breaks[j])))
		}
		which <- which(KLs %in% max(KLs))[1]
		final_KLs[i] <- KLs[which]
		breaks <- sort(c(breaks,proposed_breaks[which]))
		breaks_sets[[i]] <- breaks
	}
	return_list <- NULL
	return_list$breaks <- breaks_sets
	return_list$KLs <- final_KLs
	return(return_list)
}

KL <- function(data1,data2,breaks) {
	epsilon <- 1/length(c(data1,data2)) # we need to add epsilon to prevent KL inflation due to 0s
	hist1 <- hist(data1,breaks=breaks,plot=FALSE)
	hist1$counts <- (hist1$counts+epsilon)/sum(hist1$counts+epsilon)
	hist2 <- hist(data2,breaks=breaks,plot=FALSE)
	hist2$counts <- (hist2$counts+epsilon)/sum(hist2$counts+epsilon)
	KL <- sum(log(hist1$counts/hist2$counts)*hist1$counts)
	return(KL)
}

greedy_breaks <- find_greedy_breaks(rnorm(100, mean=0),rnorm(100, mean=1),noBreaks=19,resolution=0.01,maxBreak=5)

KL2 <- function(breaks) {
  epsilon <- 1/length(c(data1,data2))
	breaks <- c(min(c(data1,data2))-epsilon,breaks,max(c(data1,data2))+epsilon)
	hist1 <- hist(data1,breaks=breaks,plot=FALSE)
	hist1$density <- (hist1$density+epsilon)/sum(hist1$density+epsilon)
	hist2 <- hist(data2,breaks=breaks,plot=FALSE)
	hist2$density <- (hist2$density+epsilon)/sum(hist2$density+epsilon)
	KL <- sum(log(hist1$density/hist2$density)*hist1$density)
	return(-KL)
}

greedy_optim <- function(data1,data2,noBreaks,maxit=500,method="Nelder-Mead") {
  breaks <- c(min(c(data1,data2)),max(c(data1,data2)))
  final_KLs <- vector(length=noBreaks,mode="numeric")
  return <- NULL
  
  KL2 <- function(newBreak) {
    epsilon <- 1/length(c(data1,data2))
    breaks <- unique(sort(c(breaks,newBreak)))
    hist1 <- hist(data1,breaks=breaks,plot=FALSE)
    hist1$counts <- (hist1$counts+epsilon)/sum(hist1$counts+epsilon)
    hist2 <- hist(data2,breaks=breaks,plot=FALSE)
    hist2$counts <- (hist2$counts+epsilon)/sum(hist2$counts+epsilon)
    KL <- sum(log(hist1$counts/hist2$counts)*hist1$counts)
    return(-KL)
  }
  
  control <- NULL
  control$maxit <- maxit
  breaks_sets <- NULL
  for (i in 1:noBreaks) {
    breakProposals <- vector(length=100,mode="numeric")
    proposal_KLs <- vector(length=100,mode="numeric")
    for (j in 1:100) {
      breakProposals[j] <- runif(1,min=min(c(data1,data2)),max=max(c(data1,data2)))
      proposal_KLs[j] <- -KL2(breakProposals[j])
    }
    breakProposal <- breakProposals[which(proposal_KLs %in% max(proposal_KLs))[1]]
    temp <- suppressWarnings(optim(breakProposal,lower=min(c(data1,data2)),upper=max(c(data1,data2)),KL2,method=method,control=control))
    final_KLs[i] <- -temp$value
    breaks <- sort(c(breaks,temp$par))
    breaks_sets[[i]] <- breaks
  }
  return <- NULL
  return$breaks <- breaks_sets
  return$KLs <- final_KLs
  return(return)
}


optim_breaks <- function(data1,data2,noBreaks,maxit=500) {
	epsilon <- 1e-04
	final_KLs <- vector(length=noBreaks,mode="numeric")
	breaks_sets <- NULL
	control <- NULL
	control$maxit <- maxit
	for (j in 1:noBreaks) {
		proposal <- seq(min(c(data1,data2))-epsilon,max(c(data1,data2))+epsilon,(max(c(data1,data2))-min(c(data1,data2))+2*epsilon)/(j+1))[-c(1,j+2)]
		temp <- optim(proposal,KL2,method="Nelder-Mead",control=control)
		final_KLs[j] <- -temp$value
		breaks_sets[[j]] <- c(min(c(data1,data2))-epsilon,temp$par,max(c(data1,data2))+epsilon)
	}
	return <- NULL
	return$breaks <- breaks_sets
	return$KLs <- final_KLs
	return(return)
}

greedy_breaks <- greedy_breaks(data1,data2,noBreaks=20,resolution=0.01,maxBreak=2.5)
breaks <- greedy_breaks$breaks[[which(greedy_breaks$KLs %in% max(greedy_breaks$KLs))]]
for(i in 1:length(breaks)) cat(breaks[i],"\n")



