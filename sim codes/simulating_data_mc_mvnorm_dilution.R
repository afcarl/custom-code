require(multicore)
require(MASS)

transform <- function(vector,targetMean,targetSd) {
	s1 <- sd(vector)
	m1 <- mean(vector)
	s2 <- targetSd
	m2 <- targetMean
	return(m2 + (vector-m1) * (s2/s1))
}

findMean <- function(n,mean,sd) {
	list_r <- NULL
	diff <- vector(mode="numeric",length=100)
	for (try in 1:100) {
		list_r[[try]] <- rnorm(n=n,mean=mean,sd=sd)
		diff[try] <- abs(diff(c(mean(list_r[[try]]),mean)))
	}
	return(list_r[[which(diff == min(diff))]])
}

error <- 0.01
read_counts <- 500
error_terms <- 0.1
diff <- 0.04
correlations_normals <- 0.5
correlations_tumours <- 0
dilution <- seq(0,1,0.1)
sim.params <- expand.grid(dilution,read_counts,correlations_normals,correlations_tumours,diff,error_terms)
sim.params <- as.data.frame(sim.params)
colnames(sim.params) <- c("dilution","read_counts","correlations_normals","correlations_tumours","difference_fc","error_terms")


sim.data <- local({
	f <- fifo(tempfile(), open="w+b", blocking=T)
	if (inherits(fork(), "masterProcess")) {
		# Child
		progress <- 0.0
		ptm <- proc.time()[3]
		while (progress < 1 && !isIncomplete(f)) {
			msg <- readBin(f, "double")
			progress <- progress + as.numeric(msg)
			if (round(100*progress,digits=2)%%0.25==0) {
				eta <- ((proc.time()[3] - ptm) / progress) - (proc.time()[3] - ptm)
				cat(sprintf("Progress: %.2f%% ETA: %.2fmin\n", progress * 100, eta/60))
			}
		} 
		exit()
	}
	numJobs <- nrow(sim.params)*1000
	result <- mclapply(1:numJobs, mc.cores=47, mc.cleanup = TRUE, mc.preschedule = TRUE, function(i) {
		current_param <- ceiling(i / 1000)
		dilution <- sim.params[current_param,1]
		read_count <- sim.params[current_param,2]
		correlation_normals <- sim.params[current_param,3]
		correlation_tumours <- sim.params[current_param,4]
		diff <- sim.params[current_param,5]
		error_term <- sim.params[current_param,6]
		
		Sigma <- matrix(c(3,1.5,-1.5,1.5,3,-3,-1.5,-3,3),3,3,byrow=TRUE)
		temp <- mvrnorm(n=100, c(0,0,0), Sigma)
		temp_counts <- temp[,1]
		temp_gb <- temp[,2]
		temp_gb <- transform(vector=temp_gb,targetMean=0,targetSd=2*error_term)
		temp_pr <- temp[,3]
		temp_pr <- transform(vector=temp_pr,targetMean=0,targetSd=2*error_term)
		temp_counts <- round(transform(vector=temp_counts,targetMean=read_count,targetSd=error_term*read_count),digits=0)
		temp_an <- cbind(temp_counts,temp_gb,temp_pr)
		temp_pr_cpg <- matrix(ncol=11,nrow=nrow(temp_an))
		temp_gb_cpg <- matrix(ncol=11,nrow=nrow(temp_an))
		gb_cpg_means <- findMean(n=11,mean=0,sd=1.5*0.14)
		pr_cpg_means <- findMean(n=11,mean=0,sd=1.5*0.14)
		for (cpg in 1:11) {
			temp_gb_cpg[,cpg] <- transform(vector=temp_an[,2],targetMean=gb_cpg_means[cpg],targetSd=1.5*0.14)
			temp_pr_cpg[,cpg] <- transform(vector=temp_an[,3],targetMean=pr_cpg_means[cpg],targetSd=1.5*0.14)
		}
		temp_an <- cbind(temp_an,temp_gb_cpg,temp_pr_cpg)
		colnames(temp_an) <- c("count","gb","pr",paste("gb",seq(1,11,1),sep=""),paste("pr",seq(1,11,1),sep=""))
		
		# True Ts
		temp_t <- matrix(ncol=25,nrow=100,data=rep(0,25*100))
		n = round(100*(1-dilution))
		if (n > 0) {
			Sigma <- matrix(c(3,1.5-3*correlation_normals,-1.5+3*correlation_normals,1.5-3*correlation_normals,3,-3,-1.5+3*correlation_normals,-3,3),3,3,byrow=TRUE)
			temp <- mvrnorm(n=n, c(0,0,0), Sigma)
			temp_counts <- temp[,1]
			temp_gb <- temp[,2]
			temp_gb <- transform(vector=temp_gb,targetMean=diff*2.5,targetSd=2*error_term)
			temp_pr <- temp[,3]
			temp_pr <- transform(vector=temp_pr,targetMean=-diff*2.5,targetSd=2*error_term)
			temp_counts <- transform(vector=temp_counts,targetMean=read_count*(1+diff),targetSd=error_term*read_count*(1+diff))
			temp_t[1:n,1:3] <- cbind(temp_counts,temp_gb,temp_pr) 
			temp_t[1:n,1] <- round(temp_t[1:n,1],digits=0)
			temp_pr_cpg <- matrix(ncol=11,nrow=nrow(temp_t[1:n,]))
			temp_gb_cpg <- matrix(ncol=11,nrow=nrow(temp_t[1:n,]))
			for (cpg in 1:11) {
				temp_gb_cpg[,cpg] <- transform(vector=temp_t[1:n,2],targetMean=gb_cpg_means[cpg]+diff*2.5,targetSd=1.5*0.14)
				temp_pr_cpg[,cpg] <- transform(vector=temp_t[1:n,3],targetMean=pr_cpg_means[cpg]-diff*2.5,targetSd=1.5*0.14)
			}
			temp_t[1:n,4:25] <- cbind(temp_gb_cpg,temp_pr_cpg)
		}
		# Fake Ts
		n = 100-n
		if (n > 0) {
			Sigma <- matrix(c(3,1.5,-1.5,1.5,3,-3,-1.5,-3,3),3,3,byrow=TRUE)
			temp <- mvrnorm(n=n, c(0,0,0), Sigma)
			temp_counts <- temp[,1]
			temp_gb <- temp[,2]
			temp_gb <- transform(vector=temp_gb,targetMean=0,targetSd=2*error_term)
			temp_pr <- temp[,3]
			temp_pr <- transform(vector=temp_pr,targetMean=0,targetSd=2*error_term)
			temp_counts <- round(transform(vector=temp_counts,targetMean=read_count,targetSd=error_term*read_count),digits=0)
			temp_fake <- cbind(temp_counts,temp_gb,temp_pr)
			temp_pr_cpg <- matrix(ncol=11,nrow=nrow(temp_fake))
			temp_gb_cpg <- matrix(ncol=11,nrow=nrow(temp_fake))
			gb_cpg_means <- findMean(n=11,mean=0,sd=1.5*0.14)
			pr_cpg_means <- findMean(n=11,mean=0,sd=1.5*0.14)
			for (cpg in 1:11) {
				temp_gb_cpg[,cpg] <- transform(vector=temp_fake[,2],targetMean=gb_cpg_means[cpg],targetSd=1.5*0.14)
				temp_pr_cpg[,cpg] <- transform(vector=temp_fake[,3],targetMean=pr_cpg_means[cpg],targetSd=1.5*0.14)
			}
			temp_fake <- cbind(temp_fake,temp_gb_cpg,temp_pr_cpg)
			temp_t[(100-n+1):100,] <- temp_fake
		}
		colnames(temp_t) <- c("count","gb","pr",paste("gb",seq(1,11,1),sep=""),paste("pr",seq(1,11,1),sep=""))
		writeBin(1/numJobs, f)
		list(temp_an,temp_t)
		})
	close(f)
	cat("Done\n")
	result
})