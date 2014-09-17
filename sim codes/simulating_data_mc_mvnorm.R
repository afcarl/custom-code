require(parallel)
require(MASS)

transform <- function(vector,targetMean,targetSd) {
	s1 <- sd(vector)
	m1 <- mean(vector)
	s2 <- targetSd
	m2 <- targetMean
	return(m2 + (vector-m1) * (s2/s1))
}


error <- 0.05
read_counts <- 500
error_terms <- 0.14
diff <- 0.005
correlations_normals <- seq(0,1,0.05)
sim.params <- expand.grid(read_counts,correlations_normals,diff,error_terms)
sim.params <- cbind(sim.params,correlations_normals)
sim.params <- sim.params[,c(1,2,5,3,4)]
sim.params <- as.data.frame(sim.params)
colnames(sim.params) <- c("read_counts","correlations_normals","correlations_tumours","difference_fc","error_terms")
numJobs <- nrow(sim.params)*1000
sim.data <- mclapply(1:numJobs, mc.cores=47, mc.cleanup = TRUE, mc.preschedule = TRUE, function(i) {
	current_param <- ceiling(i / 1000)
	read_count <- sim.params[current_param,1]
	correlation_normals <- sim.params[current_param,2]
	correlation_tumours <- sim.params[current_param,3]
	diff <- sim.params[current_param,4]
	error_term <- sim.params[current_param,5]
	
	Sigma <- matrix(c(3,1.5,-1.5,1.5,3,-3,-1.5,-3,3),3,3,byrow=TRUE)
	temp <- mvrnorm(n=100, c(0,0,0), Sigma)
	temp_counts <- temp[,1]
	temp_gb <- temp[,2]
	temp_gb <- transform(vector=temp_gb,targetMean=0,targetSd=2*error_term)
	temp_pr <- temp[,3]
	temp_pr <- transform(vector=temp_pr,targetMean=0,targetSd=2*error_term)
	temp_counts <- trunc(transform(vector=temp_counts,targetMean=read_count,targetSd=0.1*read_count)+ rnorm(100,0,4))
	temp_an <- cbind(temp_counts,temp_gb,temp_pr)
	temp_pr_cpg <- matrix(ncol=11,nrow=nrow(temp_an))
	temp_gb_cpg <- matrix(ncol=11,nrow=nrow(temp_an))
	
	temp_pr_cpg <- matrix(ncol=11,nrow=nrow(temp_an))
	temp_gb_cpg <- matrix(ncol=11,nrow=nrow(temp_an))
	for (an in 1:nrow(temp_an)) {
		temp_pr_cpg[an,] <- rnorm(n=11,mean=temp_an[an,3],sd=2*error_term) + rnorm(11,0,error)
		temp_gb_cpg[an,] <- rnorm(n=11,mean=temp_an[an,2],sd=2*error_term) + rnorm(11,0,error)
	}
	
	temp_an <- cbind(temp_an,temp_gb_cpg,temp_pr_cpg)
	colnames(temp_an) <- c("count","gb","pr",paste("gb",seq(1,11,1),sep=""),paste("pr",seq(1,11,1),sep=""))
	
	# True Ts
	n = 100
	Sigma <- matrix(c(3,1.5-3*correlation_normals,-1.5+3*correlation_normals,1.5-3*correlation_normals,3,-3,-1.5+3*correlation_normals,-3,3),3,3,byrow=TRUE)
	temp <- mvrnorm(n=n, c(0,0,0), Sigma)
	temp_counts <- temp[,1]
	temp_gb <- temp[,2]
	temp_gb <- transform(vector=temp_gb,targetMean=diff*2.5,targetSd=2*error_term)
	temp_pr <- temp[,3]
	temp_pr <- transform(vector=temp_pr,targetMean=-diff*2.5,targetSd=2*error_term)
	temp_counts <- trunc(transform(vector=temp_counts,targetMean=read_count*(1+diff),targetSd=0.125*read_count*(1+diff))+ rnorm(n,0,8))
	temp_t <- cbind(temp_counts,temp_gb,temp_pr)
	temp_pr_cpg <- matrix(ncol=11,nrow=nrow(temp_t))
	temp_gb_cpg <- matrix(ncol=11,nrow=nrow(temp_t))
	for (t in 1:nrow(temp_t)) {
		temp_pr_cpg[t,] <- rnorm(n=11,mean=temp_t[t,3],sd=2.5*error_term) + rnorm(11,0,error)
		temp_gb_cpg[t,] <- rnorm(n=11,mean=temp_t[t,2],sd=2.5*error_term) + rnorm(11,0,error)
	}
	temp_t <- cbind(temp_t,temp_gb_cpg,temp_pr_cpg)
	colnames(temp_t) <- c("count","gb","pr",paste("gb",seq(1,11,1),sep=""),paste("pr",seq(1,11,1),sep=""))
	list(temp_an,temp_t)
	})
