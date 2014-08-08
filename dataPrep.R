#start here
#updated list of suitable genes is located in variable workingList_BRCA4
#load("./essentials_FG.RData")

IDs_promoter_all <- unique(c(unlist(dataf_BRCA_logit_PC@TSS1500Ind$SID),unlist(dataf_BRCA_logit_PC@TSS200Ind$SID),unlist(dataf_BRCA_logit_PC@UTR5Ind$SID),unlist(dataf_BRCA_logit_PC@EXON1Ind$SID)))
d_promoter <- density(dataf_BRCA_logit_PC@bmatrix[IDs_promoter_all,ANs2],bw=0.14,from=-7,to=7,n=1401)
d_promoter$y <- d_promoter$y/sum(d_promoter$y)
IDs_body_all <- unique(c(unlist(dataf_BRCA_logit_PC@GENEBODYInd$SID),unlist(dataf_BRCA_logit_PC@UTR3Ind$SID)))
d_body <- density(dataf_BRCA_logit_PC@bmatrix[IDs_body_all,ANs2],bw=0.14,from=-7,to=7,n=1401)
d_body$y <- d_body$y/sum(d_body$y)

#promoter_CpGs <- paste(promoterVars,".likelihood",sep="")
#geneBody_CpGs <- paste(geneBodyVars,".likelihood",sep="")
#require(Brobdingnag)
#length(workingList_BRCA4)
list_out <- NULL
for (i in 1:10)
{
  system(command=paste('mkdir',i,sep=" "))
  system(command=paste('mkdir ./',i,'/AN_model',sep=""))
  system(command=paste('mkdir ./',i,'/T_model',sep=""))
  IDs_promoter <- unique(c(eval(parse(text = paste('dataf_BRCA_logit_PC@TSS1500Ind$SID$','"',workingList_BRCA4[i],'"',sep=""))),eval(parse(text = paste('dataf_BRCA_logit_PC@TSS200Ind$SID$','"',workingList_BRCA4[i],'"',sep=""))), eval(parse(text = paste('dataf_BRCA_logit_PC@UTR5Ind$SID$','"',workingList_BRCA4[i],'"',sep="")))),eval(parse(text = paste('dataf_BRCA_logit_PC@EXON1Ind$SID$','"',workingList_BRCA4[i],'"',sep=""))))
  IDs_body <- unique(c(eval(parse(text = paste('dataf_BRCA_logit_PC@GENEBODYInd$SID$','"',workingList_BRCA4[i],'"',sep=""))), eval(parse(text = paste('dataf_BRCA_logit_PC@UTR3Ind$SID$','"',workingList_BRCA4[i],'"',sep="")))))
  # generate "missing" Var data
  ncol = length(IDs_promoter)+length(IDs_body)+1
  tempVar <- matrix(rep(".",75*ncol),nrow=75,ncol=ncol)
  colnames(tempVar) <- c("NAME:\tEXPR",promoterVars[1:length(IDs_promoter)],geneBodyVars[1:length(IDs_body)])
  rownames(tempVar) <- ANs2
  eval(parse(text = paste('write.table(', paste('tempVar,file = "./',i,'/AN_model/AN_VarData.tab",row.names=TRUE,col.names=TRUE,quote=FALSE,sep="\t",append=FALSE)', sep = ""))))  
  eval(parse(text = paste('write.table(', paste('tempVar,file = "./',i,'/T_model/AN_VarData.tab",row.names=TRUE,col.names=TRUE,quote=FALSE,sep="\t",append=FALSE)', sep = ""))))  
  rownames(tempVar) <- Ts2
  eval(parse(text = paste('write.table(', paste('tempVar,file = "./',i,'/AN_model/T_VarData.tab",row.names=TRUE,col.names=TRUE,quote=FALSE,sep="\t",append=FALSE)', sep = ""))))
  eval(parse(text = paste('write.table(', paste('tempVar,file = "./',i,'/T_model/T_VarData.tab",row.names=TRUE,col.names=TRUE,quote=FALSE,sep="\t",append=FALSE)', sep = ""))))
  ########################################################
  # dynamic generation of model specification files here #
  ########################################################
  stateMaps <- file(paste("./",i,"/stateMaps.txt",sep=""),"w")
  cat("NAME:\tfiveState\nSYMBOLS:\t1 2 3 4 5\nMETA_SYMBOLS:\t.=1 2 3 4 5; *=1 2 3 4 5;\n",file=stateMaps)
  close(stateMaps)
  
  variables <- file(paste("./",i,"/variables.txt",sep=""),"w")
  cat("STATE_MAP_NAME:\tfiveState\nVAR_NAMES:\tEXPR",promoterVars[1:length(IDs_promoter)],geneBodyVars[1:length(IDs_body)],sep=" ",file=variables)
  cat("\n",file=variables)
  close(variables)
  
  promoterPots <- paste("\nNAME:\t\tpot_",promoterVars,"\nTYPE:\t\trowNorm\nPOT_MAT:\t[1,5]((1,1,1,1,1))\nPC_MAT:\t\t[1,5]((0.1,0.1,0.1,0.1,0.1))\n",sep="")
  geneBodyPots <- paste("\nNAME:\t\tpot_",geneBodyVars,"\nTYPE:\t\trowNorm\nPOT_MAT:\t[1,5]((1,1,1,1,1))\nPC_MAT:\t\t[1,5]((0.1,0.1,0.1,0.1,0.1))\n",sep="")
  
  potentials <- file(paste("./",i,"/factorPotentials.txt",sep=""),"w")
  cat("NAME:\t\tpot_EXPR.CpG_P\nTYPE:\t\trowNorm\nPOT_MAT:\t[5,5]((1,1,1,1,1),\n\t\t\t(1,1,1,1,1),\n\t\t\t(1,1,1,1,1),\n\t\t\t(1,1,1,1,1),\n\t\t\t(1,1,1,1,1))\nPC_MAT:\t\t[5,5]((0.1,0.1,0.1,0.1,0.1),\n\t\t\t(0.1,0.1,0.1,0.1,0.1),\n\t\t\t(0.1,0.1,0.1,0.1,0.1),\n\t\t\t(0.1,0.1,0.1,0.1,0.1),\n\t\t\t(0.1,0.1,0.1,0.1,0.1))\n\nNAME:\t\tpot_EXPR.CpG_GB\nTYPE:\t\trowNorm\nPOT_MAT:\t[5,5]((1,1,1,1,1),\n\t\t\t(1,1,1,1,1),\n\t\t\t(1,1,1,1,1),\n\t\t\t(1,1,1,1,1),\n\t\t\t(1,1,1,1,1))\nPC_MAT:\t\t[5,5]((0.1,0.1,0.1,0.1,0.1),\n\t\t\t(0.1,0.1,0.1,0.1,0.1),\n\t\t\t(0.1,0.1,0.1,0.1,0.1),\n\t\t\t(0.1,0.1,0.1,0.1,0.1),\n\t\t\t(0.1,0.1,0.1,0.1,0.1))\n",
      promoterPots[1:length(IDs_promoter)],geneBodyPots[1:length(IDs_body)],sep="",file=potentials)
  cat("\nNAME:\t\tpot_EXPR.likelihood\nTYPE:\t\trowNorm\nPOT_MAT:\t[1,5]((1,1,1,1,1))\nPC_MAT:\t\t[1,5]((0.1,0.1,0.1,0.1,0.1))\n",file=potentials)
  close(potentials)

  factorGraph <- file(paste("./",i,"/factorGraph.txt",sep=""),"w")
  cat("NAME:\tEXPR.likelihood\nNB1:\tEXPR\nPOT:\tpot_EXPR.likelihood\n",file=factorGraph)
  cat(paste("\nNAME:\tEXPR.",promoterVars[1:length(IDs_promoter)],"\nNB1:\tEXPR\nNB2:\t",promoterVars[1:length(IDs_promoter)],"\nPOT:\tpot_EXPR.CpG_P\n",sep="",collapse=""),file=factorGraph)
  cat(paste("\nNAME:\tEXPR.",geneBodyVars[1:length(IDs_body)],"\nNB1:\tEXPR\nNB2:\t",geneBodyVars[1:length(IDs_body)],"\nPOT:\tpot_EXPR.CpG_GB\n",sep="",collapse=""),file=factorGraph)
  cat(paste("\nNAME:\t",promoterVars[1:length(IDs_promoter)],".likelihood\nNB1:\t",promoterVars[1:length(IDs_promoter)],"\nPOT:\tpot_",promoterVars[1:length(IDs_promoter)],"\n",sep="",collapse=""),file=factorGraph)
  cat(paste("\nNAME:\t",geneBodyVars[1:length(IDs_body)],".likelihood\nNB1:\t",geneBodyVars[1:length(IDs_body)],"\nPOT:\tpot_",geneBodyVars[1:length(IDs_body)],"\n",sep="",collapse=""),file=factorGraph)
  close(factorGraph)  
  ########################################################
  
  nsim <- 1 # partial expression evidence will sum up to this value
	percentiles <- c(0.2,0.4,0.6,0.8) # according to which the binning scheme will be defined
  
  ##########################################################################
  ############################## T model ###################################
  # 15-fold cross-validation begins here to obtain ANs likelihoods under AN model
  ANs_AN_likelihoods <- NULL
  for (fold in 1:15){
    current_ANs <- ANs2[-seq(from=fold*5-4,to=fold*5,by=1)]
    current_which_ANs <- which_ANs[-seq(from=fold*5-4,to=fold*5,by=1)]
    
    tempVar <- matrix(rep(".",70*ncol),nrow=70,ncol=ncol)
    colnames(tempVar) <- c("NAME:\tEXPR",promoterVars[1:length(IDs_promoter)],geneBodyVars[1:length(IDs_body)])
    rownames(tempVar) <- current_ANs
    eval(parse(text = paste('write.table(', paste('tempVar,file = "./',i,'/AN_model/',fold,'_AN_VarData.tab",row.names=TRUE,col.names=TRUE,quote=FALSE,sep="\t",append=FALSE)', sep = ""))))  

    # gene body  
    #f_body <- ecdf(dataf_BRCA_logit_PC@bmatrix[IDs_body,current_ANs])
    #breaksBODY <- quantile(f_body,probs=seq(0,1,by=0.1)) # 10 break points
    mean <- mean(dataf_BRCA_logit_PC@bmatrix[IDs_body,current_ANs])
    sd <- sd(dataf_BRCA_logit_PC@bmatrix[IDs_body,current_ANs])
    breaksBODY <- round(sort(c(-7.01,qnorm(percentiles,mean=mean,sd=sd),7.01)),digits=2)
    
    # promoter	
    #f_promoter <- ecdf(dataf_BRCA_logit_PC@bmatrix[IDs_promoter,current_ANs])
    #breaksPROMOTER <- quantile(f_promoter,probs=seq(0,1,by=0.1)) # 20 break points
    mean <- mean(dataf_BRCA_logit_PC@bmatrix[IDs_promoter,current_ANs])
    sd <- sd(dataf_BRCA_logit_PC@bmatrix[IDs_promoter,current_ANs])
    breaksPROMOTER <- round(sort(c(-7.01,qnorm(percentiles,mean=mean,sd=sd),7.01)),digits=2)
    # expression
    lambda <- counts_BRCA_plusOne[workingList_BRCA4[i],current_ANs[1]]
    X <- seq(round(max(lambda-(4*lambda*lambda^(-1/2)),1)),round(lambda+(4*lambda*lambda^(-1/2))))
    list <- NULL
    list[[1]] <- cbind(X/factors_ls[current_which_ANs[1]],dpois(X,lambda=lambda)*factors_ls[current_which_ANs[1]])
    for (j in 2:length(current_ANs)) {
      lambda <- counts_BRCA_plusOne[workingList_BRCA4[i],current_ANs[j]]
      X <- seq(round(max(lambda-(4*lambda*lambda^(-1/2)),1)),round(lambda+(4*lambda*lambda^(-1/2))))
      current <- factors_ls[current_which_ANs[j]]
      list[[j]] <- cbind(X/current,dpois(X,lambda=lambda)*current)
    }
    
    mean <- mean(cpm_BRCA_plusOne[workingList_BRCA4[i],current_ANs])
    sd <- sd(cpm_BRCA_plusOne[workingList_BRCA4[i],current_ANs])
    breaksEXPRESSION <- sort(c(0,qnorm(percentiles,mean=mean,sd=sd),1000000))
    
    # generate data for ANs
    tempS <- matrix(ncol=ncol)
    for (current_sample in 1:length(current_ANs)) {
      
      all_labels <- as.character(seq(1,5,1))
      
      # expression
      AN <- list[[current_sample]]
      AN <- as.data.frame(AN)
      colnames(AN) <- c("cpm","density")
      AN$interval <- findInterval(AN[,1],breaksEXPRESSION)
      AN <- aggregate(cbind(cpm,density) ~ interval, data=AN, sum) #sum up the repetitive intervals' densities
      AN$interval <- as.factor(AN$interval)
      labels <- colnames(t(summary(AN$interval)))
      frequencies <- AN$density
      frequencies <- AN$density * nsim/sum(AN$density)
      frequencies_expr <- rep(0,5) # EXPR.likelihood
      for (label in 1:5) {if (!is.na(match(all_labels[label], labels))) {
        frequencies_expr[label] <- frequencies[match(all_labels[label], labels)]
      }}  	
      
      # gene body
      cpg_list_gb <- NULL # CpG_GB_X.likelihood
      for (gb_cpg in 1:length(IDs_body)) {
        x <- dataf_BRCA_logit_PC@bmatrix[IDs_body[gb_cpg],current_ANs[current_sample]]
        if (x < -7) x <- -7
        dens <- density(x,bw=0.14,from=-7,to=7,n=1401)
        dens$y <- dens$y/sum(dens$y)
        dens$y <- dens$y * d_body$y # P(M|d)*P(d)
        dens <- as.data.frame(cbind(dens$x,dens$y))
        colnames(dens) <- c("M","density")
        dens$interval <- findInterval(dens[,1],breaksBODY)
        dens <- aggregate(cbind(M,density) ~ interval, data=dens, sum)
        for (bin in 2:length(breaksBODY)) { # (P(M|d)*P(d))/P(M)
          dens[bin-1,3] <- dens[bin-1,3] / sum(d_body$y[(intersect(which(d_body$x >= breaksBODY[bin-1]),which(d_body$x < breaksBODY[bin])))])
        }
        dens$interval <- as.factor(dens$interval)
        labels <- colnames(t(summary(dens$interval)))
        frequencies <- dens$density
        frequencies_gb <- rep(0,5)
        for (label in 1:5) {
          if (!is.na(match(all_labels[label], labels)))
            frequencies_gb[label] <- frequencies[match(all_labels[label], labels)]
        }
        cpg_list_gb[[gb_cpg]] <- frequencies_gb
      }
      
      # promoter
      cpg_list_pr <- NULL # CpG_P_X.likelihood
      for (pr_cpg in 1:length(IDs_promoter)) {
        x <- dataf_BRCA_logit_PC@bmatrix[IDs_promoter[pr_cpg],current_ANs[current_sample]]
        if (x < -7) x <- -7
        dens <- density(x,bw=0.14,from=-7,to=7,n=1401)
        dens$y <- dens$y/sum(dens$y)
        dens$y <- dens$y * d_promoter$y # P(M|d)*P(d)
        dens <- as.data.frame(cbind(dens$x,dens$y))
        colnames(dens) <- c("M","density")
        dens$interval <- findInterval(dens[,1],breaksPROMOTER)
        dens <- aggregate(cbind(M,density) ~ interval, data=dens, sum)
        for (bin in 2:length(breaksPROMOTER)) { # (P(M|d)*P(d))/P(M)
          dens[bin-1,3] <- dens[bin-1,3] / sum(d_promoter$y[(intersect(which(d_promoter$x >= breaksPROMOTER[bin-1]),which(d_promoter$x < breaksPROMOTER[bin])))])
        }
        dens$interval <- as.factor(dens$interval)
        labels <- colnames(t(summary(dens$interval)))
        frequencies <- dens$density
        frequencies_pr <- rep(0,5)
        for (label in 1:5) {
          if (!is.na(match(all_labels[label], labels)))
            frequencies_pr[label] <- frequencies[match(all_labels[label], labels)]
        }
        cpg_list_pr[[pr_cpg]] <- frequencies_pr
      }
      
      tempS_formated <- matrix()
      tempS_formated[1,1] <- paste('[1,5]((',paste(frequencies_expr,sep="",collapse=","),'))',sep="")
      for (element in 1:length(cpg_list_pr)) {
        tempS_formated <- cbind(tempS_formated,paste('[1,5]((',paste(cpg_list_pr[[element]],sep="",collapse=","),'))',sep=""))
      }
      for (element in 1:length(cpg_list_gb)) {
        tempS_formated <- cbind(tempS_formated,paste('[1,5]((',paste(cpg_list_gb[[element]],sep="",collapse=","),'))',sep=""))
      }
      tempS <- rbind(tempS,tempS_formated)
    }
    tempFac <- tempS[2:71,]
    colnames(tempFac) <- c("NAME:\tEXPR.likelihood",promoter_CpGs[1:length(IDs_promoter)],geneBody_CpGs[1:length(IDs_body)])
    rownames(tempFac) <- current_ANs  
    eval(parse(text = paste('write.table(', paste('tempFac,file ="./',i,'/AN_model/',fold,'_AN_FacData.tab",row.names=TRUE,col.names=TRUE,quote=FALSE,sep="\t",append=FALSE)', sep = ""))))
    
    # build and query the current x-val fold model
    system(command=paste('mkdir ./',i,'/AN_model/',fold,'/',sep=""))
    system(command=paste('../phy/src/dfgTrain --emTrain --maxIter 1000 --dfgSpecPrefix=./',i,'/ --outSpecPrefix=./',i,'/AN_model/',fold,'/ ./',i,'/AN_model/',fold,'_AN_VarData.tab ./',i,'/AN_model/',fold,'_AN_FacData.tab',sep=""))
    l5o_ANs <- ANs2[seq(from=fold*5-4,to=fold*5,by=1)]
    l5o_which_ANs <- which_ANs[seq(from=fold*5-4,to=fold*5,by=1)]
    
    # expression
    lambda <- counts_BRCA_plusOne[workingList_BRCA4[i],l5o_ANs[1]]
    X <- seq(round(max(lambda-(4*lambda*lambda^(-1/2)),1)),round(lambda+(4*lambda*lambda^(-1/2))))
    list <- NULL
    list[[1]] <- cbind(X/factors_ls[l5o_which_ANs[1]],dpois(X,lambda=lambda)*factors_ls[l5o_which_ANs[1]])
    for (j in 2:length(l5o_ANs)) {
      lambda <- counts_BRCA_plusOne[workingList_BRCA4[i],l5o_ANs[j]]
      X <- seq(round(max(lambda-(4*lambda*lambda^(-1/2)),1)),round(lambda+(4*lambda*lambda^(-1/2))))
      current <- factors_ls[l5o_which_ANs[j]]
      list[[j]] <- cbind(X/current,dpois(X,lambda=lambda)*current)
    }
    
    tempS <- matrix(ncol=ncol)
    for (current_sample in 1:length(l5o_ANs)) {     
      all_labels <- as.character(seq(1,5,1))
      
      # expression
      AN <- list[[current_sample]]
      AN <- as.data.frame(AN)
      colnames(AN) <- c("cpm","density")
      AN$interval <- findInterval(AN[,1],breaksEXPRESSION)
      AN <- aggregate(cbind(cpm,density) ~ interval, data=AN, sum) #sum up the repetitive intervals' densities
      AN$interval <- as.factor(AN$interval)
      labels <- colnames(t(summary(AN$interval)))
      frequencies <- AN$density
      frequencies <- AN$density * nsim/sum(AN$density)
      frequencies_expr <- rep(0,5) # EXPR.likelihood
      for (label in 1:5) {if (!is.na(match(all_labels[label], labels))) {
        frequencies_expr[label] <- frequencies[match(all_labels[label], labels)]
      }}    
      
      # gene body
      cpg_list_gb <- NULL # CpG_GB_X.likelihood
      for (gb_cpg in 1:length(IDs_body)) {
        x <- dataf_BRCA_logit_PC@bmatrix[IDs_body[gb_cpg],l5o_ANs[current_sample]]
        if (x < -7) x <- -7
        dens <- density(x,bw=0.14,from=-7,to=7,n=1401)
        dens$y <- dens$y/sum(dens$y)
        dens$y <- dens$y * d_body$y # P(M|d)*P(d)
        dens <- as.data.frame(cbind(dens$x,dens$y))
        colnames(dens) <- c("M","density")
        dens$interval <- findInterval(dens[,1],breaksBODY)
        dens <- aggregate(cbind(M,density) ~ interval, data=dens, sum)
        for (bin in 2:length(breaksBODY)) { # (P(M|d)*P(d))/P(M)
          dens[bin-1,3] <- dens[bin-1,3] / sum(d_body$y[(intersect(which(d_body$x >= breaksBODY[bin-1]),which(d_body$x < breaksBODY[bin])))])
        }
        dens$interval <- as.factor(dens$interval)
        labels <- colnames(t(summary(dens$interval)))
        frequencies <- dens$density
        frequencies_gb <- rep(0,5)
        for (label in 1:5) {
          if (!is.na(match(all_labels[label], labels)))
            frequencies_gb[label] <- frequencies[match(all_labels[label], labels)]
        }
        cpg_list_gb[[gb_cpg]] <- frequencies_gb
      }
      
      # promoter
      cpg_list_pr <- NULL # CpG_P_X.likelihood
      for (pr_cpg in 1:length(IDs_promoter)) {
        dens <- density(dataf_BRCA_logit_PC@bmatrix[IDs_promoter[pr_cpg],l5o_ANs[current_sample]],bw=0.14,from=-7,to=7,n=1401)
        dens$y <- dens$y/sum(dens$y)
        dens$y <- dens$y * d_promoter$y # P(M|d)*P(d)
        dens <- as.data.frame(cbind(dens$x,dens$y))
        colnames(dens) <- c("M","density")
        dens$interval <- findInterval(dens[,1],breaksPROMOTER)
        dens <- aggregate(cbind(M,density) ~ interval, data=dens, sum)
        for (bin in 2:length(breaksPROMOTER)) { # (P(M|d)*P(d))/P(M)
          dens[bin-1,3] <- dens[bin-1,3] / sum(d_promoter$y[(intersect(which(d_promoter$x >= breaksPROMOTER[bin-1]),which(d_promoter$x < breaksPROMOTER[bin])))])
        }
        dens$interval <- as.factor(dens$interval)
        labels <- colnames(t(summary(dens$interval)))
        frequencies <- dens$density
        frequencies_pr <- rep(0,5)
        for (label in 1:5) {
          if (!is.na(match(all_labels[label], labels)))
            frequencies_pr[label] <- frequencies[match(all_labels[label], labels)]
        }
        cpg_list_pr[[pr_cpg]] <- frequencies_pr
      }
      
      tempS_formated <- matrix()
      tempS_formated[1,1] <- paste('[1,5]((',paste(frequencies_expr,sep="",collapse=","),'))',sep="")
      for (element in 1:length(cpg_list_pr)) {
        tempS_formated <- cbind(tempS_formated,paste('[1,5]((',paste(cpg_list_pr[[element]],sep="",collapse=","),'))',sep=""))
      }
      for (element in 1:length(cpg_list_gb)) {
        tempS_formated <- cbind(tempS_formated,paste('[1,5]((',paste(cpg_list_gb[[element]],sep="",collapse=","),'))',sep=""))
      }
      tempS <- rbind(tempS,tempS_formated)
    }
    tempFac <- tempS[2:6,]
    colnames(tempFac) <- c("NAME:\tEXPR.likelihood",promoter_CpGs[1:length(IDs_promoter)],geneBody_CpGs[1:length(IDs_body)])
    rownames(tempFac) <- l5o_ANs  
    eval(parse(text = paste('write.table(', paste('tempFac,file ="./',i,'/AN_model/',fold,'_evalAN_FacData.tab",row.names=TRUE,col.names=TRUE,quote=FALSE,sep="\t",append=FALSE)', sep = ""))))
    tempVar <- matrix(rep(".",5*ncol),nrow=5,ncol=ncol)
    colnames(tempVar) <- c("NAME:\tEXPR",promoterVars[1:length(IDs_promoter)],geneBodyVars[1:length(IDs_body)])
    rownames(tempVar) <- l5o_ANs
    eval(parse(text = paste('write.table(', paste('tempVar,file = "./',i,'/AN_model/',fold,'_evalAN_VarData.tab",row.names=TRUE,col.names=TRUE,quote=FALSE,sep="\t",append=FALSE)', sep = ""))))  
    
    string<-system(intern=TRUE,command=paste('../phy/src/dfgEval --dfgSpecPrefix=./',i,'/AN_model/',fold,'/ -l -n - ./',i,'/AN_model/',fold,'_evalAN_VarData.tab ./',i,'/AN_model/',fold,'_evalAN_FacData.tab',sep=""))
    ANs_AN_likelihoods <- c(ANs_AN_likelihoods,string)
  }
  ANs_AN_likelihoods <- as.numeric(substring(ANs_AN_likelihoods[-seq(1,90,6)],17))
  
  # full AN model developed from here, to obtain likelihoods of Ts
  # gene body  
  mean <- mean(dataf_BRCA_logit_PC@bmatrix[IDs_body,ANs2])
  sd <- sd(dataf_BRCA_logit_PC@bmatrix[IDs_body,ANs2])
  breaksBODY <- round(sort(c(-7.01,qnorm(percentiles,mean=mean,sd=sd),7.01)),digits=2)
  
  # promoter  
  mean <- mean(dataf_BRCA_logit_PC@bmatrix[IDs_promoter,ANs2])
  sd <- sd(dataf_BRCA_logit_PC@bmatrix[IDs_promoter,ANs2])
  breaksPROMOTER <- round(sort(c(-7.01,qnorm(percentiles,mean=mean,sd=sd),7.01)),digits=2)
  # expression
  mean <- mean(cpm_BRCA_plusOne[workingList_BRCA4[i],ANs2])
  sd <- sd(cpm_BRCA_plusOne[workingList_BRCA4[i],ANs2])
  breaksEXPRESSION <- sort(c(0,qnorm(percentiles,mean=mean,sd=sd),1000000))
  
  # generate FacData for Ts
  lambda <- counts_BRCA_plusOne[workingList_BRCA4[i],Ts2[1]]
  X <- seq(round(max(lambda-(4*lambda*lambda^(-1/2)),1)),round(lambda+(4*lambda*lambda^(-1/2))))
  list <- NULL
  list[[1]] <- cbind(X/factors_ls[which_Ts[1]],dpois(X,lambda=lambda)*factors_ls[which_Ts[1]])
  for (j in 2:length(Ts2)) {
    lambda <- counts_BRCA_plusOne[workingList_BRCA4[i],Ts2[j]]
    X <- seq(round(max(lambda-(4*lambda*lambda^(-1/2)),1)),round(lambda+(4*lambda*lambda^(-1/2))))
    current <- factors_ls[which_Ts[j]]
    list[[j]] <- cbind(X/current,dpois(X,lambda=lambda)*current)
  } 
  tempS <- matrix(ncol=ncol)
  for (i_t in 1:length(Ts2)) {
    
    all_labels <- as.character(seq(1,5,1))
    
    # expression
    T <- list[[i_t]]
    T <- as.data.frame(T)
    colnames(T) <- c("cpm","density")
    T$interval <- findInterval(T[,1],breaksEXPRESSION)
    T <- aggregate(cbind(cpm,density) ~ interval, data=T, sum) #sum up the repetitive intervals' densities
    T$interval <- as.factor(T$interval)
    labels <- colnames(t(summary(T$interval)))
    frequencies <- T$density
    frequencies <- T$density * nsim/sum(T$density)
    frequencies_expr <- rep(0,5) # EXPR.likelihood
    for (label in 1:5) {if (!is.na(match(all_labels[label], labels))) {
      frequencies_expr[label] <- frequencies[match(all_labels[label], labels)]
    }}  	
    
    # gene body
    cpg_list_gb <- NULL # CpG_GB_X.likelihood
    for (gb_cpg in 1:length(IDs_body)) {
      x <- dataf_BRCA_logit_PC@bmatrix[IDs_body[gb_cpg],Ts2[i_t]]
      if (x < -7) x <- -7
      dens <- density(x,bw=0.14,from=-7,to=7,n=1401)
      dens$y <- dens$y/sum(dens$y)
      dens$y <- dens$y * d_body$y # P(M|d)*P(d)
      dens <- as.data.frame(cbind(dens$x,dens$y))
      colnames(dens) <- c("M","density")
      dens$interval <- findInterval(dens[,1],breaksBODY)
      dens <- aggregate(cbind(M,density) ~ interval, data=dens, sum)
      for (bin in 2:length(breaksBODY)) { # (P(M|d)*P(d))/P(M)
        dens[bin-1,3] <- dens[bin-1,3] / sum(d_body$y[(intersect(which(d_body$x >= breaksBODY[bin-1]),which(d_body$x < breaksBODY[bin])))])
      }
      dens$interval <- as.factor(dens$interval)
      labels <- colnames(t(summary(dens$interval)))
      frequencies <- dens$density
      frequencies_gb <- rep(0,5)
      for (label in 1:5) {
        if (!is.na(match(all_labels[label], labels)))
          frequencies_gb[label] <- frequencies[match(all_labels[label], labels)]
      }
      cpg_list_gb[[gb_cpg]] <- frequencies_gb
    }
    
    # promoter
    cpg_list_pr <- NULL # CpG_P_X.likelihood
    for (pr_cpg in 1:length(IDs_promoter)) {
      dens <- density(dataf_BRCA_logit_PC@bmatrix[IDs_promoter[pr_cpg],Ts2[i_t]],bw=0.14,from=-7,to=7,n=1401)
      dens$y <- dens$y/sum(dens$y)
      dens$y <- dens$y * d_promoter$y # P(M|d)*P(d)
      dens <- as.data.frame(cbind(dens$x,dens$y))
      colnames(dens) <- c("M","density")
      dens$interval <- findInterval(dens[,1],breaksPROMOTER)
      dens <- aggregate(cbind(M,density) ~ interval, data=dens, sum)
      for (bin in 2:length(breaksPROMOTER)) { # (P(M|d)*P(d))/P(M)
        dens[bin-1,3] <- dens[bin-1,3] / sum(d_promoter$y[(intersect(which(d_promoter$x >= breaksPROMOTER[bin-1]),which(d_promoter$x < breaksPROMOTER[bin])))])
      }
      dens$interval <- as.factor(dens$interval)
      labels <- colnames(t(summary(dens$interval)))
      frequencies <- dens$density
      frequencies_pr <- rep(0,5)
      for (label in 1:5) {
        if (!is.na(match(all_labels[label], labels)))
          frequencies_pr[label] <- frequencies[match(all_labels[label], labels)]
      }
      cpg_list_pr[[pr_cpg]] <- frequencies_pr
    }
    
    tempS_formated <- matrix()
    tempS_formated[1,1] <- paste('[1,5]((',paste(frequencies_expr,sep="",collapse=","),'))',sep="")
    for (element in 1:length(cpg_list_pr)) {
      tempS_formated <- cbind(tempS_formated,paste('[1,5]((',paste(cpg_list_pr[[element]],sep="",collapse=","),'))',sep=""))
    }
    for (element in 1:length(cpg_list_gb)) {
      tempS_formated <- cbind(tempS_formated,paste('[1,5]((',paste(cpg_list_gb[[element]],sep="",collapse=","),'))',sep=""))
    }
    tempS <- rbind(tempS,tempS_formated)
  }
  tempFac <- tempS[2:76,]
  colnames(tempFac) <- c("NAME:\tEXPR.likelihood",promoter_CpGs[1:length(IDs_promoter)],geneBody_CpGs[1:length(IDs_body)])
  rownames(tempFac) <- Ts2
  eval(parse(text = paste('write.table(', paste('tempFac,file = "./',i,'/AN_model/T_FacData.tab",row.names=TRUE,col.names=TRUE,quote=FALSE,sep="\t",append=FALSE)', sep = ""))))
  
  # generate FacData for full set of ANs
  lambda <- counts_BRCA_plusOne[workingList_BRCA4[i],ANs2[1]]
  X <- seq(round(max(lambda-(4*lambda*lambda^(-1/2)),1)),round(lambda+(4*lambda*lambda^(-1/2))))
  list <- NULL
  list[[1]] <- cbind(X/factors_ls[which_ANs[1]],dpois(X,lambda=lambda)*factors_ls[which_ANs[1]])
  for (j in 2:length(ANs2)) {
    lambda <- counts_BRCA_plusOne[workingList_BRCA4[i],ANs2[j]]
    X <- seq(round(max(lambda-(4*lambda*lambda^(-1/2)),1)),round(lambda+(4*lambda*lambda^(-1/2))))
    current <- factors_ls[which_ANs[j]]
    list[[j]] <- cbind(X/current,dpois(X,lambda=lambda)*current)
  }  
  tempS <- matrix(ncol=ncol)
  for (current_sample in 1:length(ANs2)) {
    
    all_labels <- as.character(seq(1,5,1))
    
    # expression
    AN <- list[[current_sample]]
    AN <- as.data.frame(AN)
    colnames(AN) <- c("cpm","density")
    AN$interval <- findInterval(AN[,1],breaksEXPRESSION)
    AN <- aggregate(cbind(cpm,density) ~ interval, data=AN, sum) #sum up the repetitive intervals' densities
    AN$interval <- as.factor(AN$interval)
    labels <- colnames(t(summary(AN$interval)))
    frequencies <- AN$density
    frequencies <- AN$density * nsim/sum(AN$density)
    frequencies_expr <- rep(0,5) # EXPR.likelihood
    for (label in 1:5) {if (!is.na(match(all_labels[label], labels))) {
      frequencies_expr[label] <- frequencies[match(all_labels[label], labels)]
    }}    
    
    # gene body
    cpg_list_gb <- NULL # CpG_GB_X.likelihood
    for (gb_cpg in 1:length(IDs_body)) {
      x <- dataf_BRCA_logit_PC@bmatrix[IDs_body[gb_cpg],ANs2[current_sample]]
      if (x < -7) x <- -7
      dens <- density(x,bw=0.14,from=-7,to=7,n=1401)
      dens$y <- dens$y/sum(dens$y)
      dens$y <- dens$y * d_body$y # P(M|d)*P(d)
      dens <- as.data.frame(cbind(dens$x,dens$y))
      colnames(dens) <- c("M","density")
      dens$interval <- findInterval(dens[,1],breaksBODY)
      dens <- aggregate(cbind(M,density) ~ interval, data=dens, sum)
      for (bin in 2:length(breaksBODY)) { # (P(M|d)*P(d))/P(M)
        dens[bin-1,3] <- dens[bin-1,3] / sum(d_body$y[(intersect(which(d_body$x >= breaksBODY[bin-1]),which(d_body$x < breaksBODY[bin])))])
      }
      dens$interval <- as.factor(dens$interval)
      labels <- colnames(t(summary(dens$interval)))
      frequencies <- dens$density
      frequencies_gb <- rep(0,5)
      for (label in 1:5) {
        if (!is.na(match(all_labels[label], labels)))
          frequencies_gb[label] <- frequencies[match(all_labels[label], labels)]
      }
      cpg_list_gb[[gb_cpg]] <- frequencies_gb
    }
    
    # promoter
    cpg_list_pr <- NULL # CpG_P_X.likelihood
    for (pr_cpg in 1:length(IDs_promoter)) {
      dens <- density(dataf_BRCA_logit_PC@bmatrix[IDs_promoter[pr_cpg],ANs2[current_sample]],bw=0.14,from=-7,to=7,n=1401)
      dens$y <- dens$y/sum(dens$y)
      dens$y <- dens$y * d_promoter$y # P(M|d)*P(d)
      dens <- as.data.frame(cbind(dens$x,dens$y))
      colnames(dens) <- c("M","density")
      dens$interval <- findInterval(dens[,1],breaksPROMOTER)
      dens <- aggregate(cbind(M,density) ~ interval, data=dens, sum)
      for (bin in 2:length(breaksPROMOTER)) { # (P(M|d)*P(d))/P(M)
        dens[bin-1,3] <- dens[bin-1,3] / sum(d_promoter$y[(intersect(which(d_promoter$x >= breaksPROMOTER[bin-1]),which(d_promoter$x < breaksPROMOTER[bin])))])
      }
      dens$interval <- as.factor(dens$interval)
      labels <- colnames(t(summary(dens$interval)))
      frequencies <- dens$density
      frequencies_pr <- rep(0,5)
      for (label in 1:5) {
        if (!is.na(match(all_labels[label], labels)))
          frequencies_pr[label] <- frequencies[match(all_labels[label], labels)]
      }
      cpg_list_pr[[pr_cpg]] <- frequencies_pr
    }
    
    tempS_formated <- matrix()
    tempS_formated[1,1] <- paste('[1,5]((',paste(frequencies_expr,sep="",collapse=","),'))',sep="")
    for (element in 1:length(cpg_list_pr)) {
      tempS_formated <- cbind(tempS_formated,paste('[1,5]((',paste(cpg_list_pr[[element]],sep="",collapse=","),'))',sep=""))
    }
    for (element in 1:length(cpg_list_gb)) {
      tempS_formated <- cbind(tempS_formated,paste('[1,5]((',paste(cpg_list_gb[[element]],sep="",collapse=","),'))',sep=""))
    }
    tempS <- rbind(tempS,tempS_formated)
  }
  tempFac <- tempS[2:76,]
  colnames(tempFac) <- c("NAME:\tEXPR.likelihood",promoter_CpGs[1:length(IDs_promoter)],geneBody_CpGs[1:length(IDs_body)])
  rownames(tempFac) <- ANs2  
  eval(parse(text = paste('write.table(', paste('tempFac,file ="./',i,'/AN_model/AN_FacData.tab",row.names=TRUE,col.names=TRUE,quote=FALSE,sep="\t",append=FALSE)', sep = ""))))
  
  # build and query the full model with T samples
  system(command=paste('mkdir ./',i,'/AN_model/','all/',sep=""))
  system(command=paste('../phy/src/dfgTrain --emTrain --maxIter 1000 --dfgSpecPrefix=./',i,'/ --outSpecPrefix=./',i,'/AN_model/all/ ./',i,'/AN_model/AN_VarData.tab ./',i,'/AN_model/AN_FacData.tab',sep=""))
  string<-system(intern=TRUE,command=paste('../phy/src/dfgEval --dfgSpecPrefix=./',i,'/AN_model/',fold,'/ -l -n - ./',i,'/AN_model/T_VarData.tab ./',i,'/AN_model/T_FacData.tab',sep=""))
  Ts_AN_likelihoods <- as.numeric(substring(string[-1],17))
  
  ##########################################################################
  ############################## T model ###################################
  # 15-fold cross-validation begins here to obtain Ts likelihoods under T model
  Ts_T_likelihoods <- NULL
  for (fold in 1:15){
    current_Ts <- Ts2[-seq(from=fold*5-4,to=fold*5,by=1)]
    current_which_Ts <- which_Ts[-seq(from=fold*5-4,to=fold*5,by=1)]
    
    tempVar <- matrix(rep(".",70*ncol),nrow=70,ncol=ncol)
    colnames(tempVar) <- c("NAME:\tEXPR",promoterVars[1:length(IDs_promoter)],geneBodyVars[1:length(IDs_body)])
    rownames(tempVar) <- current_Ts
    eval(parse(text = paste('write.table(', paste('tempVar,file = "./',i,'/T_model/',fold,'_T_VarData.tab",row.names=TRUE,col.names=TRUE,quote=FALSE,sep="\t",append=FALSE)', sep = ""))))  
    
    # gene body  
    #f_body <- ecdf(dataf_BRCA_logit_PC@bmatrix[IDs_body,current_Ts])
    #breaksBODY <- quantile(f_body,probs=seq(0,1,by=0.1)) # 10 break points
    mean <- mean(dataf_BRCA_logit_PC@bmatrix[IDs_body,current_Ts])
    sd <- sd(dataf_BRCA_logit_PC@bmatrix[IDs_body,current_Ts])
    breaksBODY <- round(sort(c(-7.01,qnorm(percentiles,mean=mean,sd=sd),7.01)),digits=2)
    
    # promoter  
    #f_promoter <- ecdf(dataf_BRCA_logit_PC@bmatrix[IDs_promoter,current_Ts])
    #breaksPROMOTER <- quantile(f_promoter,probs=seq(0,1,by=0.1)) # 20 break points
    mean <- mean(dataf_BRCA_logit_PC@bmatrix[IDs_promoter,current_Ts])
    sd <- sd(dataf_BRCA_logit_PC@bmatrix[IDs_promoter,current_Ts])
    breaksPROMOTER <- round(sort(c(-7.01,qnorm(percentiles,mean=mean,sd=sd),7.01)),digits=2)
    # expression
    lambda <- counts_BRCA_plusOne[workingList_BRCA4[i],current_Ts[1]]
    X <- seq(round(max(lambda-(4*lambda*lambda^(-1/2)),1)),round(lambda+(4*lambda*lambda^(-1/2))))
    list <- NULL
    list[[1]] <- cbind(X/factors_ls[current_which_Ts[1]],dpois(X,lambda=lambda)*factors_ls[current_which_Ts[1]])
    for (j in 2:length(current_Ts)) {
      lambda <- counts_BRCA_plusOne[workingList_BRCA4[i],current_Ts[j]]
      X <- seq(round(max(lambda-(4*lambda*lambda^(-1/2)),1)),round(lambda+(4*lambda*lambda^(-1/2))))
      current <- factors_ls[current_which_Ts[j]]
      list[[j]] <- cbind(X/current,dpois(X,lambda=lambda)*current)
    }
    
    mean <- mean(cpm_BRCA_plusOne[workingList_BRCA4[i],current_Ts])
    sd <- sd(cpm_BRCA_plusOne[workingList_BRCA4[i],current_Ts])
    breaksEXPRESSION <- sort(c(0,qnorm(percentiles,mean=mean,sd=sd),1000000))
    
    # generate data for Ts
    tempS <- matrix(ncol=ncol)
    for (current_sample in 1:length(current_Ts)) {
      
      all_labels <- as.character(seq(1,5,1))
      
      # expression
      AN <- list[[current_sample]]
      AN <- as.data.frame(AN)
      colnames(AN) <- c("cpm","density")
      AN$interval <- findInterval(AN[,1],breaksEXPRESSION)
      AN <- aggregate(cbind(cpm,density) ~ interval, data=AN, sum) #sum up the repetitive intervals' densities
      AN$interval <- as.factor(AN$interval)
      labels <- colnames(t(summary(AN$interval)))
      frequencies <- AN$density
      frequencies <- AN$density * nsim/sum(AN$density)
      frequencies_expr <- rep(0,5) # EXPR.likelihood
      for (label in 1:5) {if (!is.na(match(all_labels[label], labels))) {
        frequencies_expr[label] <- frequencies[match(all_labels[label], labels)]
      }}    
      
      # gene body
      cpg_list_gb <- NULL # CpG_GB_X.likelihood
      for (gb_cpg in 1:length(IDs_body)) {
        x <- dataf_BRCA_logit_PC@bmatrix[IDs_body[gb_cpg],current_Ts[current_sample]]
        if (x < -7) x <- -7
        dens <- density(x,bw=0.14,from=-7,to=7,n=1401)
        dens$y <- dens$y/sum(dens$y)
        dens$y <- dens$y * d_body$y # P(M|d)*P(d)
        dens <- as.data.frame(cbind(dens$x,dens$y))
        colnames(dens) <- c("M","density")
        dens$interval <- findInterval(dens[,1],breaksBODY)
        dens <- aggregate(cbind(M,density) ~ interval, data=dens, sum)
        for (bin in 2:length(breaksBODY)) { # (P(M|d)*P(d))/P(M)
          dens[bin-1,3] <- dens[bin-1,3] / sum(d_body$y[(intersect(which(d_body$x >= breaksBODY[bin-1]),which(d_body$x < breaksBODY[bin])))])
        }
        dens$interval <- as.factor(dens$interval)
        labels <- colnames(t(summary(dens$interval)))
        frequencies <- dens$density
        frequencies_gb <- rep(0,5)
        for (label in 1:5) {
          if (!is.na(match(all_labels[label], labels)))
            frequencies_gb[label] <- frequencies[match(all_labels[label], labels)]
        }
        cpg_list_gb[[gb_cpg]] <- frequencies_gb
      }
      
      # promoter
      cpg_list_pr <- NULL # CpG_P_X.likelihood
      for (pr_cpg in 1:length(IDs_promoter)) {
        x <- dataf_BRCA_logit_PC@bmatrix[IDs_promoter[pr_cpg],current_Ts[current_sample]]
        if (x < -7) x <- -7
        dens <- density(x,bw=0.14,from=-7,to=7,n=1401)
        dens$y <- dens$y/sum(dens$y)
        dens$y <- dens$y * d_promoter$y # P(M|d)*P(d)
        dens <- as.data.frame(cbind(dens$x,dens$y))
        colnames(dens) <- c("M","density")
        dens$interval <- findInterval(dens[,1],breaksPROMOTER)
        dens <- aggregate(cbind(M,density) ~ interval, data=dens, sum)
        for (bin in 2:length(breaksPROMOTER)) { # (P(M|d)*P(d))/P(M)
          dens[bin-1,3] <- dens[bin-1,3] / sum(d_promoter$y[(intersect(which(d_promoter$x >= breaksPROMOTER[bin-1]),which(d_promoter$x < breaksPROMOTER[bin])))])
        }
        dens$interval <- as.factor(dens$interval)
        labels <- colnames(t(summary(dens$interval)))
        frequencies <- dens$density
        frequencies_pr <- rep(0,5)
        for (label in 1:5) {
          if (!is.na(match(all_labels[label], labels)))
            frequencies_pr[label] <- frequencies[match(all_labels[label], labels)]
        }
        cpg_list_pr[[pr_cpg]] <- frequencies_pr
      }
      
      tempS_formated <- matrix()
      tempS_formated[1,1] <- paste('[1,5]((',paste(frequencies_expr,sep="",collapse=","),'))',sep="")
      for (element in 1:length(cpg_list_pr)) {
        tempS_formated <- cbind(tempS_formated,paste('[1,5]((',paste(cpg_list_pr[[element]],sep="",collapse=","),'))',sep=""))
      }
      for (element in 1:length(cpg_list_gb)) {
        tempS_formated <- cbind(tempS_formated,paste('[1,5]((',paste(cpg_list_gb[[element]],sep="",collapse=","),'))',sep=""))
      }
      tempS <- rbind(tempS,tempS_formated)
    }
    tempFac <- tempS[2:71,]
    colnames(tempFac) <- c("NAME:\tEXPR.likelihood",promoter_CpGs[1:length(IDs_promoter)],geneBody_CpGs[1:length(IDs_body)])
    rownames(tempFac) <- current_Ts  
    eval(parse(text = paste('write.table(', paste('tempFac,file ="./',i,'/T_model/',fold,'_T_FacData.tab",row.names=TRUE,col.names=TRUE,quote=FALSE,sep="\t",append=FALSE)', sep = ""))))
    
    # build and query the current x-val fold model
    system(command=paste('mkdir ./',i,'/T_model/',fold,'/',sep=""))
    system(command=paste('../phy/src/dfgTrain --emTrain --maxIter 1000 --dfgSpecPrefix=./',i,'/ --outSpecPrefix=./',i,'/T_model/',fold,'/ ./',i,'/T_model/',fold,'_T_VarData.tab ./',i,'/T_model/',fold,'_T_FacData.tab',sep=""))
    l5o_Ts <- Ts2[seq(from=fold*5-4,to=fold*5,by=1)]
    l5o_which_Ts <- which_Ts[seq(from=fold*5-4,to=fold*5,by=1)]
    
    # expression
    lambda <- counts_BRCA_plusOne[workingList_BRCA4[i],l5o_Ts[1]]
    X <- seq(round(max(lambda-(4*lambda*lambda^(-1/2)),1)),round(lambda+(4*lambda*lambda^(-1/2))))
    list <- NULL
    list[[1]] <- cbind(X/factors_ls[l5o_which_Ts[1]],dpois(X,lambda=lambda)*factors_ls[l5o_which_Ts[1]])
    for (j in 2:length(l5o_Ts)) {
      lambda <- counts_BRCA_plusOne[workingList_BRCA4[i],l5o_Ts[j]]
      X <- seq(round(max(lambda-(4*lambda*lambda^(-1/2)),1)),round(lambda+(4*lambda*lambda^(-1/2))))
      current <- factors_ls[l5o_which_Ts[j]]
      list[[j]] <- cbind(X/current,dpois(X,lambda=lambda)*current)
    }
    
    tempS <- matrix(ncol=ncol)
    for (current_sample in 1:length(l5o_Ts)) {     
      all_labels <- as.character(seq(1,5,1))
      
      # expression
      AN <- list[[current_sample]]
      AN <- as.data.frame(AN)
      colnames(AN) <- c("cpm","density")
      AN$interval <- findInterval(AN[,1],breaksEXPRESSION)
      AN <- aggregate(cbind(cpm,density) ~ interval, data=AN, sum) #sum up the repetitive intervals' densities
      AN$interval <- as.factor(AN$interval)
      labels <- colnames(t(summary(AN$interval)))
      frequencies <- AN$density
      frequencies <- AN$density * nsim/sum(AN$density)
      frequencies_expr <- rep(0,5) # EXPR.likelihood
      for (label in 1:5) {if (!is.na(match(all_labels[label], labels))) {
        frequencies_expr[label] <- frequencies[match(all_labels[label], labels)]
      }}    
      
      # gene body
      cpg_list_gb <- NULL # CpG_GB_X.likelihood
      for (gb_cpg in 1:length(IDs_body)) {
        x <- dataf_BRCA_logit_PC@bmatrix[IDs_body[gb_cpg],l5o_Ts[current_sample]]
        if (x < -7) x <- -7
        dens <- density(x,bw=0.14,from=-7,to=7,n=1401)
        dens$y <- dens$y/sum(dens$y)
        dens$y <- dens$y * d_body$y # P(M|d)*P(d)
        dens <- as.data.frame(cbind(dens$x,dens$y))
        colnames(dens) <- c("M","density")
        dens$interval <- findInterval(dens[,1],breaksBODY)
        dens <- aggregate(cbind(M,density) ~ interval, data=dens, sum)
        for (bin in 2:length(breaksBODY)) { # (P(M|d)*P(d))/P(M)
          dens[bin-1,3] <- dens[bin-1,3] / sum(d_body$y[(intersect(which(d_body$x >= breaksBODY[bin-1]),which(d_body$x < breaksBODY[bin])))])
        }
        dens$interval <- as.factor(dens$interval)
        labels <- colnames(t(summary(dens$interval)))
        frequencies <- dens$density
        frequencies_gb <- rep(0,5)
        for (label in 1:5) {
          if (!is.na(match(all_labels[label], labels)))
            frequencies_gb[label] <- frequencies[match(all_labels[label], labels)]
        }
        cpg_list_gb[[gb_cpg]] <- frequencies_gb
      }
      
      # promoter
      cpg_list_pr <- NULL # CpG_P_X.likelihood
      for (pr_cpg in 1:length(IDs_promoter)) {
        dens <- density(dataf_BRCA_logit_PC@bmatrix[IDs_promoter[pr_cpg],l5o_Ts[current_sample]],bw=0.14,from=-7,to=7,n=1401)
        dens$y <- dens$y/sum(dens$y)
        dens$y <- dens$y * d_promoter$y # P(M|d)*P(d)
        dens <- as.data.frame(cbind(dens$x,dens$y))
        colnames(dens) <- c("M","density")
        dens$interval <- findInterval(dens[,1],breaksPROMOTER)
        dens <- aggregate(cbind(M,density) ~ interval, data=dens, sum)
        for (bin in 2:length(breaksPROMOTER)) { # (P(M|d)*P(d))/P(M)
          dens[bin-1,3] <- dens[bin-1,3] / sum(d_promoter$y[(intersect(which(d_promoter$x >= breaksPROMOTER[bin-1]),which(d_promoter$x < breaksPROMOTER[bin])))])
        }
        dens$interval <- as.factor(dens$interval)
        labels <- colnames(t(summary(dens$interval)))
        frequencies <- dens$density
        frequencies_pr <- rep(0,5)
        for (label in 1:5) {
          if (!is.na(match(all_labels[label], labels)))
            frequencies_pr[label] <- frequencies[match(all_labels[label], labels)]
        }
        cpg_list_pr[[pr_cpg]] <- frequencies_pr
      }
      
      tempS_formated <- matrix()
      tempS_formated[1,1] <- paste('[1,5]((',paste(frequencies_expr,sep="",collapse=","),'))',sep="")
      for (element in 1:length(cpg_list_pr)) {
        tempS_formated <- cbind(tempS_formated,paste('[1,5]((',paste(cpg_list_pr[[element]],sep="",collapse=","),'))',sep=""))
      }
      for (element in 1:length(cpg_list_gb)) {
        tempS_formated <- cbind(tempS_formated,paste('[1,5]((',paste(cpg_list_gb[[element]],sep="",collapse=","),'))',sep=""))
      }
      tempS <- rbind(tempS,tempS_formated)
    }
    tempFac <- tempS[2:6,]
    colnames(tempFac) <- c("NAME:\tEXPR.likelihood",promoter_CpGs[1:length(IDs_promoter)],geneBody_CpGs[1:length(IDs_body)])
    rownames(tempFac) <- l5o_Ts  
    eval(parse(text = paste('write.table(', paste('tempFac,file ="./',i,'/T_model/',fold,'_evalT_FacData.tab",row.names=TRUE,col.names=TRUE,quote=FALSE,sep="\t",append=FALSE)', sep = ""))))
    tempVar <- matrix(rep(".",5*ncol),nrow=5,ncol=ncol)
    colnames(tempVar) <- c("NAME:\tEXPR",promoterVars[1:length(IDs_promoter)],geneBodyVars[1:length(IDs_body)])
    rownames(tempVar) <- l5o_Ts
    eval(parse(text = paste('write.table(', paste('tempVar,file = "./',i,'/T_model/',fold,'_evalT_VarData.tab",row.names=TRUE,col.names=TRUE,quote=FALSE,sep="\t",append=FALSE)', sep = ""))))  
    
    string<-system(intern=TRUE,command=paste('../phy/src/dfgEval --dfgSpecPrefix=./',i,'/T_model/',fold,'/ -l -n - ./',i,'/T_model/',fold,'_evalT_VarData.tab ./',i,'/T_model/',fold,'_evalT_FacData.tab',sep=""))
    Ts_T_likelihoods <- c(Ts_T_likelihoods,string)    
  }
  Ts_T_likelihoods <- as.numeric(substring(Ts_T_likelihoods[-seq(1,90,6)],17))
  
  
  # full T model developed from here, to obtain likelihoods of ANs
  # gene body  
  mean <- mean(dataf_BRCA_logit_PC@bmatrix[IDs_body,Ts2])
  sd <- sd(dataf_BRCA_logit_PC@bmatrix[IDs_body,Ts2])
  breaksBODY <- round(sort(c(-7.01,qnorm(percentiles,mean=mean,sd=sd),7.01)),digits=2)
  
  # promoter  
  mean <- mean(dataf_BRCA_logit_PC@bmatrix[IDs_promoter,Ts2])
  sd <- sd(dataf_BRCA_logit_PC@bmatrix[IDs_promoter,Ts2])
  breaksPROMOTER <- round(sort(c(-7.01,qnorm(percentiles,mean=mean,sd=sd),7.01)),digits=2)
  # expression
  mean <- mean(cpm_BRCA_plusOne[workingList_BRCA4[i],Ts2])
  sd <- sd(cpm_BRCA_plusOne[workingList_BRCA4[i],Ts2])
  breaksEXPRESSION <- sort(c(0,qnorm(percentiles,mean=mean,sd=sd),1000000))
  
  # generate FacData for ANs
  lambda <- counts_BRCA_plusOne[workingList_BRCA4[i],ANs2[1]]
  X <- seq(round(max(lambda-(4*lambda*lambda^(-1/2)),1)),round(lambda+(4*lambda*lambda^(-1/2))))
  list <- NULL
  list[[1]] <- cbind(X/factors_ls[which_ANs[1]],dpois(X,lambda=lambda)*factors_ls[which_ANs[1]])
  for (j in 2:length(ANs2)) {
    lambda <- counts_BRCA_plusOne[workingList_BRCA4[i],ANs2[j]]
    X <- seq(round(max(lambda-(4*lambda*lambda^(-1/2)),1)),round(lambda+(4*lambda*lambda^(-1/2))))
    current <- factors_ls[which_ANs[j]]
    list[[j]] <- cbind(X/current,dpois(X,lambda=lambda)*current)
  }
  tempS <- matrix(ncol=ncol)
  for (i_t in 1:length(ANs2)) {
    
    all_labels <- as.character(seq(1,5,1))
    
    # expression
    T <- list[[i_t]]
    T <- as.data.frame(T)
    colnames(T) <- c("cpm","density")
    T$interval <- findInterval(T[,1],breaksEXPRESSION)
    T <- aggregate(cbind(cpm,density) ~ interval, data=T, sum) #sum up the repetitive intervals' densities
    T$interval <- as.factor(T$interval)
    labels <- colnames(t(summary(T$interval)))
    frequencies <- T$density
    frequencies <- T$density * nsim/sum(T$density)
    frequencies_expr <- rep(0,5) # EXPR.likelihood
    for (label in 1:5) {if (!is.na(match(all_labels[label], labels))) {
      frequencies_expr[label] <- frequencies[match(all_labels[label], labels)]
    }}    
    
    # gene body
    cpg_list_gb <- NULL # CpG_GB_X.likelihood
    for (gb_cpg in 1:length(IDs_body)) {
      x <- dataf_BRCA_logit_PC@bmatrix[IDs_body[gb_cpg],Ts2[i_t]]
      if (x < -7) x <- -7
      dens <- density(x,bw=0.14,from=-7,to=7,n=1401)
      dens$y <- dens$y/sum(dens$y)
      dens$y <- dens$y * d_body$y # P(M|d)*P(d)
      dens <- as.data.frame(cbind(dens$x,dens$y))
      colnames(dens) <- c("M","density")
      dens$interval <- findInterval(dens[,1],breaksBODY)
      dens <- aggregate(cbind(M,density) ~ interval, data=dens, sum)
      for (bin in 2:length(breaksBODY)) { # (P(M|d)*P(d))/P(M)
        dens[bin-1,3] <- dens[bin-1,3] / sum(d_body$y[(intersect(which(d_body$x >= breaksBODY[bin-1]),which(d_body$x < breaksBODY[bin])))])
      }
      dens$interval <- as.factor(dens$interval)
      labels <- colnames(t(summary(dens$interval)))
      frequencies <- dens$density
      frequencies_gb <- rep(0,5)
      for (label in 1:5) {
        if (!is.na(match(all_labels[label], labels)))
          frequencies_gb[label] <- frequencies[match(all_labels[label], labels)]
      }
      cpg_list_gb[[gb_cpg]] <- frequencies_gb
    }
    
    # promoter
    cpg_list_pr <- NULL # CpG_P_X.likelihood
    for (pr_cpg in 1:length(IDs_promoter)) {
      dens <- density(dataf_BRCA_logit_PC@bmatrix[IDs_promoter[pr_cpg],Ts2[i_t]],bw=0.14,from=-7,to=7,n=1401)
      dens$y <- dens$y/sum(dens$y)
      dens$y <- dens$y * d_promoter$y # P(M|d)*P(d)
      dens <- as.data.frame(cbind(dens$x,dens$y))
      colnames(dens) <- c("M","density")
      dens$interval <- findInterval(dens[,1],breaksPROMOTER)
      dens <- aggregate(cbind(M,density) ~ interval, data=dens, sum)
      for (bin in 2:length(breaksPROMOTER)) { # (P(M|d)*P(d))/P(M)
        dens[bin-1,3] <- dens[bin-1,3] / sum(d_promoter$y[(intersect(which(d_promoter$x >= breaksPROMOTER[bin-1]),which(d_promoter$x < breaksPROMOTER[bin])))])
      }
      dens$interval <- as.factor(dens$interval)
      labels <- colnames(t(summary(dens$interval)))
      frequencies <- dens$density
      frequencies_pr <- rep(0,5)
      for (label in 1:5) {
        if (!is.na(match(all_labels[label], labels)))
          frequencies_pr[label] <- frequencies[match(all_labels[label], labels)]
      }
      cpg_list_pr[[pr_cpg]] <- frequencies_pr
    }
    
    tempS_formated <- matrix()
    tempS_formated[1,1] <- paste('[1,5]((',paste(frequencies_expr,sep="",collapse=","),'))',sep="")
    for (element in 1:length(cpg_list_pr)) {
      tempS_formated <- cbind(tempS_formated,paste('[1,5]((',paste(cpg_list_pr[[element]],sep="",collapse=","),'))',sep=""))
    }
    for (element in 1:length(cpg_list_gb)) {
      tempS_formated <- cbind(tempS_formated,paste('[1,5]((',paste(cpg_list_gb[[element]],sep="",collapse=","),'))',sep=""))
    }
    tempS <- rbind(tempS,tempS_formated)
  }
  tempFac <- tempS[2:76,]
  colnames(tempFac) <- c("NAME:\tEXPR.likelihood",promoter_CpGs[1:length(IDs_promoter)],geneBody_CpGs[1:length(IDs_body)])
  rownames(tempFac) <- ANs2
  eval(parse(text = paste('write.table(', paste('tempFac,file = "./',i,'/T_model/AN_FacData.tab",row.names=TRUE,col.names=TRUE,quote=FALSE,sep="\t",append=FALSE)', sep = ""))))
  
  # generate FacData for full set of Ts
  lambda <- counts_BRCA_plusOne[workingList_BRCA4[i],Ts2[1]]
  X <- seq(round(max(lambda-(4*lambda*lambda^(-1/2)),1)),round(lambda+(4*lambda*lambda^(-1/2))))
  list <- NULL
  list[[1]] <- cbind(X/factors_ls[which_Ts[1]],dpois(X,lambda=lambda)*factors_ls[which_Ts[1]])
  for (j in 2:length(Ts2)) {
    lambda <- counts_BRCA_plusOne[workingList_BRCA4[i],Ts2[j]]
    X <- seq(round(max(lambda-(4*lambda*lambda^(-1/2)),1)),round(lambda+(4*lambda*lambda^(-1/2))))
    current <- factors_ls[which_Ts[j]]
    list[[j]] <- cbind(X/current,dpois(X,lambda=lambda)*current)
  }
  tempS <- matrix(ncol=ncol)
  for (current_sample in 1:length(Ts2)) {
    
    all_labels <- as.character(seq(1,5,1))
    
    # expression
    AN <- list[[current_sample]]
    AN <- as.data.frame(AN)
    colnames(AN) <- c("cpm","density")
    AN$interval <- findInterval(AN[,1],breaksEXPRESSION)
    AN <- aggregate(cbind(cpm,density) ~ interval, data=AN, sum) #sum up the repetitive intervals' densities
    AN$interval <- as.factor(AN$interval)
    labels <- colnames(t(summary(AN$interval)))
    frequencies <- AN$density
    frequencies <- AN$density * nsim/sum(AN$density)
    frequencies_expr <- rep(0,5) # EXPR.likelihood
    for (label in 1:5) {if (!is.na(match(all_labels[label], labels))) {
      frequencies_expr[label] <- frequencies[match(all_labels[label], labels)]
    }}    
    
    # gene body
    cpg_list_gb <- NULL # CpG_GB_X.likelihood
    for (gb_cpg in 1:length(IDs_body)) {
      x <- dataf_BRCA_logit_PC@bmatrix[IDs_body[gb_cpg],Ts2[current_sample]]
      if (x < -7) x <- -7
      dens <- density(x,bw=0.14,from=-7,to=7,n=1401)
      dens$y <- dens$y/sum(dens$y)
      dens$y <- dens$y * d_body$y # P(M|d)*P(d)
      dens <- as.data.frame(cbind(dens$x,dens$y))
      colnames(dens) <- c("M","density")
      dens$interval <- findInterval(dens[,1],breaksBODY)
      dens <- aggregate(cbind(M,density) ~ interval, data=dens, sum)
      for (bin in 2:length(breaksBODY)) { # (P(M|d)*P(d))/P(M)
        dens[bin-1,3] <- dens[bin-1,3] / sum(d_body$y[(intersect(which(d_body$x >= breaksBODY[bin-1]),which(d_body$x < breaksBODY[bin])))])
      }
      dens$interval <- as.factor(dens$interval)
      labels <- colnames(t(summary(dens$interval)))
      frequencies <- dens$density
      frequencies_gb <- rep(0,5)
      for (label in 1:5) {
        if (!is.na(match(all_labels[label], labels)))
          frequencies_gb[label] <- frequencies[match(all_labels[label], labels)]
      }
      cpg_list_gb[[gb_cpg]] <- frequencies_gb
    }
    
    # promoter
    cpg_list_pr <- NULL # CpG_P_X.likelihood
    for (pr_cpg in 1:length(IDs_promoter)) {
      dens <- density(dataf_BRCA_logit_PC@bmatrix[IDs_promoter[pr_cpg],ANs2[current_sample]],bw=0.14,from=-7,to=7,n=1401)
      dens$y <- dens$y/sum(dens$y)
      dens$y <- dens$y * d_promoter$y # P(M|d)*P(d)
      dens <- as.data.frame(cbind(dens$x,dens$y))
      colnames(dens) <- c("M","density")
      dens$interval <- findInterval(dens[,1],breaksPROMOTER)
      dens <- aggregate(cbind(M,density) ~ interval, data=dens, sum)
      for (bin in 2:length(breaksPROMOTER)) { # (P(M|d)*P(d))/P(M)
        dens[bin-1,3] <- dens[bin-1,3] / sum(d_promoter$y[(intersect(which(d_promoter$x >= breaksPROMOTER[bin-1]),which(d_promoter$x < breaksPROMOTER[bin])))])
      }
      dens$interval <- as.factor(dens$interval)
      labels <- colnames(t(summary(dens$interval)))
      frequencies <- dens$density
      frequencies_pr <- rep(0,5)
      for (label in 1:5) {
        if (!is.na(match(all_labels[label], labels)))
          frequencies_pr[label] <- frequencies[match(all_labels[label], labels)]
      }
      cpg_list_pr[[pr_cpg]] <- frequencies_pr
    }
    
    tempS_formated <- matrix()
    tempS_formated[1,1] <- paste('[1,5]((',paste(frequencies_expr,sep="",collapse=","),'))',sep="")
    for (element in 1:length(cpg_list_pr)) {
      tempS_formated <- cbind(tempS_formated,paste('[1,5]((',paste(cpg_list_pr[[element]],sep="",collapse=","),'))',sep=""))
    }
    for (element in 1:length(cpg_list_gb)) {
      tempS_formated <- cbind(tempS_formated,paste('[1,5]((',paste(cpg_list_gb[[element]],sep="",collapse=","),'))',sep=""))
    }
    tempS <- rbind(tempS,tempS_formated)
  }
  tempFac <- tempS[2:76,]
  colnames(tempFac) <- c("NAME:\tEXPR.likelihood",promoter_CpGs[1:length(IDs_promoter)],geneBody_CpGs[1:length(IDs_body)])
  rownames(tempFac) <- Ts2  
  eval(parse(text = paste('write.table(', paste('tempFac,file ="./',i,'/T_model/T_FacData.tab",row.names=TRUE,col.names=TRUE,quote=FALSE,sep="\t",append=FALSE)', sep = ""))))
  
  # build and query the full model with AN samples
  system(command=paste('mkdir ./',i,'/T_model/','all/',sep=""))
  system(command=paste('../phy/src/dfgTrain --emTrain --maxIter 1000 --dfgSpecPrefix=./',i,'/ --outSpecPrefix=./',i,'/T_model/all/ ./',i,'/T_model/T_VarData.tab ./',i,'/T_model/T_FacData.tab',sep=""))
  string<-system(intern=TRUE,command=paste('../phy/src/dfgEval --dfgSpecPrefix=./',i,'/T_model/',fold,'/ -l -n - ./',i,'/T_model/AN_VarData.tab ./',i,'/T_model/AN_FacData.tab',sep=""))
  ANs_T_likelihoods <- as.numeric(substring(string[-1],17))
  ##########################################################################
  require(Brobdingnag)
  ANs_T_likelihoods <- as.brob(exp(1)^-ANs_T_likelihoods)
  ANs_AN_likelihoods <- as.brob(exp(1)^-ANs_AN_likelihoods)
  Ts_T_likelihoods <- as.brob(exp(1)^-Ts_T_likelihoods)
  Ts_AN_likelihoods <- as.brob(exp(1)^-Ts_AN_likelihoods)
  list_out[[i]] <- list(ANs_AN_likelihoods,ANs_T_likelihoods,Ts_T_likelihoods,Ts_AN_likelihoods)
  # calculate P(class=AN|D)
  #P_t <- 0.1
  #P_an <- 0.9
  #P_cAN_dAN <- as.numeric((ANs_AN_likelihoods*P_an) /(ANs_AN_likelihoods*P_an+ANs_T_likelihoods*P_t))
  #P_cAN_dT <- as.numeric((Ts_AN_likelihoods*P_an) /(Ts_AN_likelihoods*P_an+Ts_T_likelihoods*P_t))
  # calculate P(class=T|D)
  #P_cT_dAN <- as.numeric((ANs_T_likelihoods*P_t) /(ANs_AN_likelihoods*P_an+ANs_T_likelihoods*P_t))
  #P_cT_dT <- as.numeric((Ts_T_likelihoods*P_t) /(Ts_AN_likelihoods*P_an+Ts_T_likelihoods*P_t))
  
}