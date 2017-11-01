####
# Analysis for Cai (2010) polytomous data
# Get libraries for analysis
#####

# install.packages("rlecuyer",dependencies = T)
# install.packages("snow",dependencies = T)
# install.packages("GPArotation",dependencies = T)
# install.packages("mvnfast",dependencies = T)
# install.packages("psych",dependencies = T)
# install.packages("foreach",dependencies = T)


library(rlecuyer)		# for rand num generation
library(snow)			# for parallel processing
library(GPArotation)	# for rotations
library(mvnfast)		# for function mvrnorm
library(psych)			# for ML factor analysis
library(foreach)

setwd("//unraidtower/storage/Documents/SAEM/SAEM_IRT")
J=200;N=10000;Q=10;K=4
#J=20;N=1000;Q=2;K=4
niter<-500;nEM<-200	#number of burnin EM cycles < niter
RMwindow<-100;estgain=1;eps=0.000001
# nproc 	<- 2

csvfile = paste0("Poly_Q",Q,"_K",K,"_N",N,"_ResponseData.csv")
datafile = paste0("Poly_J",J,"_N",N,"_Q",Q,"_K",K,".rds")
GEN.DATA = readRDS(datafile)
ncat <- K
IQ <- diag(Q) 
TrueTheta = data.matrix(read.csv(csvfile,header=TRUE)[,1:Q])
# Read data
# y.a = scan("form3_social.dat", what = "numeric",sep = "\n")
y.a = data.matrix(read.csv(csvfile,header=TRUE)[,-(1:Q)]-1)
# y.a = data.matrix(do.call(rbind,lapply(strsplit(y.a,"\\s"), function(x) as.numeric(x[-1]))))
# Set up n x m x (ncat-1) data array of binary values for
m = ncol(y.a); n = nrow(y.a)

# DrawA generates A factor coefficients w/ eigenanalysis
drawA <- function(covZ,Q,m,a)
{
  out 		<- eigen( covZ, symmetric=TRUE)
  if (!is.na(a)) {
    if (Q>1) {
      Aload<-as.matrix(apply(rbind(out$vectors[,1:Q],a),2,
                             function(x) (if (mean(x[1:(length(x)/2)]*x[length(x)/2+1:(length(x)/2)])>0 & sum(x[1:(length(x)/2)])>0) {return(x[1:(length(x)/2)])} 
                                          else {return(-1*x[1:(length(x)/2)])})),nrow(covZ),Q)
      Atemp = out$vectors[,1:Q]%*%sqrt(diag(out$values[1:Q]))
    } else {
      ifelse(mean(out$vectors[,1]*aa)>0,Aload<-out$vectors[,1],Aload<-(-1)*out$vectors[,1])
      Atemp = as.matrix(Aload*sqrt(out$values[1]))
    }
  } else {
    Atemp 	<- out$vectors[,1:Q]%*%sqrt(diag(out$values[1:Q]))
  }
  return(Atemp)}

# cumulative option indicators
y.b = array(NA,c(m,ncat-1,n))
for (i in 1:n) 
  for (j in 1:m)
    # for missing responses, item propensities = 9 
    if (y.a[i,j]==9)	{y.b[j,,i] <- rep(9,(ncat-1))} else
    {
      for (k in 0:(ncat-2))
      {		
        y.b[j,(k+1),i] <- ifelse(y.a[i,j]<=k,0,1) 
      }}
Y = y.b

# nproc 	<- 24

procs = nproc = parallel::detectCores()
#for (procs in 1:nproc) {
  parallelCluster <- parallel::makeCluster(procs,type="SOCK")
  print(parallelCluster)
  clusterSetupRNG(parallelCluster, seed = round(runif(6)*1001))
  clusterEvalQ(parallelCluster,library(mvnfast))
  
  # wrapper for snowfall called by drawX. For each person,
  # m x (ncat-1) matrix is returned for each person
  n1cat = (ncat-1)
  MY = lapply(1:(dim(Y)[1]),function(x) (which(Y[x,,]==9,arr.ind = TRUE)))
  MN = lapply(1:(dim(Y)[1]),function(x) (which(Y[x,,]!=9,arr.ind = TRUE)))
  R  = lapply(1:(dim(Y)[1]),function(x) Y[x,,])
  missList = lapply(MY,function(x) (list(miss = nrow(x),nmis = n*n1cat-nrow(x),mcol = (n*n1cat-nrow(x))/n1cat)))
  
  # wrapper for snowfall called by drawX. For each person,
  # m x (ncat-1) matrix is returned for each person
  wrapX <- function(j,A,b,d,theta) {
    # ker is m x (ncat-1) matrix of item kernels
    ker <- apply(theta%*%matrix(A[j,],Q,1)-b[j],1,
                 function(x) x - c(d[j,]))
    Xi	= matrix(NA,n1cat,n)
    # Generate missing option propensities
    Xi[MY[[j]]] <- rnorm(missList[[j]]$miss)
    # Generate nonmissing option propensities
    r1 		<- matrix(R[[j]][MN[[j]]]           ,n1cat,missList[[j]]$mcol)
    U  		<- matrix(runif(missList[[j]]$nmis) ,n1cat,missList[[j]]$mcol)
    P		  <- matrix(pnorm(-ker)[MN[[j]]]      ,n1cat,missList[[j]]$mcol)
    Xi[MN[[j]]]  <-  qnorm( r1*U + P*(r1 + U - 2*r1*U) )
    return(ker + Xi)}
  
  # wrapT is wrapper function for drawT
  wrapT <- function(i,A,Z,BTB_INV,b)             
  { return(rmvn(1,BTB_INV%*%(t(A)%*%(Z[i,] + b)),BTB_INV))}
  clusterExport(parallelCluster,c("Y","n","m","ncat","Q","MY","MN","R","n1cat","missList"))
  clusterEvalQ(parallelCluster,c("wrapX","wrapT"))
  
  # Starting values 
  theta = matrix(rnorm(n*Q,0,1),n,Q)
  A		= matrix(runif(m*Q)-.5,m,Q)
  b		= matrix(rnorm(m),m,1);	
  d		= matrix(seq(-1.0,1.0,length.out = n1cat),m,n1cat, byrow=TRUE)
  
  X2		<- simplify2array(parSapply(parallelCluster,1:m,wrapX,A=A,b=b,d=d,theta=theta,simplify=FALSE), higher=TRUE)	
  Z		<- apply(X2,c(2,3),mean) 
  covZ  <- cov(Z)
  A	<- drawA(covZ,Q,m,a=NA)
  ATA 		<- t(A)%*%A #*4
  BTB_INV	<- solve(IQ + ATA)
  theta	<- t(parSapply(parallelCluster,1:n,wrapT,A=A,Z = Z,BTB_INV=BTB_INV,b=b))
  meanX	<- 0; meanD <- 0; meanB <- 0
  # alpha <- ifelse(1:niter < (nEM+1),rep(1,niter),1/(1:niter-nEM)^(2/3))
  
  alpha<-c(rep(1,nEM-RMwindow),
           runif(rep(1,ceiling(RMwindow/2)),min = 1.0/(1:(RMwindow/2))^estgain, max = 1.0),1.0,
           runif(rep(1,ceiling(RMwindow/2)),min = 1.0/((RMwindow/2+1):RMwindow)^estgain, max = 1.0/(1:(RMwindow/2))^estgain),
           1/(ceiling(RMwindow/2):(niter-nEM+ceiling(RMwindow/2)))^estgain)  
  
  # Now iterate
  Aiter = array(GEN.DATA$XI[,1:Q],dim=c(m,Q,niter+1))
  Aiter[,,1]<-A
  #Rprof()
  Tstart <- Sys.time()
  for (i in 1:niter) {
    #if (i%%100==0) cat(nproc,i,"|")
    cat(c(".",":","\n")[i%%c(10,100,1000)==0])
    X2		<- simplify2array(parSapply(parallelCluster,1:m,wrapX,simplify=FALSE,A=A,b=b,d=d,theta=theta), higher=TRUE)	
    X3		<- t(apply(X2,c(1,3),mean))
    meanB <- meanB + alpha[i]*(t(t(rowMeans(X3)))-meanB)
    meanD <- meanD + alpha[i]*(t(apply(X3, 1, scale, scale=FALSE)) - meanD)
    b         <- -meanB;    d        <- -meanD   
    #Z		<- apply(X2,c(2,3),mean)
    Z		<- colMeans(X2)
    covZ	<- covZ + alpha[i]*(cov(Z)-covZ)
    if (i<20) {prevA <- NA} else {prevA <- A}
    A <- Anew		<- drawA(covZ-diag(m)/n1cat,Q,m,a=prevA)
    ATA 		<- t(A)%*%A
    BTB_INV	<- solve(IQ + ATA)
    theta		<- thetanew	<- t(parSapply(parallelCluster,1:n,wrapT,A=A,Z=Z,BTB_INV=BTB_INV,b=b))
    # Aiter[,,i+1]<-A
    # if (i%%100==0) {
    #   par(mfrow=c(1,Q))
    #   for (qq in 1:Q) {
    #     plot(c(1,niter),c(-2,2))
    #     abline(h=GEN.DATA$XI[,qq],col=1:m)
    #     for (jj in 1:m) {lines(1:i,Aiter[jj,qq,1:i],col=jj)}
    #   }
    # }
  }
  Time <- Sys.time()-Tstart;	print(paste(nproc,"processors:",Time))
  #summaryRprof()
  if(!is.null(parallelCluster)) {
    stopCluster(parallelCluster)
    parallelCluster <- c()
  }
#}

  
  # set up bifactor target matrix
  cor(as.vector(d),as.vector(sweep(GEN.DATA$TAU,1,GEN.DATA$XI[,Q+1],"-")))
  lm(as.vector(sweep(GEN.DATA$TAU,1,GEN.DATA$XI[,Q+1],"-"))~as.vector(d))
  cor(as.vector(d),as.vector(GEN.DATA$TAU))
  cor(as.vector(b),as.vector(GEN.DATA$XI[,Q+1]))
  cor(as.vector(A),as.vector(GEN.DATA$XI[,1:Q]))
  cor(as.vector(GEN.DATA$THETA),as.vector(theta))
  
  ############## GREG's reconstruction
  # B.TMS <- A
  # B.TMS <- data.matrix(B.TMS)
  # #B  <- Varimax(B.TMS,Tmat=diag(ncol(B.TMS)),normalize=TRUE,eps=1e-5, maxit=10000)
  # B <- oblimin(B.TMS,Tmat=diag(ncol(B.TMS)),normalize=TRUE,eps=1e-5, maxit=10000,gam=0)
  
  mla <- fa(covZ,nfactors=Q,covar=TRUE,rotate="oblimin",fm="ml")
  covar     <- covZ - diag(m)
  out         <- eigen( covar, symmetric=TRUE)
  Avec     <- out$values
  
  # sweep average (b) out of gen.tau
  LL=sweep(GEN.DATA$TAU,1,GEN.DATA$XI[,Q+1])
  
  # examine plot of estimated thresholds (d) again swept gen.tau
  for ( i in 1:(K-1)) print(summary(lm(d[,i]~LL[,i])))
  
  # get eigenvectors for A
  out         <- eigen(A[,1:Q]%*%t(A[,1:Q]))
  Avec     <- out$values
  Aload     <- out$vectors
  TA     <- Aload[,1:Q]%*%sqrt(diag(Avec[1:Q]))
  
  # set up bifactor target matrix
  Target = matrix(NA,m,Q)
  Target[GEN.DATA$XI[,-(Q+1)]==0] = 0
  Target = matrix(Target,m,Q)
  RotA = targetQ(TA, Tmat=diag(Q), Target=Target, normalize=FALSE, 
                 eps=1e-4, maxit=10000)
  
  par(mfrow=c(2,Q/2))
  for (qq in 1:Q) {
    RotA$loadings[,qq]<-ifelse(sum(RotA$loadings[,qq])<0,-1,1)*RotA$loadings[,qq]
    plot(cbind(GEN.DATA$XI[,qq],RotA$loadings[,qq]),xlab="Actual",ylab="Fit",main=paste0("Q",qq))
    abline(0,1)
    cor(-RotA$loadings[,qq],GEN.DATA$XI[,qq])
    lm(-RotA$loadings[,qq]~GEN.DATA$XI[,qq])
  }