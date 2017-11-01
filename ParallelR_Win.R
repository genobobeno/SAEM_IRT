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

# DrawA generates A factor coefficients w/ eigenanalysis
drawA <- function(covZ,Q,m)
{
  # Subtract diag(m)
  out 		<- eigen(covZ - diag(m)/4, symmetric=TRUE)
  return(out$vectors[,1:Q]%*%sqrt(diag(out$values[1:Q])))}

# Read data

y.a = scan("form3_social.dat", what = "numeric",sep = "\n")
y.a = data.matrix(do.call(rbind,lapply(strsplit(y.a,"\\s"), function(x) as.numeric(x[-1]))))
# y.a = read.csv("form3_social.csv", header=FALSE)
# y.a = data.matrix(y.a)
# Number of items (m), schools (n), dimensions (Q)
m <- ncol(y.a); n <- nrow(y.a) ; Q <- 3; ncat <- 5; IQ <- diag(Q); n1cat = ncat-1
# Set up n x m x (ncat-1) data array of binary values for
# cumulative option indicators
y.b = array(NA,c(m,n1cat,n))
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

niter 	<- 200	
nEM		<- 100	#number of burnin EM cycles < niter
# nproc 	<- 24

nproc = parallel::detectCores()
for (procs in 1:nproc) {
  parallelCluster <- parallel::makeCluster(procs,type="SOCK")
  print(parallelCluster)
  clusterSetupRNG(parallelCluster, seed = round(runif(6)*1001))
  clusterEvalQ(parallelCluster,library(mvnfast))
  
  # wrapper for snowfall called by drawX. For each person,
  # m x (ncat-1) matrix is returned for each person
  MY = lapply(1:(dim(Y)[1]),function(x) (which(Y[x,,]==9,arr.ind = TRUE)))
  MN = lapply(1:(dim(Y)[1]),function(x) (which(Y[x,,]!=9,arr.ind = TRUE)))
  R  = lapply(1:(dim(Y)[1]),function(x) Y[x,,])
  missList = lapply(MY,function(x) (list(miss = nrow(x),nmis = n*n1cat-nrow(x),mcol = (n*n1cat-nrow(x))/4)))
  
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
  { return(rmvn(1,BTB_INV%*%(t(A)%*%(Z[i,] + b)*4),BTB_INV))}
  clusterExport(parallelCluster,c("Y","n","m","ncat","Q","MY","MN","R","n1cat","missList"))
  clusterEvalQ(parallelCluster,c("wrapX","wrapT"))
  
  # Starting values 
  theta = matrix(rnorm(n*Q,0,1),n,Q)
  A		= matrix(.25*(runif(m*Q)-.5),m,Q)
  b		= matrix(rnorm(m),m,1);	
  d		= matrix(c(-1.0, -0.5, 0.5, 1.0),m,n1cat, byrow=TRUE)
  
  X2		<- simplify2array(parSapply(parallelCluster,1:m,wrapX,A=A,b=b,d=d,theta=theta,simplify=FALSE), higher=TRUE)	
  Z		<- apply(X2,c(2,3),mean) 
  covZ  <- cov(Z)
  A	<- drawA(covZ,Q,m)
  ATA 		<- t(A)%*%A*4
  BTB_INV	<- solve(IQ + ATA)
  theta	<- t(parSapply(parallelCluster,1:n,wrapT,A=A,Z = Z,BTB_INV=BTB_INV,b=b))
  meanX	<- 0; meanD <- 0; meanB <- 0
  alpha <- ifelse(1:niter < (nEM+1),rep(1,niter),1/(1:niter-nEM)^(2/3))
  
  # Now iterate
  Tstart <- Sys.time()
  for (i in 1:niter) {
    X2		<- simplify2array(parSapply(parallelCluster,1:m,wrapX,simplify=FALSE,A=A,b=b,d=d,theta=theta), higher=TRUE)	
    X3		<- t(apply(X2,c(1,3),mean))
    b <- meanB <- meanB + alpha[i]*(t(t(rowMeans(X3)))-meanB)
    d <- meanD <- meanD + alpha[i]*(t(apply(X3, 1, scale, scale=FALSE)) - meanD)
    Z		<- colMeans(X2)
    covZ	<- covZ + alpha[i]*(cov(Z)-covZ)
    A <- Anew		<- drawA(covZ,Q,m)
    ATA 		<- t(A)%*%A*4
    BTB_INV	<- solve(IQ + ATA)
    theta		<- thetanew	<- t(parSapply(parallelCluster,1:n,wrapT,A=A,Z=Z,BTB_INV=BTB_INV,b=b))
  }
  Time <- Sys.time()-Tstart;	print(paste(nproc,"processors:",Time))
  if(!is.null(parallelCluster)) {
    stopCluster(parallelCluster)
    parallelCluster <- c()
  }
}

