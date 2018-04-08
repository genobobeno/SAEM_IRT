####
# Analysis for Cai (2010) polytomous data. Get libraries for analysis
# This is working unbelievably well. Need to validate. Big time. Ply version 4.
#

rm(list=ls(all=TRUE))
library(rlecuyer)		# for rand num generation
library(snow)		# for parallel processing
library(GPArotation)	# for rotations
library(mvnfast)		# for function mvrnorm
library(psych)		# for ML factor analysis

drawAEigen <- function(covZ,Q,m,ncat)
{
  # Subtract diag(m)
  covar 	<- covZ - diag(m)
  out 		<- eigen( covar, symmetric=TRUE)
  Avec 		<- out$values #; print(round(Avec[1:Q],2))
  Aload 	<- out$vectors
  Atemp 	<- Aload[,1:Q]%*%sqrt(diag(Avec[1:Q]))
  return(list(Atemp=Atemp,Avec=Avec))}

# Number of items (m), schools (n), dimensions (Q)
m <- 100; n <- 100000 ; Q <- 10; ncat <- 4; IQ <- diag(Q)

# Read data
load("C:\\Users\\gregory.camilli\\Desktop\\Eugene\\Poly_J100_N1e+05_Q10_K4.rda")
y.a <- gen.rp
y.a = data.matrix(y.a)
# Set up n x m x (ncat-1) data array of binary values for cumulative option indicators
y.b = array(NA,c(m,ncat-1,n))
for (i in 1:n) 
  for (j in 1:m)
    # for missing responses, item propensities = 9 
    if (y.a[i,j]==9)	{y.b[j,,i] <- rep(9,(ncat-1))} else
    {
      for (k in 0:(ncat-2))
      {		
        if (y.a[i,j]<=k) 	{y.b[j,(k+1),i] <- 0} 
        else 			{y.b[j,(k+1),i] <- 1} 
      }}
Y = y.b

niter 	<-   150	
nEM		<-   50	#number of burnin EM cycles < niter
nproc 	<-   10
cl <- makeCluster(nproc,type="SOCK")
clusterSetupRNG(cl, seed = round(runif(6)*1001))
clusterEvalQ(cl,library(mvnfast))
clusterExport(cl,c("Y","n","m","ncat","Q"))
nn = rep(1:n)
clusterExport(cl,c("y.a","nn"))





# wrapper for Z. Returns single Z (instead of mean Z) for estimating theta.
wrapZ <- function(j) {
  
  Zj = matrix(NA,n,1)
  eta  	<- theta%*%t(t(A[j,]))
  bd	<- matrix(b[j]+d[j,],n,ncat-1,byrow=TRUE)
  hold <- sweep(bd,1,eta,"-")
  pp <- cbind(-Inf,hold,Inf)
  U  <- matrix(runif(n),n,1)
  
  # Trim extreme values
  U[U>0.99999] <- 0.99999; U[U<0.00001] <- 0.00001
  
  # Item propensities for all examinees from truncated normal
  yL <- matrix(y.a[,j]+1,n,1) 
  indL <- cbind(nn,yL); 		indU <- cbind(nn,yL + 1)
  pL <- pnorm( pp[indL] ) ; 	pU <- pnorm( pp[indU] )
  
  Zj <- eta + qnorm( pL + U*(pU - pL)  ) 
  c(Zj)
}
clusterEvalQ(cl,"wrapZ")

# wrapper for snow. m x (ncat-1) matrix is returned for each person, 
# for estimating item thresholds
wrapX <- function(j) {
  
  # ker is m x (ncat-1) matrix of item kernels
  ker <- apply(theta%*%matrix(A[j,],Q,1)-b[j],1,
               function(x) x - c(d[j,]))
  pn 	<- pnorm(-ker)
  r  	= Y[j,,]
  Xi	= matrix(NA,ncat-1,n)
  
  # Count missing and nonmissing values (9)
  miss	 	<- sum(r==9)
  nmis	 	<- n*(ncat-1)- miss
  mcol 		<- nmis/(ncat-1)
  
  # extract missing & nonmissing reponses
  my <- which(r == 9,arr.ind=TRUE)
  mn <- which(r != 9,arr.ind=TRUE)
  
  # Generate missing option propensities
  Xi[my] <- rnorm(miss)
  
  # Generate nonmissing option propensities
  r1 		<- matrix(r[mn]      ,ncat-1,mcol)
  U  		<- matrix(runif(nmis),ncat-1,mcol)
  P		<- matrix(pn[mn]     ,ncat-1,mcol)
  Xi[mn] 	<- qnorm( r1*U + P*(r1 + U - 2*r1*U) )
  Xi		<- ker + Xi
  return(Xi)}
clusterEvalQ(cl,"wrapX")



# wrapT is wrapper function for estimating theta
wrapT <- function(i)             
{
  IQ 		<- diag(Q)
  ATZ   	<- t(A)%*%(Z[i,])
  That    	<- BTB_INV%*%ATZ
  ttempi 	<- rmvn(1,That,BTB_INV)
  return(ttempi)}
clusterEvalQ(cl,"wrapT")

# Starting values 
theta 	= matrix(rnorm(n*Q,0,1),n,Q)
b		= matrix(.5*rnorm(m),m,1);	
d		= matrix(c(-0.50, 0, 0.50),m,ncat-1, byrow=TRUE)
A		= matrix(.2,m,Q); Aold1 <- matrix(runif(m*Q),m,Q)

AA=gen.xi[,1:10]; bb=gen.xi[,11]; dd=gen.tau-rowMeans(gen.tau); ttheta=gen.theta
clusterExport(cl,c("A","b","d","theta"))		
meanD <- 0; meanB <- 0; covZ <- 0; plot(1,1); oldvec <- rep(0,Q)
Aold2 <- matrix(0,m,Q)

# Now iterate
Tstart <- Sys.time()
for (i in 1:niter) {
  if (i < (nEM+1)) alpha <- 1 else alpha <- 1/(i-nEM)
  print(c(i,alpha))
  
  X1 		<- parSapply(cl,1:m,wrapX,simplify=FALSE)
  X2		<- simplify2array(X1, higher=TRUE)	
  X3		<- t(apply(X2,c(1,3),mean))
  meanBMC	<- t(t(rowMeans(X3)))
  meanDMC 	<- t(apply(X3, 1, scale, scale=FALSE)) 
  
  # meanB, meanD are suff stats for item difficulty & thresholds
  meanB	<-  meanB + alpha*(meanBMC-meanB)
  meanD	<-  meanD + alpha*(meanDMC-meanD)
  
  b 	<- -meanB;	d		<- -meanD
  
  clusterExport(cl,c("b","d"))
  # Get Z to estimate theta
  X1 		<- parSapply(cl,1:m,wrapZ,simplify=FALSE)
  Z		<- simplify2array(X1, higher=TRUE)	
  clusterExport(cl,"Z") 
  
  #Compute covariances
  covZZMC	<- cov(Z)
  covZ		<- covZ  + alpha*(covZZMC-covZ)
  out <- drawAEigen(covZ,Q,m,ncat)
  A <- out$Atemp; Anew <- A
  newvec <- out$Avec
  dif <- abs(newvec-oldvec); dif = max(dif); oldvec <- newvec
  cat (round(dif,5),"Eigen Change","\n")
  
  clusterExport(cl,c("A","b","d"))
  
  
  #Draw theta
  ATA 		<- t(A)%*%A
  BTB_INV	<- solve(IQ + ATA)
  clusterExport(cl,"BTB_INV")
  
  thetanew	<- t(parSapply(cl,1:n,wrapT))
  theta		<- thetanew
  clusterExport(cl,"theta")
  
  ##################################
  # Following code is only for diagnostics – doesn’t affect estimation
  #
  
  for (ii in 1:Q) if (i > 4){
    if (cor(Aold1[,ii],Anew[,ii]) < -.25) Anew[,ii] <- -Anew[,ii]}
  Aold1 <- Anew
  
  RotT = targetT(Anew,Tmat=diag(Q),Target=AA,normalize=TRUE,eps=1e-4, maxit=1000)
  A2   <- RotT$loadings
  
  RotV <- Varimax(A2,Tmat=diag(Q),normalize=TRUE,eps=1e-4,maxit=100)
  A2   = RotV$loadings
  
  for (ii in 1:Q) if (abs(min(A2[,ii]))> max(A2[,ii]) ) A2[,ii]=-A2[,ii]
  
  MSE = sqrt((sum((A2-AA)^2))/(Q*m))
  cat (round(MSE,5),"MSE_A  ","\n")
  MSE = (sum(abs(bb-b)))/(Q*m); cat (round(MSE,5),"MSE_b  ","\n")
  MSE = (sum(abs(dd-d)))/(Q*m); cat (round(MSE,5),"MSE_d  ","\n")
  dif <- A2 - Aold2
  Aold2 <- A2
  dif <- max(abs(dif)); cat (round(dif,5),"MAX_dif","\n")
  
  
  dev.off()
  dev.new(width=12, height=10)
  # 10 figures arranged in 2 rows and 5 columns
  par(mfrow=c(2,6))
  for (k in 1:Q)
  {
    plot(AA[,k], A2[,k],ylab=paste("Est A",k),
         xlab = paste("True A",k)); abline(a=0,b=1)
  }
  plot(bb,b, ylab=paste("Est b"), xlab=paste("True b")) 
  abline(a=0,b=1)
  
  plot(dd[,2],d[,2], ylab=paste("Est d",2),	xlab = paste("True d",2))
  abline(a=0,b=1)
  cat ("\n")
  #
  # End diagnostics
  #########################################################
  
  print(c(max(theta),min(theta),max(A2),min(A2),max(Z),min(Z)))	
}

Time <- Sys.time()-Tstart;	Time
stopCluster(cl)


dev.off()
dev.new(width=12, height=10)
flex = read.csv("C:/Users/gregory.camilli/Desktop/eugene/new100k.csv",header=FALSE)
flex=data.matrix(flex)

RotV <- Varimax(flex/1.749,Tmat=diag(Q),normalize=TRUE,eps=1e-4,maxit=1000)
A2a   = RotV$loadings

RotV <- Varimax(-flex/1.749,Tmat=diag(Q),normalize=TRUE,eps=1e-4,maxit=1000)
A2b   = RotV$loadings

RotT = targetT(A2a,Tmat=diag(Q),Target=AA,normalize=TRUE,eps=1e-4, maxit=1000)
A2a   <- RotT$loadings

RotT = targetT(A2b,Tmat=diag(Q),Target=AA,normalize=TRUE,eps=1e-4, maxit=1000)
A2b   <- RotT$loadings

sa = sum(abs(A2a) < .05); sb = sum(abs(A2b) < .05)

if (sa > sb) {A2 = A2a} else {Anew=-Anew; A2 = A2b}


# 10 figures arranged in 2 rows and 5 columns
par(mfrow=c(2,5))
for (k in 1:10)
{
  plot(AA[,k], A2[,k],ylab=paste("Est A",k),
       xlab = paste("True A",k)); abline(a=0,b=1)
  abline(a=0,b=1)
}

RotV <- Varimax(Anew/1.749,Tmat=diag(Q),normalize=TRUE,eps=1e-4,maxit=1000)
A2a   = RotV$loadings

RotV <- Varimax(-Anew/1.749,Tmat=diag(Q),normalize=TRUE,eps=1e-4,maxit=1000)
A2b   = RotV$loadings

RotT = targetT(A2a,Tmat=diag(Q),Target=AA,normalize=TRUE,eps=1e-4, maxit=1000)
A2a   <- RotT$loadings

RotT = targetT(A2b,Tmat=diag(Q),Target=AA,normalize=TRUE,eps=1e-4, maxit=1000)
A2b   <- RotT$loadings

sa = sum(abs(A2a) < .05); sb = sum(abs(A2b) < .05)

if (sa > sb) {A2 = A2a} else {Anew=-Anew; A2 = A2b}

dev.off()
dev.new(width=12, height=10)
# 10 figures arranged in 2 rows and 5 columns
par(mfrow=c(2,5))
for (k in 1:Q)
{
  plot(AA[,k], A2[,k],ylab=paste("Est A",k),
       xlab = paste("True A",k)); abline(a=0,b=1)
}






