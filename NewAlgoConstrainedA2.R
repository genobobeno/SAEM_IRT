####
# Analysis for Cai (2010) polytomous data. Get libraries for analysis

rm(list=ls(all=TRUE))
library(rlecuyer)		# for rand num generation
library(snow)		# for parallel processing
library(GPArotation)	# for rotations
library(mvnfast)		# for function mvrnorm
library(psych)		# for ML factor analysis

# DrawA generates A factor coefficients w/ eigenanalysis
drawALowerDiag <- function(covT,covTZ,Q,m,n)
{
  Atemp <- matrix(0,Q,Q)
  for (k in 1:(Q-1)) { 
    E 	= k; IQ	= diag(k)
    TT 	= covT[1:k,1:k]
    TZ 	= covTZ[1:k,k] 
    Atemp[k,1:E] <- solve((2*m)*IQ + n*TT)%*%TZ*n
  }
  E 	= Q
  IQ 	= diag(Q)  #  Q X Q
  TT	= covT     #  Q X Q
  Atemp2<-t(sapply(Q:m,function(k) t(solve((2*m)*IQ + n*TT)%*%(covTZ[,k])*n) ))
                        # [Q X Q] * [Q X 1] = Q X 1
  rbind(Atemp,Atemp2)
}
# Number of items (m), schools (n), dimensions (Q)
m <- 100; n <- 100000 ; Q <- 10; ncat <- 4; IQ <- diag(Q)

# Read data
load("C:/Users/egeis_admin/Documents/Rutgers/SAEM_IRT/Poly_J100_N1e+05_Q10_K4.rda")
y.a <- gen.rp
y.a = data.matrix(y.a)

# comps<-princomp(y.a) # 1.182
# which(abs(diff(comps$sdev))>2*sd(diff(comps$sdev)))
# lmcomp<-lm(comps$sdev~c(1:100))
# lmcomp<-lm(comps$sdev[40:60]~c(40:60))
# plot(1:length(comps$sdev),comps$sdev)
# abline(lmcomp)
# MPcomps<-princomp(matrix(as.numeric(cut(rnorm(1000000),
#                                         breaks=qnorm(p = c(0,0.2,0.4,0.6,0.8,1)))),
#                          ncol=100,nrow=10000))
# plot(lm(MPcomps$sdev~c(1:100)))

# Set up n x m x (ncat-1) data array of binary values for
# cumulative option indicators
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
nEM		<-   125	#number of burnin EM cycles < niter
nproc 	<-    10
cl <- makeCluster(nproc,type="SOCK")
clusterSetupRNG(cl, seed = round(runif(6)*1001))
clusterEvalQ(cl,library(mvnfast))
clusterExport(cl,c("Y","n","m","ncat","Q","y.a"))

# wrapper for Z. Returns single Z (instead of mean Z) for est theta.
# wrapper for Z. Returns single Z (instead of mean Z) for est theta.
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
  yL <- matrix(y.a[,j]+1,n,1); yU <- yL + 1
  nn = rep(1:n)
  indL <- cbind(nn,yL); 		indU <- cbind(nn,yU)
  pL <- pnorm( pp[indL] ) ; 	pU <- pnorm( pp[indU] )
  
  Zj <- eta + qnorm( pL + U*(pU - pL)  ) 
  
  return(c(Zj))}
clusterEvalQ(cl,"wrapZ")

# wrapper for snow. m x (ncat-1) matrix is returned for each person
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
  mcol 	<- nmis/(ncat-1)
  
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

# wrapT is wrapper function for drawT
wrapT <- function(i)             
{
  IQ 		<- diag(Q)
  ATZ   	<- t(A)%*%(Z[i,])
  That    	<- BTB_INV%*%ATZ
  ttempi 	<- rmvn(1,That,BTB_INV)
  return(ttempi)}
clusterEvalQ(cl,"wrapT")

# Starting values 
A		= matrix(.2,m,Q)
theta 	= matrix(rnorm(n*Q,0,1),n,Q)
b		= matrix(.5*rnorm(m),m,1);	
d		= matrix(c(-0.25, 0, 0.25),m,ncat-1, byrow=TRUE)
for (L in 1:9) A[L,] <- c(rep(.2,L),rep(0,10-L))

AA=gen.xi[,1:10]; bb=gen.xi[,11]; dd=gen.tau-rowMeans(gen.tau); ttheta=gen.theta
clusterExport(cl,c("A","b","d","theta"))		
covT <- 0; covTZ	<- 0; meanD <- 0; meanB <- 0

# Now iterate
Tstart <- Sys.time()
for (i in 1:niter) {
  if (i < (nEM+1)) alpha <- 1 else alpha <- 1/(i-nEM)
  print(c(i,alpha)); print(i)
  X1 		<- parSapply(cl,1:m,wrapX,simplify=FALSE)
  X2		<- simplify2array(X1, higher=TRUE)	
  X3		<- t(apply(X2,c(1,3),mean))
  meanBMC	<- t(t(rowMeans(X3)))
  meanDMC 	<- t(apply(X3, 1, scale, scale=FALSE)) 
  
  # meanB, meanD are suff stats for item difficulty & thresholds
  meanB	<-  meanB + alpha*(meanBMC-meanB)
  meanD	<-  meanD + alpha*(meanDMC-meanD)
  
  b 	<- -meanB;	d		<- -meanD	
  #b 	<- bb; d <- dd	# to fix b and d to true values
  
  # Get Z to estimate theta
  X1 		<- parSapply(cl,1:m,wrapZ,simplify=FALSE)
  Z		<- simplify2array(X1, higher=TRUE)	
  clusterExport(cl,"Z") 
  
  #compue covariances
  covTMC	<- cov(theta)
  covTZMC	<- cov(theta,Z)
  covT		<- covT  + alpha*(covTMC-covT)
  covTZ		<- covTZ + alpha*(covTZMC-covTZ)
  
  #get ML estimate of A
  Anew		<- drawA(covT,covTZ,Q,m,n)
  print(max(Anew - A))
  A <- Anew
  #A <- AA	# to fix A
  clusterExport(cl,c("A","b","d"))
  
  ##################################
  # Following code is only for diagnostics – doesn’t affect estimation
  #
  RotV <- oblimin(Anew,Tmat=diag(Q),normalize=TRUE,eps=1e-5,maxit=1000)
  A2   = RotV$loadings
  for (ii in 1:Q) if (max(A2[,ii])<0) A2[,ii]=-A2[,ii]
  
  RotT = targetQ(A2,Tmat=diag(Q),Target=AA,normalize=TRUE,eps=1e-5, maxit=1000)
  A2   <- RotT$loadings
  for (ii in 1:Q) if (max(A2[,ii])<0) A2[,ii]=-A2[,ii]
  
  MSE = sqrt((sum((A2-AA)^2))/(Q*m))
  cat (round(MSE,5),"MSE_A  ","\n")
  MSE = (sum(abs(bb-b)))/(Q*m); cat (round(MSE,5),"MSE_b  ","\n")
  MSE = (sum(abs(dd-d)))/(Q*m); cat (round(MSE,5),"MSE_d  ","\n")
  
  dev.off()
  dev.new(width=12, height=10)
  # 10 figures arranged in 2 rows and 5 columns
  par(mfrow=c(2,6))
  for (k in 1:10)
  {
    plot(A2[,k],AA[,k]); abline(a=0,b=1)
  }
  plot(b,bb); abline(a=0,b=1)
  plot(d[,2],dd[,2]); abline(a=0,b=1)
  #
  # End diagnostics
  #########################################################
  ATA 		<- t(A)%*%A
  BTB_INV	<- solve(IQ + ATA)
  clusterExport(cl,"BTB_INV")
  
  #Draw theta
  thetanew	<- t(parSapply(cl,1:n,wrapT))
  print (max(thetanew-theta))
  theta		<- thetanew
  #theta <- ttheta # to fix theta at true values
  clusterExport(cl,"theta")
  
  ########################
  # More diagnostics
  MSE = sqrt((sum((theta-ttheta)^2))/(n*Q))
  cat (round(MSE,5),"MSE_theta","\n")
  cat ("\n")
  # End diagnostics 
  # #####################
  print(c(max(theta), min(theta),max(A2),min(A2),max(Z),min(Z)))
  print(A2[1,1])	
}
Time <- Sys.time()-Tstart;	Time
stopCluster(cl)

