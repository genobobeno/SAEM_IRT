####
# Analysis for Cai (2010) polytomous data
# Get libraries for analysis
#####

rm(list=ls(all=TRUE))
library(rlecuyer)        # for rand num generation
library(snow)            # for parallel processing
library(GPArotation)    # for rotations
library(mvnfast)        # for function mvrnorm
library(psych)            # for ML factor analysis

# DrawA generates A factor coefficients w/ eigenanalysis
drawA <- function(covT,covTZ,Q,J,n) {
  Abtemp <- matrix(0,J,Q)
  # TT      <- covT
  for (k in 1:J)
  { 
    if (k<=10) {    
      E = k
      IQ    = diag(k)
      TT     = covT[1:E,1:E]
      TZ     = covTZ[1:E,k] 
    } else {    
      E     = Q
      IQ     = diag(Q)
      TT    = covT
      TZ    = covTZ[,k]
    }
    Abtemp[k,1:E] <- solve(n*TT + 0*IQ)%*%TZ*n
  }
  return(Abtemp)
}
# Number of items (J), schools (n), dimensions (Q)
J <- 100; n <- 100000 ; Q <- 10; ncat <- 4; IQ <- diag(Q)

# Read data
load("C:\\Users\\gregory.camilli\\Desktop\\Eugene\\Poly_J100_N1e+05_Q10_K4.rda")
y.a <- gen.rp
y.a = data.matrix(y.a)
# Set up n x J x (ncat-1) data array of binary values for
# cumulative option indicators
y.b = array(NA,c(J,ncat-1,n))
for (i in 1:n) 
  for (j in 1:J)
    # for missing responses, item propensities = 9 
    if (y.a[i,j]==9)    {y.b[j,,i] <- rep(9,(ncat-1))} else
    {
      for (k in 0:(ncat-2))
      {        
        if (y.a[i,j]<=k)     {y.b[j,(k+1),i] <- 0} 
        else             {y.b[j,(k+1),i] <- 1} 
      }}
Y = y.b


niter     <-  15    
nEM        <-  15    #number of burnin EM cycles < niter
nproc     <-  10

cl <- makeCluster(nproc,type="SOCK")
clusterSetupRNG(cl, seed = round(runif(6)*1001))
clusterEvalQ(cl,library(mvnfast))
clusterExport(cl,c("Y","n","J","ncat","Q"))

# wrapper for snow called by drawX. For each person,
# J x (ncat-1) matrix is returned for each person
wrapX <- function(j) {
  
  # ker is J x (ncat-1) matrix of item kernels
  ker <- apply(theta%*%matrix(A[j,],Q,1)-b[j],1,
               function(x) x - c(d[j,]))
  
  pn     <- pnorm(-ker)
  r      = Y[j,,]
  Xi    = matrix(NA,ncat-1,n)
  
  # Count missing and nonmissing values (9)
  miss         <- sum(r==9)
  nmis         <- n*(ncat-1)- miss
  mcol     <- nmis/(ncat-1)
  
  # extract missing & nonmissing reponses
  my <- which(r == 9,arr.ind=TRUE)
  mn <- which(r != 9,arr.ind=TRUE)
  
  # Generate missing option propensities
  Xi[my] <- rnorm(miss)
  
  # Generate nonmissing option propensities
  r1         <- matrix(r[mn]      ,ncat-1,mcol)
  U          <- matrix(runif(nmis),ncat-1,mcol)
  P        <- matrix(pn[mn]     ,ncat-1,mcol)
  Xi[mn]     <- qnorm( r1*U + P*(r1 + U - 2*r1*U) )
  Xi        <- ker + Xi
  return(Xi)}

clusterEvalQ(cl,"wrapX")

# wrapT is wrapper function for drawT
wrapT <- function(i)             
{
  IQ         <- diag(Q)
  ATZ       <- t(A)%*%(Z[i,] + b)*(ncat-1)
  That        <- BTB_INV%*%ATZ
  ttempi     <- rmvn(1,That,BTB_INV)
  return(ttempi)}

clusterEvalQ(cl,"wrapT")




# Starting values 
theta     = matrix(rnorm(n*Q,0,1),n,Q)
clusterExport(cl,"theta")

A        = matrix(.25*(runif(J*Q)-.5),J,Q)
b        = matrix(rnorm(J),J,1);    
d        = matrix(c(-0.50, -0.25, 0.25, 0.50),J,ncat-1, byrow=TRUE)
#clusterExport(cl,c("A","b","d","theta"))

AA=gen.xi[,1:10]
bb=gen.xi[,11]
dd=gen.tau-rowMeans(gen.tau)
ttheta=gen.theta
clusterExport(cl,c("A","b","d","theta"))        
covT <- 0; covTZ    <- 0; meanD <- 0; meanB <- 0


# Now iterate
Tstart <- Sys.time()
for (i in 1:niter) {
  if (i < (nEM+1)) alpha <- 1 else alpha <- 1/(i-nEM)
  print(c(i,alpha)); print(i)
  X1         <- parSapply(cl,1:J,wrapX,simplify=FALSE)
  X2        <- simplify2array(X1, higher=TRUE)    
  X3        <- t(apply(X2,c(1,3),mean))
  meanBMC    <- t(t(rowMeans(X3)))
  meanDMC     <- t(apply(X3, 1, scale, scale=FALSE)) 
  
  # meanB is sufficient stat for overall item difficulty
  # meanD is sufficient stat for cumulative thresholds
  meanB    <-  meanB + alpha*(meanBMC-meanB)
  meanD    <-  meanD + alpha*(meanDMC-meanD)
  
  b     <- -meanB;    d        <- -meanD    
  #b     <- bb; d <- dd    # to fiex b and d to true values
  
  # use avg(X) to estimate A
  Z        <- colMeans(X2); clusterExport(cl,"Z")
  covTMC    <- cov(theta)
  covTZMC    <- cov(theta,Z)
  
  covT        <- covT  + alpha*(covTMC-covT)
  covTZ        <- covTZ + alpha*(covTZMC-covTZ)
  
  #get ML estimate of A
  Anew        <- drawA(covT,covTZ,Q,J,n)
  print(max(Anew - A))
  A <- Anew
  #A <- AA    # to fix A
  
  clusterExport(cl,c("A","b","d"))
  
  
  
  
  
  
  ##################################
  # Following code is only for diagnostics – doesn’t affect estimation
  #
  RotV <- oblimin(Anew,Tmat=diag(Q),normalize=FALSE,eps=1e-5,maxit=1000)
  Anew = RotV$loadings
  for (ii in 1:Q) if (sum(Anew[,ii])<0) Anew[,ii]=-Anew[,ii]
  
  RotT = targetT(Anew, Tmat=diag(Q), Target=AA, normalize=FALSE, 
                 eps=1e-6, maxit=1000)
  Anew <- RotT$loadings
  for (ii in 1:Q) if (sum(Anew[,ii])<0) Anew[,ii]=-Anew[,ii]
  
  MSE = sqrt((sum((Anew-AA)^2))/(Q*m))
  cat (round(MSE,5),"MSE_A  ","\n")
  MSE = (sum(abs(bb-b)))/(Q*m); cat (round(MSE,5),"MSE_b  ","\n")
  MSE = (sum(abs(dd-d)))/(Q*m); cat (round(MSE,5),"MSE_d  ","\n")
  
  dev.off()
  dev.new(width=12, height=10)
  # 10 figures arranged in 2 rows and 5 columns
  par(mfrow=c(2,6))
  for (k in 1:10)
  {
    plot(Anew[,k],AA[,k])
    abline(a=0,b=1)
  }
  plot(b,bb)
  abline(a=0,b=1)
  plot(d[,2],dd[,2])
  abline(a=0,b=1)
  #
  # End diagnostics
  #########################################################
  
  ATA         <- t(A)%*%A*(ncat-1)
  BTB_INV    <- solve(IQ + ATA)
  clusterExport(cl,"BTB_INV")
  
  #Draw theta
  thetanew    <- t(parSapply(cl,1:n,wrapT))
  print(max(thetanew-theta))
  theta       <- thetanew
  theta <- ttheta # to fix theta at true values
  clusterExport(cl,"theta")
  
  ########################
  # More diagnostics
  MSE = sqrt((sum((theta-ttheta)^2))/(n*Q))
  cat (round(MSE,5),"MSE_theta","\n")
  cat ("\n")
  # End diagnostics 
  # ######################
  
  print(c(max(theta), min(theta),max(Anew),min(Anew),max(Z),min(Z)))    
}
Time <- Sys.time()-Tstart;    Time
stopCluster(cl)
