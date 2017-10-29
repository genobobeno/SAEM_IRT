####
# Analysis for Cai (2010) polytomous data
# Get libraries for analysis
#####
setwd("~/Documents/Documents/Rutgers_GSE/Camilli/ParSAEM/SAEM_IRT/")
install.packages("rlecuyer",dependencies = T)
install.packages("parallel",dependencies = T)
install.packages("GPArotation",dependencies = T)
install.packages("mvnfast",dependencies = T)
install.packages("psych",dependencies = T)
install.packages("foreach",dependencies = T)
install.packages("doParallel",dependencies = T)
install.packages("Rcpp",dependencies = T)
install.packages("RcppArmadillo",dependencies = T)
install.packages("doRNG",dependencies = T)
install.packages("compiler",dependencies = T)
install.packages("devtools",dependencies = T)
if ("lineprof" %in% rownames(installed.packages()) == FALSE) devtools::install_github("hadley/lineprof")

#Rprof()

library(devtools)
library(rlecuyer)		# for rand num generation
library(parallel)			# for parallel processing
library(GPArotation)	# for rotations
library(mvnfast)		# for function mvrnorm
library(psych)			# for ML factor analysis
suppressMessages(library(foreach))
library(doParallel)
library(doRNG)
library(compiler)
# library(Rcpp)
# library(RcppArmadillo)
enableJIT(1)

J=200;N=10000;Q=10;K=4
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
#y.a = data.matrix(do.call(rbind,lapply(strsplit(y.a,"\\s"), function(x) as.numeric(x[-1]))))
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
    if (y.a[i,j]==9)	{y.b[j,,i] <- 9} else
    {
      for (k in 0:(ncat-2))
      {		
        y.b[j,(k+1),i] <- ifelse(y.a[i,j]<=k, 0,1) 
      }}
Y = y.b

#for (nproc in 1:8) {
nproc = parallel::detectCores()
parallelCluster <- parallel::makeCluster(nproc,type="FORK")
print(parallelCluster)
registerDoParallel(parallelCluster)
generateRandom <- function(rng='default',n) {
  setRNG(rng)
  runif(n)
}

clusterEvalQ(parallelCluster,library(mvnfast))
clusterSetRNGStream(parallelCluster,iseed = round(runif(6)*1001))

# clusterSetupRNG(parallelCluster, seed = round(runif(6)*1001))

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
{ return(rmvn(1,BTB_INV%*%(t(A)%*%(Z[i,] + b)),BTB_INV))} #*4
clusterExport(parallelCluster,c("Y","n","m","ncat","Q","MY","MN","R","n1cat","missList","wrapX","wrapT"))

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

alpha<-c(rep(1,nEM-RMwindow),
        runif(rep(1,ceiling(RMwindow/2)),min = 1.0/(1:(RMwindow/2))^estgain, max = 1.0),1.0,
        runif(rep(1,ceiling(RMwindow/2)),min = 1.0/((RMwindow/2+1):RMwindow)^estgain, max = 1.0/(1:(RMwindow/2))^estgain),
        1/(ceiling(RMwindow/2):(niter-nEM+ceiling(RMwindow/2)))^estgain)  

# Now iterate
Aiter = array(GEN.DATA$XI[,1:Q],dim=c(m,Q,niter+1))
Aiter[,,1]<-A
Tstart <- Sys.time()
for (i in 1:niter) {
  cat(c(".",":","\n")[i%%c(10,100,1000)==0])
  X2		<- simplify2array(parSapply(parallelCluster,1:m,wrapX,simplify=FALSE,A=A,b=b,d=d,theta=theta), higher=TRUE)	
  X3		<- t(apply(X2,c(1,3),mean))
  meanB <- meanB + alpha[i]*(t(t(rowMeans(X3)))-meanB)
  meanD <- meanD + alpha[i]*(t(apply(X3, 1, scale, scale=FALSE)) - meanD)
  b         <- -meanB;    d        <- -meanD   
  Z		<- apply(X2,c(2,3),mean)
  covZ	<- covZ + alpha[i]*(cov(Z)-covZ)
  if (i<20) {prevA <- NA} else {prevA <- A}
  A <- Anew		<- drawA(covZ-diag(m)/n1cat,Q,m,a=prevA)
  ATA 		<- t(A)%*%A #*4
  BTB_INV	<- solve(IQ + ATA)
  theta		<- thetanew	<- t(parSapply(parallelCluster,1:n,wrapT,A=A,Z=Z,BTB_INV=BTB_INV,b=b))
  Aiter[,,i+1]<-A
  if (i%%100==0) {
    par(mfrow=c(1,Q))
    for (qq in 1:Q) {
      plot(c(1,niter),c(-2,2))
      abline(h=GEN.DATA$XI[,qq],col=1:m)
      for (jj in 1:m) {lines(1:i,Aiter[jj,qq,1:i],col=jj)}
    }
  }
}
Time <- Sys.time()-Tstart;	print(paste(nproc,"processors:",Time))
if(!is.null(parallelCluster)) {
  parallel::stopCluster(parallelCluster)
  parallelCluster <- c()
}
#summaryRprof()



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

# round(cbind(-RotA$loadings,GEN.DATA$XI[,1:Q]),2)
# plot(cbind(-RotA$loadings[,1],GEN.DATA$XI[,1]))
# abline(0,1)
# cor(-RotA$loadings[,1],GEN.DATA$XI[,1])
# lm(-RotA$loadings[,1]~GEN.DATA$XI[,1])
# plot(cbind(RotA$loadings[,2],GEN.DATA$XI[,2]))
# abline(0,1)
# cor(RotA$loadings[,2],GEN.DATA$XI[,2])
# lm(RotA$loadings[,2]~GEN.DATA$XI[,2])




Target = matrix(NA,J,Q)
for (i in 2:Q) Target[i*J/Q+1:(J/Q),i] = 0
WR <- matrix(0,J,2)
WR[which(GEN.DATA$XI[,1:Q]==0)] <- 1

AR<-targetT(A, Tmat=diag(ncol(A)), Target=Target, normalize=FALSE, eps=1e-5, maxit=1000)

cor(as.vector(AR$loadings),as.vector(GEN.DATA$XI[,1:Q]))

ARList = list()
for (i in 1:(2*Q)) {
  Tmat <- matrix(-1,Q,Q)
  Tmat[i] <-  1
  AR1 <- pstT(A, Tmat=Tmat, W=WR, Target=GEN.DATA$XI[,1:Q], normalize=FALSE, eps=1e-8, maxit=1000)
  AR2 <- bifactorT(L = A, Tmat=Tmat, normalize=TRUE, eps=1e-8, maxit=1000)
  AR3 <- targetT(A, Tmat=Tmat, Target=GEN.DATA$XI[,1:Q], normalize=FALSE, eps=1e-8, maxit=1000)
  AR4 <- oblimin(L = A, Tmat=Tmat, gam=0, normalize=FALSE, eps=1e-8, maxit=1000)
  AR = c(cor(as.vector(AR1$loadings),as.vector(GEN.DATA$XI[,1:Q])),
         cor(as.vector(AR2$loadings),as.vector(GEN.DATA$XI[,1:Q])),
         cor(as.vector(AR3$loadings),as.vector(GEN.DATA$XI[,1:Q])),
         cor(as.vector(AR4$loadings),as.vector(GEN.DATA$XI[,1:Q])))
  ARList[[i]] <- c(which.max(AR),max(AR))
}  

ARList
par(mfrow=c(2,1))
plot(as.vector(GEN.DATA$XI[,1:2]),as.vector(AR1$loadings),main="All Slopes",xlab="Generated Slopes",ylab="Reconstructed Slopes")
abline(0,1)
plot(rbind(AR1$loadings,GEN.DATA$XI[,1:2]),type="n",main="Arrows",xlab="A1",ylab="A2")
arrows(x0=GEN.DATA$XI[,1],y0 = GEN.DATA$XI[,2],x1 = AR1$loadings[,1],y1 = AR1$loadings[,2],col = rep(c(1,2),c(15,15)))
points(rbind(AR1$loadings,GEN.DATA$XI[,1:2]),col=rep(c(1,2),c(30,30)))



#       if (settings$rmethod=="targetT" & tolower(settings$fm)!="pca") {
#         for (i in 1:length(RTS)) {
#           ATest<-Apar[,RTS[[i]]]
#           rtest<-c(rtest,sum(abs(ATest-targetT(A, Tmat=diag(ncol(A)), Target=ATest, normalize=FALSE, eps=1e-5, maxit=1000)$loadings)))
#           print("Rotating A, permuting:")
#           print(RTS[[i]])
#           print("Generated:")
#           print(ATest)
#           print("Rotated A:")
#           print(targetT(A, Tmat=diag(ncol(A)), Target=ATest, normalize=FALSE, eps=1e-5, maxit=1000)$loadings)
#           print(rtest)
#         }    
#         Fctr<-RTS[[which.min(rtest)]]
#         AR<-targetT(A, Tmat=diag(ncol(A)), Target=Apar[,Fctr], normalize=FALSE, eps=1e-5, maxit=1000)
#         AR$APermute<-Fctr
#         A<-AR$loadings
#         #Rotate Theta via %*%t(Th)
#       } else if (settings$rmethod=="pstT" & tolower(settings$fm)!="pca") {
#         print("starting pstT rotation")
#         for (i in 1:length(RTS)) {
#           # A is A_gen
#           # B is estimated loading matrix
#           # W is a weight matrix. The rotation target is the bifactor 0’s
#           # pstT is partially specified target orthogonal rotation
#           ATest<-Apar[,RTS[[i]]]
#           WR <- matrix(0,J,settings$Adim)
#           WR[which(gen.xi[,1:settings$Adim]==0)] <- 1
#           Tmat <- matrix(-1,settings$Adim,settings$Adim)
#           Tmat[1,1] <-  1
#           print("Created WR and Tmat")
#           print(WR)
#           print(Tmat)
#           rtest<-c(rtest,sum(abs(ATest-abs(pstT(A, Tmat=Tmat, W=WR, Target=as.matrix(ATest), normalize=TRUE, eps=1e-8, maxit=1000)$loadings))))
#           # examine mean square residual of loadings. Not too shabby.
#           print("Partially specified target rotation, MSE:")
#         }    
#         Fctr<-RTS[[which.min(rtest)]]
#         AR <- pstT(A, Tmat=Tmat, W=WR, Target=as.matrix(Apar[,Fctr]), normalize=TRUE, eps=1e-8, maxit=1000)
#         rits<-1
#         while (min(apply((AR$loadings>0)+0,2,mean))<0.5) {
#           rots<-rep(-1,settings$Adim^2)
#           sr<-sample(1:settings$Adim^2)
#           rots[sr[1:(rits%%(settings$Adim^2)+1)]]<-1
#           Tmat<-matrix(rots,settings$Adim,settings$Adim)
#           AR <- pstT(A, Tmat=Tmat, W=WR, Target=as.matrix(Apar[,Fctr]), normalize=TRUE, eps=1e-8, maxit=1000)
#           rits<-rits+1
#         }