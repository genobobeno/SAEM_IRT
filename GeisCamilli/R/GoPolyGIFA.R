
GoPolyGIFA <- function(rp,init=Init,settings=settings,TargetA=NA) {
  library(rlecuyer)		# for rand num generation
  library(snow)			# for parallel processing
  library(GPArotation)	# for rotations
  library(mvnfast)		# for function mvrnorm
  library(psych)			# for ML factor analysis
  library(foreach)
  clusterEvalQ(cl,library(mvnfast))
  
  J=ncol(rp);N=nrow(rp);Q=settings$Adim
  # Starting values 
  theta = init$THat
  A		= init$XI[,1:settings$Adim]
  b		= init$XI[,ncol(init$XI)-settings$guess]	
  d		= init$D
  ncat <- settings$ncat
  #J=20;N=1000;Q=2;K=4
  niter<-ceiling(2/settings$eps)
  nEM<-settings$burnin	#number of burnin EM cycles < niter
  RMwindow<-ceiling(settings$burnin*(0.2))
  estgain=settings$estgain
  eps=settings$eps
  # nproc 	<- 2
  alpha<-GainConstant(settings=settings)
  prior<-list(tmu=settings$tmu,tsigma=settings$tsigma)
    # c(rep(1,nEM-RMwindow),
    #        runif(rep(1,ceiling(RMwindow/2)),min = 1.0/(1:(RMwindow/2))^estgain, max = 1.0),1.0,
    #        runif(rep(1,ceiling(RMwindow/2)),min = 1.0/((RMwindow/2+1):RMwindow)^estgain, max = 1.0/(1:(RMwindow/2))^estgain),
    #        1/(ceiling(RMwindow/2):(niter-nEM+ceiling(RMwindow/2)))^estgain)  

  # csvfile = paste0("Poly_Q",Q,"_K",K,"_N",N,"_ResponseData.csv")
  # datafile = paste0("Poly_J",J,"_N",N,"_Q",Q,"_K",K,".rds")
  # GEN.DATA = readRDS(datafile)
  IQ <- diag(Q) 
  #TrueTheta = data.matrix(read.csv(csvfile,header=TRUE)[,1:Q])
  # Read data
  # y.a = scan("form3_social.dat", what = "numeric",sep = "\n")
  # y.a = data.matrix(read.csv(csvfile,header=TRUE)[,-(1:Q)]-1)
  # y.a = data.matrix(do.call(rbind,lapply(strsplit(y.a,"\\s"), function(x) as.numeric(x[-1]))))
  # Set up n x m x (ncat-1) data array of binary values for
  y.a = data.matrix(rp)
  # cumulative option indicators
  y.b = array(NA,c(J,ncat-1,N))
  for (i in 1:N) { 
    for (j in 1:J) {
      # for missing responses, item propensities = 9 
      if (y.a[i,j]==9)	{
        y.b[j,,i] <- rep(9,(ncat-1))
      } else {
        for (k in 0:(ncat-2))
          {		
            y.b[j,(k+1),i] <- ifelse(y.a[i,j]<=k,0,1) 
          }
      }
    }
  }
  Y <- y.b
  # nproc 	<- 24
  #  procs = nproc = parallel::detectCores()
  # for (procs in 1:nproc) {
  # if (settings$parallel) {
  #   parallelCluster <- parallel::makeCluster(settings$cores,type="SOCK")
  #   print(parallelCluster)
  # clusterSetupRNG(parallelCluster, seed = round(runif(6)*1001))
  
  # wrapper for snowfall called by drawX. For each person,
  # m x (ncat-1) matrix is returned for each person
  n1cat = (ncat-1)
  MY = lapply(1:(dim(Y)[1]),function(x) (which(Y[x,,]==9,arr.ind = TRUE)))
  MN = lapply(1:(dim(Y)[1]),function(x) (which(Y[x,,]!=9,arr.ind = TRUE)))
  R  = lapply(1:(dim(Y)[1]),function(x) Y[x,,])
  missList = lapply(MY,function(x) (list(miss = nrow(x),nmis = N*n1cat-nrow(x),mcol = (N*n1cat-nrow(x))/n1cat)))
  
  refList <- list(MY=MY,MN=MN,R=R,missList=missList)
  # wrapper for snowfall called by drawX. For each person,
  # m x (ncat-1) matrix is returned for each person
  # wrapT is wrapper function for drawT
  clusterExport(cl,c("Y","Q","n1cat","N","J","MY","MN","R","missList"),envir=environment())
  clusterEvalQ(cl,c("WrapX","WrapT"))
  
  X2		<- simplify2array(parSapply(cl,1:J,WrapX,A=A,b=b,d=d,theta=theta,simplify=FALSE), higher=TRUE)	
  Z		<- apply(X2,c(2,3),mean) 
  covZ  <- cov(Z)
  A	<- DrawA(covZ,Q,a=NA)
  ATA 		<- t(A)%*%A #*4
  BTB_INV	<- solve(IQ + ATA)
  theta	<- t(parSapply(cl,1:N,WrapT,A=A,Z = Z,BTB_INV=BTB_INV,b=b))
  meanX	<- 0; meanD <- 0; meanB <- 0
  # alpha <- ifelse(1:niter < (nEM+1),rep(1,niter),1/(1:niter-nEM)^(2/3))
  
  # Now iterate
  if (settings$record) {
    Aiter = array(0,dim=c(J,Q,niter+1))
    Aiter[,,1]<-A
  }
  Tstart <- Sys.time()
  i<-1
  prevA=mat.or.vec(J,Q)
  while (max(A-prevA)>eps) {
    #if (i%%100==0) cat(nproc,i,"|")
    if (i%%10==1) cat(".")
    if (i%%100==1) cat(":")
    if (i%%500==1) cat("\n",i,"\t : ")
    X2		<- simplify2array(parSapply(cl,1:J,WrapX,simplify=FALSE,A=A,b=b,d=d,theta=theta), higher=TRUE)	
    X3		<- t(apply(X2,c(1,3),mean))
    meanB <- meanB + alpha[i]*(t(t(rowMeans(X3)))-meanB)
    meanD <- meanD + alpha[i]*(t(apply(X3, 1, scale, scale=FALSE)) - meanD)
    b         <- -meanB;    d        <- -meanD   
    #Z		<- apply(X2,c(2,3),mean)
    Z		<- colMeans(X2)
    covZ	<- covZ + alpha[i]*(cov(Z)-covZ)
    if (i>20) prevA <- A
    A <- Anew	<- DrawA(covZ-diag(J)/n1cat,Q,a=prevA)
    ATA 		<- t(A)%*%A
    BTB_INV	<- solve(IQ + ATA)
    theta		<- thetanew	<- t(parSapply(cl,1:N,WrapT,A=A,Z=Z,BTB_INV=BTB_INV,b=b))
    if (settings$record) Aiter[,,i+1]<-A
    if (settings$plots) {
      if (i%%100==0) {
        par(mfrow=c(1,Q))
        for (qq in 1:Q) {
          plot(c(1,niter),c(-2,2))
          abline(h=GEN.DATA$XI[,qq],col=1:J)
          for (jj in 1:J) {lines(1:i,Aiter[jj,qq,1:i],col=jj)}
        }
      }
    }
    i<-i+1
  }
  Time <- Sys.time()-Tstart;	print(paste(settings$cores,"processors:",Time))
  
  if (settings$Adim>1 & !is.na(TargetA)[1]) {
    if (length(TargetA)==length(A)) {
      # get eigenvectors for A
      out    <- eigen(A[,1:Q]%*%t(A[,1:Q]))
      Avec   <- out$values ; Aload     <- out$vectors
      TA     <- Aload[,1:Q]%*%sqrt(diag(Avec[1:Q]))
      # set up bifactor target matrix
      # RTS<-permn(1:settings$Adim)
      # rtest<-c()
      # for (i in 1:length(RTS)) {
      Target = matrix(NA,J,Q)
      Target[TargetA==0] = 0
      Target = matrix(Target,J,Q)
      #   rtest<-c(rtest,sum(abs(TargetA[,RTS[[i]]]-abs(pstT(A, Tmat=diag(Q), W=WR, Target=as.matrix(gen.xi[,RTS[[i]]]), normalize=TRUE, eps=1e-8, maxit=1000)$loadings))))
      #   # }    
      #   # Fctr<-RTS[[which.min(rtest)]]
      Rtest<-c()
      TestMat = diag(Q)
      RMAT = list()
      SRot<-sample(1:Q,8,replace = TRUE)
      SAng<-sample(2*pi/(0:360),8)
      for (i in 1:8) {
        if (SRot[i]!=SRot[i%%8+1]) {
          TestMat[SRot[i],SRot[i]]<-sin(SAng[i])*TestMat[SRot[i],SRot[i]]
          TestMat[SRot[i%%8+1],SRot[i%%8+1]]<-sin(SAng[i])*TestMat[SRot[i%%8+1],SRot[i%%8+1]]
          TestMat[SRot[i],SRot[i%%8+1]]<-ifelse(TestMat[SRot[i],SRot[i%%8+1]]==0,cos(SAng[i]),cos(SAng[i])*TestMat[SRot[i],SRot[i%%8+1]])
          TestMat[SRot[i%%8+1],SRot[i]]<-ifelse(TestMat[SRot[i%%8+1],SRot[i]]==0,cos(-SAng[i]),cos(-SAng[i])*TestMat[SRot[i%%8+1],SRot[i]])
        }
        RMAT[[i]]<-TestMat
        RotA = targetQ(TA, Tmat=TestMat, Target=Target, normalize=FALSE, 
                       eps=1e-4, maxit=10000)
        Rtest<-c(Rtest,sd(RotA$loadings[!is.na(Target) & Target==0]))
      }
      AR<- targetQ(TA, Tmat=RMAT[[which.min(Rtest)]], 
                     Target=Target, normalize=FALSE, 
                     eps=1e-4, maxit=10000)
      C<-NA;W<-NA
      THAT<-GetThetaHat(aa=A,bb=b,cc=C,rp=rp,tHat=theta,zHat=Z,w=W,
                        prior=prior,setting=settings,R=AR,
                        TAU=d+matrix(rep(b,ncat-1),length(b),ncat-1))

    } else {AR<-A}
  } else {
    AR <- A
  }
  #summaryRprof()
  C<-NA
  list(xi=cbind(AR,b),A=A,AR=AR,B=b,C=C,TAU=d+matrix(rep(b,ncat-1),length(b),ncat-1),RP=rp,
       xiError=NA,iError=NA,oError=NA,gain=alpha,EZ=Z,EZZ=covZ,
       That=theta,Tmap=NA,Tmaprot=NA,TRmap=NA,Trot=NA,
       EmpSE=NA,ThetaFix=NA,settings=settings)
}

  # for (i in 1:length(RTS)) {
  #   # A is A_gen
  #   # B is estimated loading matrix
  #   # W is a weight matrix. The rotation target is the bifactor 0’s
  #   # pstT is partially specified target orthogonal rotation
  #   WR <- matrix(0,J,settings$Adim)
  #   WR[which(gen.xi[,1:settings$Adim]==0)] <- 1
  #   Tmat <- matrix(-1,settings$Adim,settings$Adim)
  #   Tmat[1,1] <-  1
  #   rtest<-c(rtest,sum(abs(gen.xi[,RTS[[i]]]-abs(pstT(A, Tmat=Tmat, W=WR, Target=as.matrix(gen.xi[,RTS[[i]]]), normalize=TRUE, eps=1e-8, maxit=1000)$loadings))))
  # }    
  # Fctr<-RTS[[which.min(rtest)]]
  # AR <- pstT(A, Tmat=Tmat, W=WR, Target=as.matrix(gen.xi[,Fctr]), normalize=TRUE, eps=1e-8, maxit=1000)
  # rits<-1
  # while (min(apply((AR$loadings>0)+0,2,mean))<0.5) {
  #   rots<-rep(-1,settings$Adim^2)
  #   sr<-sample(1:settings$Adim^2)
  #   rots[sr[1:(rits%%(settings$Adim^2)+1)]]<-1
  #   Tmat<-matrix(rots,settings$Adim,settings$Adim)
  #   AR <- pstT(A, Tmat=Tmat, W=WR, Target=as.matrix(gen.xi[,Fctr]), normalize=TRUE, eps=1e-8, maxit=1000)
  #   rits<-rits+1
  # }
  # AR$APermute<-Fctr

  ############## GREG's reconstruction
  # B.TMS <- A
  # B.TMS <- data.matrix(B.TMS)
  # #B  <- Varimax(B.TMS,Tmat=diag(ncol(B.TMS)),normalize=TRUE,eps=1e-5, maxit=10000)
  # B <- oblimin(B.TMS,Tmat=diag(ncol(B.TMS)),normalize=TRUE,eps=1e-5, maxit=10000,gam=0)
  # mla <- fa(covZ,nfactors=Q,covar=TRUE,rotate="oblimin",fm="ml")
  # covar     <- covZ - diag(m)
  # out         <- eigen( covar, symmetric=TRUE)
  # Avec     <- out$values
  
  # par(mfrow=c(2,Q/2))
  # for (qq in 1:Q) {
  #   RotA$loadings[,qq]<-ifelse(sum(RotA$loadings[,qq])<0,-1,1)*RotA$loadings[,qq]
  #   plot(cbind(GEN.DATA$XI[,qq],RotA$loadings[,qq]),xlab="Actual",ylab="Fit",main=paste0("Q",qq))
  #   abline(0,1)
  #   cor(-RotA$loadings[,qq],GEN.DATA$XI[,qq])
  #   lm(-RotA$loadings[,qq]~GEN.DATA$XI[,qq])
  # }