
GoPolyGIFA <- function(rp,init=Init,settings=settings,TargetA=NA) {
  library(rlecuyer)		# for rand num generation
  if (get_os()=="windows") library(snow)			# for parallel processing
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
  Aold1 <- matrix(runif(J*Q),J,Q)
  # Aold2 <- matrix(0,J,Q)
  niter<-ceiling(2/settings$eps)
  nEM<-settings$burnin	#number of burnin EM cycles < niter
  RMwindow<-ceiling(settings$burnin*(0.2))
  estgain=settings$estgain
  eps=settings$eps
  alpha<-GainConstant(settings=settings)
  prior<-list(tmu=settings$tmu,tsigma=settings$tsigma)
  IQ <- diag(Q) 
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
  n1cat = (ncat-1)
  MY = lapply(1:(dim(Y)[1]),function(x) (which(Y[x,,]==9,arr.ind = TRUE)))
  MN = lapply(1:(dim(Y)[1]),function(x) (which(Y[x,,]!=9,arr.ind = TRUE)))
  R  = lapply(1:(dim(Y)[1]),function(x) Y[x,,])
  missList = lapply(MY,function(x) (list(miss = nrow(x),nmis = N*n1cat-nrow(x),mcol = (N*n1cat-nrow(x))/n1cat)))

  # Values for the double truncated WrapZ() function  
  indL<-list();indU<-list()
  for (j in 1:J) {
    yA<-y.a[,j]
    ### This missing handling needs to be revisited... 
    # Best solution is a sample from a different kind of augmented Z... thresholds off?
    if (sum(yA==9)>0) {
      yA[yA==9]<-floor(mean(yA[yA!=9]))
    }
    yL <- matrix(yA+1,N,1)
    indL[[j]] <- cbind(1:N,yL); 		indU[[j]] <- cbind(1:N,yL + 1)
  }
  clusterExport(cl,c("Y","Q","n1cat","N","J","MY","MN","R","missList","indL","indU"),envir=environment())
  clusterEvalQ(cl,c("WrapX","WrapT","WrapZ"))
  
  if (!settings$dbltrunc) {
    X2		<- simplify2array(parSapply(cl,1:J,WrapX,A=A,b=b,d=d,theta=theta,simplify=FALSE), higher=TRUE)	
    Z		<- apply(X2,c(2,3),mean) 
  } else {
    Z 		<- simplify2array(parSapply(cl,1:J,WrapZ,A=A,b=b,d=d,theta=theta,simplify=FALSE), 
                          higher=TRUE)
  }
  covZ  <- cov(Z)
  if (settings$exploreAdim) {
    Q<-Q0<-sum(eigen(covZ)$values>MarcenkoPastur(J=J,N=N))
  } else {
    Q0<-Q
  }
  atemp	<- DrawA(covZ,Q,a=NA)
  A<-atemp$Atemp
  ATA 		<- t(A)%*%A #*4
  BTB_INV	<- solve(IQ + ATA)
  theta	<- t(parSapply(cl,1:N,WrapT,A=A,Z = Z,BTB_INV=BTB_INV,b=b,dbltrunc=settings$dbltrunc))
  meanX	<- 0; meanD <- 0; meanB <- 0
  # alpha <- ifelse(1:niter < (nEM+1),rep(1,niter),1/(1:niter-nEM)^(2/3))
  
  # Now iterate
  if (settings$record) {
    Aiter = array(0,dim=c(J,Q,niter+1))
    Aiter[,,1]<-A
    Viter = array(0,dim=c(Q,niter+1)) 
    Viter[,1]<-eigen(covZ)$values[1:Q]
    Biter = array(0,dim=c(J,niter+1)) 
    Biter[,1]<-b
    Diter = array(0,dim=c(J,n1cat,niter+1))
    Diter[,,1]<-d
  }
  #Tstart<-Sys.time()
  i<-1
  prevA=mat.or.vec(J,Q)
  Avec0<-rep(0,Q)
  itest<-(A-prevA)
  while (max(abs(itest))>eps) {
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
    if (!settings$dbltrunc) {
      Z		<- colMeans(X2)
    } else {
      X1 		<- parSapply(cl,1:J,WrapZ,simplify=FALSE,A=A,b=b,d=d,theta=theta)
      Z		<- simplify2array(X1, higher=TRUE)	
    }
    covZ	<- covZ + alpha[i]*(cov(Z)-covZ)
    # print(i)
    # print("covZ")
    # print(covZ)
    # out<-eigen( covZ-diag(J), symmetric=TRUE)
    # print(out)
    # print(out$vectors[,1:Q]%*%sqrt(diag(out$values[1:Q])))
    #if (i>20) {
      prevA <- A 
      Avec0 <- atemp$Avec[1:Q]
    #}
    if (settings$drawA=="lowertriangular") {
      #compue covariances
      covTMC	<- cov(theta)
      covTZMC	<- cov(theta,Z)
      covT		<- covT  + alpha[i]*(covTMC-covT)
      covTZ		<- covTZ + alpha[i]*(covTZMC-covTZ)
      atemp	<- DrawALowerDiag(covT=covT,covTZ=covTZ,Q,J,N)
    } else if (settings$drawA=="eigen") {
      atemp	<- DrawAEigen(covZ = covZ-diag(J),Q)
    } else {
      atemp	<- DrawA(covZ-diag(J)/n1cat,Q,a=prevA)
    }
    A <- Anew <- atemp$Atemp
    # print("A")
    # print(A)
    ATA 		<- t(A)%*%A
    BTB_INV	<- solve(IQ + ATA)
    # print("BTB_INV")
    # print(BTB_INV)
    theta		<- thetanew	<- t(parSapply(cl,1:N,WrapT,A=A,Z=Z,BTB_INV=BTB_INV,b=b,dbltrunc=settings$dbltrunc))
    # print("theta")
    # print(theta)
    
    if (settings$record) {
      Aiter[,,i+1]<-A
      Viter[,i+1]<-atemp$Avec[1:Q]
      Biter[,i+1]<-b
      Diter[,,i+1]<-d
      if (settings$plots) {
        if (i%%settings$plotiter==0) {
          par(mfrow=c(ceiling((Q+2+ifelse(!is.na(TargetA)[1],Q,0))/5),5),mar=c(4,4,1,1))
          plot(c(1,100*ceiling(i/100)),c(0,2*Q),xlab="iteration",ylab="sqrt(eigenvalues)",type="n")
          abline(h=rowMeans(sqrt(Viter[,max(1,i-20):i])),col=1:J,lwd=0.5,lty=3)
          for (qq in 1:Q) {lines(1:i,sqrt(Viter[qq,1:i]),col=qq)}
          for (qq in 1:Q) {
            plot(c(1,100*ceiling(i/100)),c(-0.1,ifelse(qq==1,2.5,1)),xlab="iteration",ylab="slopes",type="n")
            if (!is.na(TargetA)[1]) {
              abline(h=abs(TargetA[,qq]),col=1,lwd=0.5,lty=3)
            } else {
              abline(h=rowMeans(abs(Aiter[,qq,max(1,i-20):i])),col=1:J,lwd=0.5,lty=3)
            } 
            for (jj in 1:J) {lines(1:i,abs(Aiter[jj,qq,1:i]),col=jj)}
          }
          plot(c(1,100*ceiling(i/100)),c(-3,3),xlab="iteration",ylab="intercepts",type="n")
          abline(h=rowMeans(Biter[,max(1,i-20):i]),col=1:J,lwd=0.5,lty=3)
          for (jj in 1:J) {lines(1:i,Biter[jj,1:i],col=jj)}
          if (!is.na(TargetA)[1]) { # TargetA
            #plot(c(1,100*ceiling(i/100)),c(0,2*Q),xlab="iteration",ylab="eigenvalues",type="n")
            for (ii in 1:Q) {
              if (i > 4) {if (cor(Aold1[,ii],Anew[,ii]) < -.25) Anew[,ii] <- -Anew[,ii]}
            }
            Aold1 <- Anew
            if (atemp$Avec[1]/atemp$Avec[2]>3) {
              Target <- matrix(NA,J,Q)
              Target[TargetA==0] <- 0
              Target <- matrix(Target,J,Q)
              RotT <- targetT(Anew,Tmat=diag(Q),Target=Target,normalize=TRUE,eps=1e-4, maxit=1000)
              #RotV <- Varimax(RotT$loadings,Tmat=diag(Q),normalize=TRUE,eps=1e-4,maxit=100)
              A2 <- RotT$loadings
            } else {
              RotT <- targetT(Anew,Tmat=diag(Q),Target=TargetA,normalize=TRUE,eps=1e-4, maxit=1000)
              RotV <- Varimax(RotT$loadings,Tmat=diag(Q),normalize=TRUE,eps=1e-4,maxit=100)
              A2 <- RotV$loadings
            }
            for (ii in 1:Q) {
              if (abs(min(A2[,ii])) > max(A2[,ii])) A2[,ii]=-A2[,ii]
              plot(TargetA[,ii], A2[,ii],ylab=paste("Est A",ii),
                   xlab = paste0("True A",ii)); abline(a=0,b=1)
            }
          }
        }
      }
    }
    i<-i+1
    if (tolower(settings$converge)=="a"|grepl("slop",tolower(settings$converge))) {
      itest<-(A-prevA)
    } else if (grepl("eig.+val",tolower(settings$converge))|grepl("e?(.+)val",tolower(settings$converge))) {
      itest<-(atemp$Avec[1:Q] - Avec0)
    } else {
      itest<-(A-prevA)
    }
  }
  # Time <- Sys.time()-Tstart;	print(paste(settings$cores,"processors:",Time))
  
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