RestartChainForSE <- function(FitDATA,IT=2000,estgain=NA,thinA=NA,thinB=NA,TargetA=TRUE) {
  # FitDATA<-list(RP=rp,xi=cbind(ARot,b),A=A,AR=AR,B=b,C=C,
  #               tau=d+matrix(rep(b,ncat-1),length(b),ncat-1),
  #               EZ=Z,EZZ=covZ,settings=settings,That=THAT$THETA,Trot=TROT)
  cat("\n")
  library(rlecuyer)		# for rand num generation
  if (get_os()=="windows") library(snow)			# for parallel processing
  library(GPArotation)	# for rotations
  library(mvnfast)		# for function mvrnorm
  library(psych)			# for ML factor analysis
  library(foreach)
  if (FitDATA$settings$parallel) { # & settings$cores>1) {
    if (get_os()=="windows") {
      cl <<- parallel::makePSOCKcluster(FitDATA$settings$cores)
      if (is.na(FitDATA$settings$ncat)[1] | FitDATA$settings$ncat==2) {
        clusterCall(cl, cppInit)
      }
      parallel::clusterSetRNGStream(cl, iseed = round(runif(6)*1001))
      #parallel::setDefaultCluster(cl)
      prl<-TRUE
    } else if (get_os()=="linux") {
      cl <<- parallel::makeCluster(FitDATA$settings$cores,type="FORK",outfile="ClusterOutput.txt")
      doParallel::registerDoParallel(cl)
      parallel::clusterSetRNGStream(cl,iseed = round(runif(6)*1001))
      #parallel::setDefaultCluster(cl)
      prl<-TRUE
    } else {
      cat("\nI have only been testing for 'linux' and 'windows' operating systems. 
        Please feel free to edit the code in StartParallel() to add yours.\n")
      prl<-FALSE
      FitDATA$settings$cores<-1
    } 
  } else {
    prl<-FALSE
    FitDATA$settings$cores<-1
  } 
  if (prl) {
    clusterEvalQ(cl,library(mvnfast))
  }
  
  print(paste("Starting",IT,"iterations of the MCMC Chain"))
  rp=FitDATA$RP
  settings=FitDATA$settings
  ncat = settings$ncat
  Q<-settings$Adim
  if (Q>1) {
    if (!is.na(FitDATA$AR)) {
      A=as.matrix(FitDATA$AR$loadings)
      THat=as.matrix(FitDATA$Trot[,1:settings$Adim])
    } else {
      A=as.matrix(FitDATA$A)
      THat=as.matrix(FitDATA$That[,1:settings$Adim])
    }
    B=FitDATA$B
    D=FitDATA$tau - matrix(rep(B,ncat-1),length(B),ncat-1)
    if (settings$fm %in% c("new","pca")) settings$fm <- "camilli"
  } else {
    A=FitDATA$A
    B=FitDATA$B
    D=FitDATA$tau - matrix(rep(B,ncat-1),length(B),ncat-1)
    THat=FitDATA$That[,1]
  }
  # Alines<-A
  # Blines<-B
  prior<-list(tmu=settings$tmu,tsigma=settings$tsigma)
  J<-ncol(rp)
  N<-nrow(rp)
  
  ASEiter<-array(A, dim=c(J,settings$Adim,1))    
  BSEiter<-matrix(B, nrow=J, ncol=1)
  DSEiter<-array(D, dim=c(J,ncat-1,1))
  LLSEiter<-vector()
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
  if (settings$parallel) {
    clusterExport(cl,c("Y","Q","n1cat","N","J","MY","MN","R","missList","indL","indU"),envir=environment())
    clusterEvalQ(cl,c("WrapX","WrapT","WrapZ"))
  }
  ATA 		<- t(A)%*%A #*4
  BTB_INV	<- solve(IQ + ATA)
  d<-D; b<-B
  atemp<-list()
  atemp$Avec<-Avec0<-rep(0,Q)
  atemp$Atemp<-A
  prevA=mat.or.vec(J,Q)
  pRotT<-NA
  pRotV<-NA
  for (i in 1:IT) {
    if (i%%10==1) cat(".")
    if (i%%100==1) cat(":")
    if (i%%500==1) cat("\n",i,"\t : ")
    
    if (settings$parallel) {
      X2		<- simplify2array(parSapply(cl,1:J,WrapX,simplify=FALSE,A=A,b=b,d=d,theta=THat), higher=TRUE)	
    } else {
      X2	<- lapply(1:J,function (x) (WrapX1(x,A=A,b=b,d=d,theta=THat,
                                            settings=settings,R=R,MN=MN,
                                            missList=missList,MY=MY)))
      X2  <- simplify2array(X2,higher=TRUE)	
    }
    X3		<- t(apply(X2,c(1,3),mean))
    b <- -(t(t(rowMeans(X3))))
    d <- -(t(apply(X3, 1, scale, scale=FALSE)))
    if (!settings$dbltrunc) {
      Z		<- colMeans(X2)
    } else {
      if (settings$parallel) {
        X1 		<- parSapply(cl,1:J,WrapZ,simplify=FALSE,A=A,b=b,d=d,theta=THat)
        Z		<- simplify2array(X1, higher=TRUE)	
      } else {
        Z 		<- simplify2array(sapply(1:J,function (x) (WrapZ1(x,A=A,b=b,d=d,theta=THat,
                                                              settings=settings,indL=indL,
                                                              indU=indU))),higher=TRUE)
      } 
    }
    covZ	<- cov(Z)
    prevA <- A 
    Avec0 <- atemp$Avec[1:Q]
    #}
    if (settings$drawA=="lowertriangular") {
      #compute covariances
      covTMC	<- cov(THat)
      covTZMC	<- cov(THat,Z)
      covT		<- covTMC
      covTZ		<- covTZMC
      atemp	<- DrawALowerDiag(covT=covT,covTZ=covTZ,Q,J,N)
    } else if (settings$drawA=="eigen") {
      atemp	<- tryCatch({
        DrawAEigen(covZ = covZ-diag(J),Q)
      }, error = function(e) {
        atemp
      })
    } else {
      atemp	<- DrawA(covZ-diag(J)/n1cat,Q,a=prevA)
    }
    A <- Anew <- as.matrix(atemp$Atemp)
    ATA 		<- t(A)%*%A
    BTB_INV	<- solve(IQ + ATA)
    
    if (settings$parallel) {
      if (Q==1) {
        THat		<- thetanew	<- as.matrix(parSapply(cl,1:N,WrapT,A=A,Z = Z,BTB_INV=BTB_INV,
                                                 b=b,dbltrunc=settings$dbltrunc))
      } else {
        THat		<- thetanew	<- t(parSapply(cl,1:N,WrapTmv,A=A,Z = Z,BTB_INV=BTB_INV,
                                         b=b,dbltrunc=settings$dbltrunc))
      }
    } else {
      if (Q==1) {
        THat	<- thetanew	<- as.matrix(sapply(1:N,function (x) (WrapT(x,A=A,Z = Z,BTB_INV=BTB_INV,
                                                                      b=b,dbltrunc=settings$dbltrunc))))
      } else {
        THat	<- thetanew	<- t(sapply(1:N,function (x) (WrapTmv(x,A=A,Z = Z,BTB_INV=BTB_INV,
                                                                b=b,dbltrunc=settings$dbltrunc))))
      }
    }
    if (TargetA) {
      Target <- FitDATA$XI[,1:settings$Adim]
      Target[Target!=0] <- NA
      Target <- matrix(Target,J,settings$Adim)
      RotT <- targetT(Anew,Tmat=diag(settings$Adim),Target=Target,normalize=TRUE,eps=1e-4, maxit=1000)
      if (sum(colSums(RotT$loadings)<0)>0) {
        for (q in which(colSums(RotT$loadings)<0)) {
          RotT$loadings[,q]<-(-1)*RotT$loadings[,q]
        }
      }
      ASEiter<-abind(ASEiter,as.matrix(RotT$loadings),along=3)    
    } else {
      ASEiter<-abind(ASEiter,as.matrix(A),along=3)    
    }
    BSEiter<-cbind(BSEiter,b)
    DSEiter<-abind(DSEiter,as.matrix(d),along=3)
  }
  if (prl) {
    # cl <- parallel::getDefaultCluster()
    parallel::stopCluster(cl)
  }
  
  return(list(Aiter=ASEiter,Biter=BSEiter,Diter=DSEiter))#,
}
