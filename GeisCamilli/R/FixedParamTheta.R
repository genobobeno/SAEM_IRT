FixedParamTheta <-
  function(FitDATA,rp,IT=settings$nesttheta) {
    cat("\n Running Fixed Parameter Theta Estimates \n")
    cat(paste(IT,"Draws \n"))
    atemp	<- list()
    if (settings$Adim>1) {
      if (!is.na(FitDATA$AR)) {
        atemp$Atemp<-A<-as.matrix(FitDATA$AR$loadings)
        THat=as.matrix(FitDATA$Trot[,1:settings$Adim])
      } else {
        atemp$Atemp<-as.matrix(FitDATA$A)
        THat=as.matrix(FitDATA$That[,1:settings$Adim])
      }
      B=FitDATA$B
      C=FitDATA$C
    } else {
      A=FitDATA$A
      B=FitDATA$B
      C=FitDATA$C
      THat=FitDATA$That[,1]
    }
    Alines<-A
    Blines<-B
    settings=FitDATA$settings
    if (!is.na(settings$ncat) & settings$ncat>2) tau<-FitDATA$tau
    if (settings$fm %in% c("new","pca","eigen")) settings$fm <- "camilli"
    prior<-list(tmu=settings$tmu,tsigma=settings$tsigma)
    J<-ncol(rp)
    N<-nrow(rp)
    AFSEiter<-array(A, dim=c(J,settings$Adim,1))    
    BFSEiter<-matrix(B, nrow=J, ncol=1)
    FTiter<-array(THat, dim=c(N,settings$Adim,1))
    ifelse(settings$guess,CFSEiter<-matrix(C, nrow=J, ncol=1),CFSEiter<-NA)
    #  ifelse(settings$guess,W<-DrawW(aa=A,bb=B,cc=C,tt=THat,rp=rp),W<-NA)  
    # Z<-SampZFast(aa=A,bb=B,that=THat,srp=rp,w=W)    
    #print(paste("Running Fixed Parameter Theta Estimates using",IT,"iterations"))  
    if (is.na(settings$ncat) | settings$ncat==2) {
      for (it in 1:IT) {
        ifelse(settings$guess,W<-DrawW(aa=A,bb=B,cc=C,tt=THat,rp=rp),W<-NA)
        #if(it%%10==0) print(paste(it,"th Iteration",sep=""))
        Z<-SampZFast(aa=A,bb=B,that=THat,srp=rp,w=W)    
        LL<-GIFAFullLL(A,B,Z,THat,prior=prior)
        if (settings$guess) {
          PSI<-GIFAEstimate(aa=A,bb=B,zz=Z,tt=THat,settings=settings,w=W,rp=rp)    
        } else {
          PSI<-GIFAEstimate(aa=A,bb=B,zz=Z,tt=THat,settings=settings)
        }
        THat<-SampT(aa=PSI$A,bb=PSI$B,zz=Z,rp=rp,prior=prior)  
        AFSEiter<-abind(AFSEiter,as.matrix(PSI$A),along=3)    
        BFSEiter<-cbind(BFSEiter,PSI$B)
        FTiter<-abind(FTiter,as.matrix(THat),along=3)
        if (settings$guess) CFSEiter<-cbind(CFSEiter,PSI$C)
      }
    } else {
      ATA 		<- t(A)%*%A #*4
      BTB_INV	<- solve(diag(settings$Adim) + ATA)
      d = FitDATA$tau-matrix(rep(B,ncol(FitDATA$tau)),
                             length(B),ncol(FitDATA$tau))
      for (it in 1:IT) {
        X2		<- simplify2array(parSapply(cl,1:length(B),WrapX,A=A,b=B,
                                        d=d,theta=THat,
                                        simplify=FALSE), higher=TRUE)	
        X3		<- t(apply(X2,c(1,3),mean))
        b         <- -t(t(rowMeans(X3)));    
        dd        <- -t(apply(X3, 1, scale, scale=FALSE))   
        zHat	<- apply(X2,c(2,3),mean) 
        # if (settings$drawA=="lowertriangular") {
        #   #compue covariances
        #   covTMC	<- cov(THat)
        #   covTZMC	<- cov(THat,Z)
        #   covT		<- covT  + alpha[i]*(covTMC-covT)
        #   covTZ		<- covTZ + alpha[i]*(covTZMC-covTZ)
        #   atemp	<- DrawALowerDiag(covT=covT,covTZ=covTZ,Q,J,N)
        # } else 
        if (settings$drawA=="eigen") {
          atemp	<- tryCatch({
            DrawAEigen(covZ = cov(zHat)-diag(J),settings$Adim)
          }, error = {
            atemp$Atemp<-atemp$Atemp+matrix(rnorm(length(atemp$Atemp),0,sd = 0.05),
                                            nrow = nrow(atemp$Atemp),
                                            ncol = ncol(atemp$Atemp))
            atemp
          })
        } else {
          atemp	<- DrawA(cov(zHat)-diag(J)/n1cat,settings$Adim,a=A)
        }
        a <- as.matrix(atemp$Atemp)
        if (settings$Adim==1) {
          THat	<- as.matrix(parSapply(cl,1:nrow(zHat),WrapT,A=A,Z = zHat,BTB_INV=BTB_INV,b=B,dbltrunc=settings$dbltrunc))
        } else {
          THat	<- t(parSapply(cl,1:nrow(zHat),WrapTmv,A=A,Z = zHat,BTB_INV=BTB_INV,b=B,dbltrunc=settings$dbltrunc))
        }
        FTiter<-abind(FTiter,as.matrix(THat),along=3)
        #if(it%%10==0) print(paste(it,"th Iteration",sep=""))
        AFSEiter<-abind(AFSEiter,a,along=3)    
        BFSEiter<-cbind(BFSEiter,b)
      }
    }
    if (settings$Adim>1) {
      FTErr<-as.matrix(apply(FTiter,c(1,2),sd))
      FAErr<-as.matrix(apply(AFSEiter,c(1,2),sd))
      FT<-as.matrix(apply(FTiter,c(1,2),mean))
    } else {
      FTErr<-as.vector(apply(as.matrix(FTiter[,1,]),1,sd))
      FAErr<-as.matrix(apply(AFSEiter[,1,],1,sd))
      FT<-as.vector(apply(as.matrix(FTiter[,1,]),1,mean))
    } 
    FBErr<-as.vector(apply(BFSEiter,1,sd))
    if (settings$plots) {
      if (settings$Adim>1) {
        par(mfrow=c(1,settings$Adim),mar=c(3,3,1,1))
        for (i in 1:settings$Adim) {
          plot(FTiter[1,i,],type="n",xlim=c(0,length(FTiter[1,1,])),ylim=range(FTiter),main=paste("Theta",i,"of 10 Examinees"))
          for (j in 1:10) {
            lines(FTiter[j,i,],col=j)
          }
          abline(h=FT[1:10,i])
        }
        par(mfrow=c(1,settings$Adim+1),mar=c(3,3,1,1))
        for (i in 1:settings$Adim) {
          plot(AFSEiter[1,i,],type="n",xlim=c(0,length(AFSEiter[1,1,])),ylim=range(AFSEiter),main=paste("A",i))
          for (j in 1:J) {
            lines(AFSEiter[j,i,],col=j)
          }
          abline(h=Alines[,i])
        }
        plot(BFSEiter[1,],type="n",xlim=c(0,length(BFSEiter[1,])),ylim=range(BFSEiter),main="B")
        for (j in 1:J) {
          lines(BFSEiter[j,],col=j)
        }
        abline(h=Blines)
      } else {
        par(mfrow=c(1,1),mar=c(3,3,1,1))
        plot(FTiter[1,1,],type="n",xlim=c(0,length(FTiter[1,1,])),ylim=range(FTiter),main=paste("Simulating Theta of 10 Examinees"))
        for (j in 1:10) {
          lines(FTiter[j,1,],col=j)
        }
        abline(h=FT[1:10])
        par(mfrow=c(1,2),mar=c(3,3,1,1))
        plot(AFSEiter[1,1,],type="n",xlim=c(0,length(AFSEiter[1,1,])),ylim=range(AFSEiter),main=paste("A",1))
        for (j in 1:J) {
          lines(AFSEiter[j,1,],col=j)
        }
        abline(h=Alines)
        plot(BFSEiter[1,],type="n",xlim=c(0,length(BFSEiter[1,])),ylim=range(BFSEiter),main="B")
        for (j in 1:J) {
          lines(BFSEiter[j,],col=j)
        }
        abline(h=Blines)
      }
    }
    return(list(FT=FT,FTErr=FTErr,FAErr=FAErr,FBErr=FBErr))
  }
