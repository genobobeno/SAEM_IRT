FixedParamTheta <-
function(FitDATA,rp,IT=100) {
  if (settings$Adim>1) {
    if (!is.na(FitDATA$AR)) {
      A=as.matrix(FitDATA$AR$loadings)
      THat=as.matrix(FitDATA$Trot[,1:settings$Adim])
    } else {
      A=as.matrix(FitDATA$A)
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
  if (settings$fm %in% c("new","pca","eigen")) settings$fm <- "camilli"
  prior<-list(tmu=settings$tmu,tsigma=settings$tsigma)
  J<-ncol(rp)
  N<-nrow(rp)
  ifelse(settings$guess,W<-DrawW(aa=A,bb=B,cc=C,tt=THat,rp=rp),W<-NA)  
  Z<-SampZ(aa=A,bb=B,that=THat,rp=rp,w=W)    
  AFSEiter<-array(A, dim=c(J,settings$Adim,1))    
  BFSEiter<-matrix(B, nrow=J, ncol=1)
  ifelse(settings$guess,CFSEiter<-matrix(C, nrow=J, ncol=1),CFSEiter<-NA)
  FTiter<-array(THat, dim=c(N,settings$Adim,1))
  #print(paste("Running Fixed Parameter Theta Estimates using",IT,"iterations"))  
  for (it in 1:IT) {
    #if(it%%10==0) print(paste(it,"th Iteration",sep=""))
    Z<-SampZ(aa=A,bb=B,that=THat,rp=rp,w=W)    
    LL<-GIFAFullLL(A,B,Z,THat,prior=prior)
    FTiter<-abind(FTiter,as.matrix(THat),along=3)
    if (settings$guess) {
      PSI<-GIFAEstimate(aa=A,bb=B,zz=Z,tt=THat,settings=settings,w=W,rp=rp)    
    } else {
      PSI<-GIFAEstimate(aa=A,bb=B,zz=Z,tt=THat,settings=settings)
    }
    THat<-SampT(aa=PSI$A,bb=PSI$B,zz=Z,rp=rp,prior=prior)  
    AFSEiter<-abind(AFSEiter,as.matrix(PSI$A),along=3)    
    BFSEiter<-cbind(BFSEiter,PSI$B)
    if (settings$guess) CFSEiter<-cbind(CFSEiter,PSI$C)
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
