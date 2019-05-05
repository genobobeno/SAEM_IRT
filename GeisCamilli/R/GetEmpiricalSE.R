GetEmpiricalSE <-
function(FitDATA,rp,IT=settings$EmpIT,estgain=settings$estgain,thinA=settings$thinA,thinB=settings$thinB,
         recordChain=FALSE) {
  cat("\n")
  print(paste("Starting",IT,"iterations of Empirical SEs, Thinning A:",thinA,"; Thinning B:",thinB))
  settings=FitDATA$settings
  settings$fm<-"camilli"
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
    if (settings$fm %in% c("new","pca")) settings$fm <- "camilli"
  } else {
    A=FitDATA$A
    B=FitDATA$B
    C=FitDATA$C
    THat=FitDATA$That[,1]
  }
  Alines<-A
  Blines<-B
  prior<-list(tmu=settings$tmu,tsigma=settings$tsigma)
  J<-ncol(rp)
  N<-nrow(rp)

  iDnew = rep(0,(settings$Adim+1)*J) # Jacobian Delta
  iGnew = mat.or.vec((settings$Adim+1)*J,(settings$Adim+1)*J) # Hessian Covariance G
  iJnew = mat.or.vec((settings$Adim+1)*J,(settings$Adim+1)*J) # Hessian D
  iDeltak0<-rep(0,(settings$Adim+1)*J) # Jacobian Delta
  iDk0<-mat.or.vec((settings$Adim+1)*J,(settings$Adim+1)*J) # Hessian D
  iGk0<-mat.or.vec((settings$Adim+1)*J,(settings$Adim+1)*J) # Hessian Covariance G

  oDnew = rep(0,(settings$Adim+1)*J) # Jacobian Delta
  oGnew = mat.or.vec((settings$Adim+1)*J,(settings$Adim+1)*J) # Hessian Covariance G
  oJnew = mat.or.vec((settings$Adim+1)*J,(settings$Adim+1)*J) # Hessian D
  oDeltak0<-rep(0,(settings$Adim+1)*J) # Jacobian Delta
  oDk0<-mat.or.vec((settings$Adim+1)*J,(settings$Adim+1)*J) # Hessian D
  oGk0<-mat.or.vec((settings$Adim+1)*J,(settings$Adim+1)*J) # Hessian Covariance G
  XIMean0<-rep(1,J*(1+settings$Adim))  
  
  ifelse(settings$guess,W<-DrawW(aa=A,bb=B,cc=C,tt=THat,rp=rp),W<-NA)  
  #Z<-SampZ(aa=A,bb=B,that=THat,rp=rp,w=W)    
  Z<-SampZFast(aa=A,bb=B,that=THat,srp=rp,w=W)    
  ASEiter<-array(A, dim=c(J,settings$Adim,1))    
  BSEiter<-matrix(B, nrow=J, ncol=1)
  ifelse(settings$guess,CSEiter<-matrix(C, nrow=J, ncol=1),CSEiter<-NA)
  ZSEiter<-matrix(Z, nrow=N, ncol=1)
  TSEiter<-array(THat, dim=c(N,settings$Adim,1))
  LLSEiter<-vector()
  #print(paste("Running Standard Errors using",IT,"iterations"))
  
  for (it in 1:IT) {
    if (it%%10==1) cat(".")
    if (it%%100==1) cat(":")
    if (it%%500==1) cat(" \n",it,"\t")
    if (settings$guess) {
      W<-DrawW(aa=A,bb=B,cc=C,tt=THat,rp=rp)
    } else {
      W<-NA  
    }
    #Z<-SampZ(aa=A,bb=B,that=THat,rp=rp,w=W)    
    Z<-SampZFast(aa=A,bb=B,that=THat,srp=rp,w=W)    
    LL<-GIFAFullLL(A,B,Z,THat,prior=prior)
    if (settings$guess) {
      PSI<-GIFAEstimate(aa=A,bb=B,zz=Z,tt=THat,settings=settings,w=W,rp=rp,EmpT=TRUE)    
    } else {
      PSI<-GIFAEstimate(aa=A,bb=B,zz=Z,tt=THat,settings=settings,EmpT=TRUE)
    }
    A<-PSI$A
    B<-PSI$B
    C<-PSI$C    
    ASEiter<-abind(ASEiter,as.matrix(A),along=3)    
    BSEiter<-cbind(BSEiter,B)
    if (settings$guess) CSEiter<-cbind(CSEiter,C)
    ZSEiter<-cbind(ZSEiter,Z)
    TSEiter<-abind(TSEiter,as.matrix(THat),along=3)
    LLSEiter<-c(LLSEiter,LL)
    gain<-1.0/it^estgain
    #     Hess<-GetHessian(A,B,Z)
    #     Jcb<-GetJacobian(A,B,Z)
    #     Dk0<-Dk0+gain*(Hess-Dk0)
    #     Gk0<-Gk0+gain*(as.matrix(Jcb)%*%t(as.matrix(Jcb))-Gk0)
    #     Deltak0<-Deltak0+gain*(Jcb-Deltak0)
    #     Hk0<-Dk0+Gk0-as.matrix(Deltak0)%*%t(as.matrix(Deltak0))
    #     ifelse(it>(settings$burnin+10),gain<-1.0/(it-settings$burnin),gain<-1) #settings$estgain
    #       Hess<-GetHessian(A,B,Z)
    #       Jcb<-GetJacobian(A,B,Z)
    #       DHess<-sechol(Hk0)
    #       Hk0<-t(DHess)%*%DHess
    #JH    <- GetJacobHess(A,B,Z)
    #JH<-GetErrorLogitApp(A=A,B=B,C=C,TH=THat,RP=rp)
    oJH<-GetErrorOgive(A=A,B=B,C=C,TH=THat,Z=Z,RP=rp)
    oJacob <- oJH$Jacob
    oHess  <- oJH$Hess
    oDk0<-oDk0+gain*(oHess-oDk0)
    oGk0<-oGk0+gain*(as.matrix(oJacob)%*%t(as.matrix(oJacob))-oGk0)
    oDeltak0<-oDeltak0+gain*(oJacob-oDeltak0)
    iJH<-GetErrorLogitApp(A=A,B=B,C=C,TH=THat,RP=rp)
    iJacob <- iJH$Jacob
    iHess  <- iJH$Hess
    iDk0<-iDk0+gain*(iHess-iDk0)
    iGk0<-iGk0+gain*(as.matrix(iJacob)%*%t(as.matrix(iJacob))-iGk0)
    iDeltak0<-iDeltak0+gain*(iJacob-iDeltak0)
    #     J2    <- matrix(Jacob,J*(1+settings$Adim),1)%*% matrix(Jacob,1,J*(1+settings$Adim))
#     Dnew = Dnew + Jacob
#     Gnew = Gnew + Hess
#     Jnew = Jnew + J2  
#     #   print(c("check out Jacob ",mean(abs(Jacob)),mean(Jacob)))
#     #   print(c("        ",max(Jacob),which(Jacob==max(Jacob),arr.ind=TRUE) ))
#     Dold = Dnew*gain
#     Gold = Gnew*gain
#     Jold = Jnew*gain
#     H1 = Gold
#     H2 = Jold - matrix(Dold,J*(settings$Adim+1),1)%*%matrix(Dold,1,J*(settings$Adim+1))
#     Hk0 = -(H1+H2)
    THat<-SampT(aa=A,bb=B,zz=Z,rp=rp,prior=prior)  
  }
  MCthin<-MCthinA<-1:floor(IT/thinA)*thinA
  MCthinB<-1:floor(IT/thinB)*thinB
  if (settings$Adim>1) {
    SEA<-as.matrix(apply(ASEiter[,,MCthinA],c(1,2),sd))
    MEA<-as.matrix(apply(ASEiter[,,MCthinA],c(1,2),mean))
    MCSA<-sqrt(as.matrix(apply(ASEiter[,,MCthinA],c(1,2),function(x) (initseq(x)$var.pos))))
  } else {
    SEA<-as.vector(apply(as.matrix(ASEiter[,1,MCthinA]),1,sd))
    MEA<-as.vector(apply(as.matrix(ASEiter[,1,MCthinA]),1,mean))
    MCSA<-sqrt(as.vector(apply(ASEiter[,1,MCthinA],1,function(x) (initseq(x)$var.pos))))
  } 
  SEB<-as.vector(apply(BSEiter[,MCthinB],1,sd))
  MEB<-as.vector(apply(BSEiter[,MCthinB],1,mean))
  MCSB<-sqrt(as.vector(apply(BSEiter[,MCthinB],1,function(x) (initseq(x)$var.pos))))
  ifelse(settings$guess,SEC<-as.vector(apply(CSEiter[,MCthin],1,sd)),SEC<-NA)
  ifelse(settings$guess,MEC<-as.vector(apply(CSEiter[,MCthin],1,mean)),MEC<-NA)
  ifelse(settings$guess,MCSC<-sqrt(as.vector(apply(CSEiter[,MCthin],1,function(x) (initseq(x)$var.pos)))),MCSC<-NA)
  if (settings$Adim>1) {
    SET<-as.matrix(apply(TSEiter[,,MCthin],c(1,2),sd))
    MET<-as.matrix(apply(TSEiter[,,MCthin],c(1,2),mean))
    MCST<-sqrt(as.matrix(apply(TSEiter[,,MCthin],c(1,2),function(x) (initseq(x)$var.pos))))
  } else {
    SET<-as.vector(apply(as.matrix(TSEiter[,1,MCthin]),1,sd))
    MET<-as.vector(apply(as.matrix(TSEiter[,1,MCthin]),1,mean))
    MCST<-sqrt(as.vector(apply(TSEiter[,1,MCthin],1,function(x) (initseq(x)$var.pos))))
  } 
  if (settings$plots) {
    par(mfrow=c(1,settings$Adim+1),mar=c(3,3,1,1))
    for (i in 1:settings$Adim) {
      plot(ASEiter[1,i,],type="n",xlim=c(0,length(ASEiter[1,1,])),ylim=range(ASEiter),main=paste("A",i))
      for (j in 1:J) {
        lines(ASEiter[j,i,],col=j)
      }
      abline(h=Alines[,i])
    }
    plot(BSEiter[1,],type="n",xlim=c(0,length(BSEiter[1,])),ylim=range(BSEiter),main="B")
    for (j in 1:J) {
      lines(BSEiter[j,],col=j)
    }
    abline(h=Blines)
  }
  oHk0<-oDk0+oGk0-as.matrix(oDeltak0)%*%t(as.matrix(oDeltak0))
  iHk0<-iDk0+iGk0-as.matrix(iDeltak0)%*%t(as.matrix(iDeltak0))
  
  VARLMIo<-diag(ginv((-1)*oHk0))
  VARLMIi<-diag(ginv((-1)*iHk0))
  VARLMIo<-matrix(VARLMIo,J,settings$Adim+1,byrow=T)[,c(2:(settings$Adim+1),1)]
  VARLMIi<-matrix(VARLMIi,J,settings$Adim+1,byrow=T)[,c(2:(settings$Adim+1),1)]
  if (!recordChain){
    return(list(SEA=SEA,MEA=MEA,MCSA=MCSA,SEB=SEB,MEB=MEB,MCSB=MCSB,SEC=SEC,MEC=MEC,MCSC=MCSC,
                SET=SET,MET=MET,MCST=MCST,VARLMIi=VARLMIi,VARLMIo=VARLMIo))
  } else {
    return(list(SEA=SEA,MEA=MEA,MCSA=MCSA,SEB=SEB,MEB=MEB,MCSB=MCSB,SEC=SEC,MEC=MEC,MCSC=MCSC,
                SET=SET,MET=MET,MCST=MCST,VARLMIi=VARLMIi,VARLMIo=VARLMIo,Aiter=ASEiter,    
                Biter=BSEiter,Citer=CSEiter,Liter=LLSEiter))
  }
}
