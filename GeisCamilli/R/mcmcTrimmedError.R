mcmcTrimmedError<-function(x,end=800,start=600,trimA=8,trimB=5,trimC=5,settings,tau=FALSE) {
  if (!tau) {
    Atrim = start+trimA*1:(floor((end-start)/trimA))
    for (iii in 1:dim(x$Aiter)[2]) {
      if (iii==1) {
        AE<-apply(x$Aiter[,1,Atrim],1,sd)
      } else {
        AE<-cbind(AE,apply(x$Aiter[,iii,Atrim],1,sd))
      }
    }
    Btrim = start+trimB*1:(floor((end-start)/trimB))
    BE<-apply(x$Biter[,Btrim],1,sd)
    if (settings$guess) {
      Ctrim = start+trimC*1:(floor((end-start)/trimC))
      CE<-apply(x$Citer[,Ctrim],1,sd)
      cbind(AE,BE,CE)
    } else {
      cbind(AE,BE)
    }
  } else {
    Dtrim = start+trimB*1:(floor((end-start)/trimB))
    for (iii in 1:dim(x$Diter)[2]) {
      if (iii==1) {
        DE<-apply(x$Diter[,1,Dtrim],1,sd)
      } else {
        DE<-cbind(DE,apply(x$Diter[,iii,Dtrim],1,sd))
      }
    }
    return(DE)
  }
}