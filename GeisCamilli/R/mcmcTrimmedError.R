mcmcTrimmedError<-function(x,end=800,start=600,trimA=8,trimB=5,trimC=5,settings) {
  for (iii in 1:dim(x$Aiter)[2]) {
    Atrim = start+trimA*1:(floor((end-start)/trimA))
    if (iii==1) {
      AE<-apply(x$Aiter[,1,Atrim],1,sd)
    } else {
      AE<-cbind(AE,apply(x$Aiter[,1,Atrim],1,sd))
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
}