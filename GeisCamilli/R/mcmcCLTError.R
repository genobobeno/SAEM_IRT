mcmcCLTError<-function(x,end=800,start=400,settings) {
  #x=MCMCDATA;end = burnin[2];start = burnin[1]
  library(mcmcse)
  # Atrim = start+trimA*1:(floor((burnin-start)/trimA))
  Trim = start:end
  MCCOV<-t(matrix(x$Aiter[,1,Trim],nrow=dim(x$Aiter)[1],ncol=length(Trim)))
  if (settings$Adim>1) {
    for (iii in 2:dim(x$Aiter)[2]) {
      MCCOV<-cbind(MCCOV,t(matrix(x$Aiter[,iii,Trim],nrow=dim(x$Aiter)[1],ncol=length(Trim))))
    }
  }
  MCCOV<-cbind(MCCOV,t(matrix(x$Biter[,Trim],nrow=dim(x$Biter)[1],ncol=length(Trim))))
  
  if (settings$guess) {
    MCCOV<-cbind(MCCOV,t(matrix(x$Citer[,Trim],nrow=dim(x$Biter)[1],ncol=length(Trim))))
  }
  
  #Btrim = start+trimB*1:(floor((burnin-start)/trimB))
  AE<-matrix(sqrt(diag(mcse.initseq(MCCOV)$cov)[1:(dim(x$Aiter)[1]*dim(x$Aiter)[2])]),
             nrow=dim(x$Aiter)[1],ncol=dim(x$Aiter)[2])
  BE<-matrix(sqrt(diag(mcse.initseq(MCCOV)$cov)[(dim(x$Aiter)[1]*dim(x$Aiter)[2])+1:dim(x$Biter)[1]]),
             nrow=dim(x$Biter)[1],ncol=1)
  if (settings$guess) {
    CE<-matrix(sqrt(diag(mcse.initseq(MCCOV)$cov)[(dim(x$Aiter)[1]*dim(x$Aiter)[2]+dim(x$Biter)[1])+1:dim(x$Citer)[1]]),
               nrow=dim(x$Citer)[1],ncol=1)
    cbind(AE,BE,CE)
  } else {
    cbind(AE,BE)
  }
}
