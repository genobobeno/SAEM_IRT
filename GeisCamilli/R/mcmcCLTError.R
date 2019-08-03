mcmcCLTError<-function(x,end=800,start=400,settings,item.ind=TRUE,tau=FALSE) {
  #x=MCMCDATA;end = burnin[2];start = burnin[1]
  library(mcmcse)
  # Atrim = start+trimA*1:(floor((burnin-start)/trimA))
  Trim = start:end
  if (!tau) {
    if (!item.ind) {  ####  IF we should run the covariance matrix of EVERYTHING (NOT assume item independence), do this.
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
      
      if ("ncat" %in% names(settings) && !is.na(settings$ncat) & settings$ncat>2) {
        for (iii in 1:dim(x$Diter)[2]) {
          MCCOV<-cbind(MCCOV,t(matrix(x$Diter[,iii,Trim],nrow=dim(x$Diter)[1],ncol=length(Trim))))
        }
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
    } else {
      AE<-matrix(0,nrow=dim(x$Aiter)[1],ncol=dim(x$Aiter)[2])
      BE<-matrix(0,nrow=dim(x$Aiter)[1],ncol=1)
      CE<-matrix(0,nrow=dim(x$Aiter)[1],ncol=1)
      for (j in 1:dim(x$Aiter)[1]) {
        MCCOV<-matrix(x$Aiter[j,1,Trim],ncol=1,nrow=length(Trim))
        if (settings$Adim>1) {
          for (iii in 2:dim(x$Aiter)[2]) {
            MCCOV<-cbind(MCCOV,matrix(x$Aiter[j,iii,Trim],ncol=1,nrow=length(Trim)))
          }
        }
        MCCOV<-cbind(MCCOV,matrix(x$Biter[j,Trim],ncol=1,nrow=length(Trim)))
        if (settings$guess) {
          MCCOV<-cbind(MCCOV,matrix(x$Citer[,Trim],ncol=1,nrow=length(Trim)))
        }
        SDE<-sqrt(diag(mcse.initseq(MCCOV)$cov))
        AE[j,]<- SDE[1:(dim(x$Aiter)[2])]
        BE[j] <- SDE[dim(x$Aiter)[2]+1]
        if (settings$guess) {
          CE[j] <- SDE[dim(x$Aiter)[2]+2]
        }
      }
      if (settings$guess) {
        cbind(AE,BE,CE)
      } else {
        cbind(AE,BE)
      }
    }
  } else {
    if (!item.ind) {  ####  IF we should run the covariance matrix of EVERYTHING (NOT assume item independence), do this.
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
      for (iii in 1:dim(x$Diter)[2]) {
        MCCOV<-cbind(MCCOV,t(matrix(x$Diter[,iii,Trim],nrow=dim(x$Diter)[1],ncol=length(Trim))))
      }
      
      if (settings$ncat>2) {
        for (iii in 1:dim(x$Diter)[2]) {
          MCCOV<-cbind(MCCOV,t(matrix(x$Diter[,iii,Trim],nrow=dim(x$Diter)[1],ncol=length(Trim))))
        }
      }
      #Btrim = start+trimB*1:(floor((burnin-start)/trimB))
      DE<-matrix(sqrt(diag(mcse.initseq(MCCOV)$cov)[dim(x$Aiter)[1]*(dim(x$Aiter)[2]+1+settings$guess)+1:(dim(x$Diter)[1]*dim(x$Diter)[2])]),
                 nrow=dim(x$Diter)[1],ncol=dim(x$Diter)[2])
      DE
    } else {
      AE<-matrix(0,nrow=dim(x$Aiter)[1],ncol=dim(x$Aiter)[2])
      DE<-matrix(0,nrow=dim(x$Diter)[1],ncol=dim(x$Diter)[2])
      BE<-matrix(0,nrow=dim(x$Aiter)[1],ncol=1)
      CE<-matrix(0,nrow=dim(x$Aiter)[1],ncol=1)
      for (j in 1:dim(x$Aiter)[1]) { #j=28
        MCCOV<-matrix(x$Aiter[j,1,Trim],ncol=1,nrow=length(Trim))
        if (settings$Adim>1) {
          for (iii in 2:dim(x$Aiter)[2]) {
            MCCOV<-cbind(MCCOV,matrix(x$Aiter[j,iii,Trim],ncol=1,nrow=length(Trim)))
          }
        }
        MCCOV<-cbind(MCCOV,matrix(x$Biter[j,Trim],ncol=1,nrow=length(Trim)))
        if (settings$guess) {
          MCCOV<-cbind(MCCOV,matrix(x$Citer[,Trim],ncol=1,nrow=length(Trim)))
        }
        if (settings$ncat>2) {
          for (iii in 1:dim(x$Diter)[2]) {
            MCCOV<-cbind(MCCOV,matrix(x$Diter[j,iii,Trim],ncol=1,nrow=length(Trim)))
          }
        }
        SDE<-sqrt(diag(mcse.initseq(MCCOV)$cov))
        DE[j,]<-SDE[dim(x$Aiter)[2]+1+settings$guess+1:(dim(x$Diter)[2])]
      }
      DE
    }
  }
}