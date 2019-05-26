InitializePrior <-
function(rp,settings=settings) {
  J<-ncol(rp)
  N<-nrow(rp)
  if (sum(rp==9)>0) {
    for (i in 1:J) {
      TF<-rp[,i]==9
      if (sum(TF)>0) rp[TF,i]<-mean(rp[!TF,i])
    }
  }
  POLY<-FALSE
  if (!is.na(settings$ncat) & settings$ncat>2) POLY<-TRUE
  if (settings$initialize=="random") {
    if (settings$Adim==1) {
      THat<-rnorm(N,0,1)
      A<-runif(J,0.5,1.5)
    } else {
      A<-mat.or.vec(J,settings$Adim)
      THat<-mat.or.vec(N,settings$Adim)
      for (i in 1:settings$Adim) {
        A[,i]<-runif(J,0.5,1.5)
        THat[,i]<-rnorm(N,0,1) # Small squeeze to prevent infinities
      }
    }
    B<-rnorm(J,0,1)  
    XI<-cbind(A,B)
    if (settings$guess) {
      C<-runif(J,0.01,0.3)
      XI<-cbind(XI,C)
    }
  } else {  #if (settings$initialize=="best") {
    if (settings$Adim==1) {
      A<-rep(1,J)
      if (!POLY) {
        THat<-qnorm(0.005+0.99*(rowSums(rp)/J)) # Small squeeze to prevent infinities
      } else {
        THat<-qnorm(0.005+0.99*((rowSums(rp)/(settings$ncat-1))/J)) # Small squeeze to prevent infinities
      }
    } else {
      A0 = matrix(c(rep(0.5,J),rep(0.1,J*(settings$Adim-1))),J,settings$Adim)+0.1*matrix(runif(J*settings$Adim),J,settings$Adim)
      A<-mat.or.vec(J,settings$Adim)+A0
      if (POLY) {
        THat<-as.matrix(qnorm(0.005+0.99*((rowSums(rp)/(settings$ncat-1))/J))) # Small squeeze to prevent infinities
      } else {
        THat<-as.matrix(qnorm(0.005+0.99*(rowSums(rp)/J))) # Small squeeze to prevent infinities
      }
      for (i in 2:settings$Adim) {
        THat<-cbind(THat,THat[,1]+rnorm(N,0,0.01)) #adding a little stochasticity...
      }
    }
    if (!is.na(settings$ncat) & settings$ncat>2) {
      B<-qnorm(1-colSums(rp)/(settings$ncat-1)/N)*2
    } else {
      B<-qnorm(1-colSums(rp)/N)*2
    }
    XI<-cbind(A,B)
    if (settings$guess) {
      if (settings$Adim==1) {
        rpGuess0 = rp[THat<quantile(THat,probs=0.12),]; rpGuess1 = rp[THat>quantile(THat,probs=0.88),]
      } else {
        rpGuess0 = rp[apply(THat,1,sum)/settings$Adim<quantile(apply(THat,1,sum)/settings$Adim,probs=0.12),]
        rpGuess1 = rp[apply(THat,1,sum)/settings$Adim>quantile(apply(THat,1,sum)/settings$Adim,probs=0.88),]
      }
      Gap = apply(rpGuess1,2,mean,na.rm=TRUE) - apply(rpGuess0,2,mean,na.rm=TRUE)
      C<-apply(rpGuess0,2,mean,na.rm=TRUE)
      C[Gap>0.5 & apply(rpGuess1,2,mean,na.rm=TRUE)>0.78]<-(apply(rpGuess0,2,mean,na.rm=TRUE)*(1-Gap))[Gap>0.5 & apply(rpGuess1,2,mean,na.rm=TRUE)>0.78]
      #C<-seq(0.35,0.01,length.out = J)[order(Gap)] #C<-runif(J,0.01,0.3)
      XI<-cbind(XI,C)
    }
  }
  if (POLY) {
    D = matrix(seq(-1.0,1.0,length.out = (settings$ncat-1)),J,(settings$ncat-1), byrow=TRUE)
  } else {D<-NA}
  # theta = matrix(rnorm(n*Q,0,1),n,Q)
  # A		= matrix(runif(m*Q)-.5,m,Q)
  # b		= matrix(rnorm(m),m,1);	
  # d		= matrix(seq(-1.0,1.0,length.out = n1cat),m,n1cat, byrow=TRUE)
  list(XI=XI,THat=THat,D=D)
}
