InitializeParams <-
function(rp,settings=settings) {
  J<-ncol(rp)
  N<-nrow(rp)
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
  } else if (settings$initialize=="best") {
    if (settings$Adim==1) {
      A<-rep(1,J)
      THat<-qnorm(0.005+0.99*(rowSums(rp)/J)) # Small squeeze to prevent infinities
    } else {
      if (settings$Adim==1) {A0 = 1} else {A0 = matrix(c(rep(0.85,J),rep(0.45,J*(settings$Adim-1))),J,settings$Adim)+0.3*matrix(runif(J*settings$Adim),J,settings$Adim)}
      A<-mat.or.vec(J,settings$Adim)+A0
      THat<-as.matrix(qnorm(0.005+0.99*(rowSums(rp)/J))) # Small squeeze to prevent infinities
      for (i in 2:settings$Adim) {
        THat<-cbind(THat,THat[,1])
      }
    }
    B<-qnorm(1-colSums(rp)/N)*2
    XI<-cbind(A,B)
    if (settings$guess) {
      C<-runif(J,0.01,0.3)
      XI<-cbind(XI,C)
    }
  }
  return(list(XI=XI,THat=THat))
}
