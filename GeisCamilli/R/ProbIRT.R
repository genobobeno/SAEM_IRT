ProbIRT <-
function(xi,theta,j=NA,Adim=settings$Adim,guess=settings$guess) { # gives back a N X J frame of probabilities... or N x 1
  if (Adim>1) print("NOTE: Multidimensional 2PL or 3PL is not implemented")
  if (Adim==1) {
    if (guess) {
      aa<-as.matrix(xi[,1])
      bb<-as.matrix(xi[,2])
      cc<-as.matrix(xi[,3])
    } else {
      aa<-as.matrix(xi[,1])
      bb<-as.matrix(xi[,2])
    }
    if (is.na(j)) {
      Z<-outer(theta,as.vector(bb),'-')
      AT<-t(aa%*%t(as.matrix(rep(1,length(theta)))))
      if (guess) C<-t(cc%*%t(as.matrix(rep(1,nrow(as.matrix(theta))))))
    } else {
      AT<-aa[j]
      Z<-(theta-bb[j])
      if (guess) C<-cc[j]
    }
  }
  ifelse(guess,p<-C+(1-C)*exp(1.7*AT*Z)/(1+exp(1.7*AT*Z)),p<-exp(1.7*AT*Z)/(1+exp(1.7*AT*Z)))
  return (p)
}
