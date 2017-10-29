ProbOgive <-
function(xi,theta,j=NA,guess=FALSE,tau=NA) { # gives back a N X J frame of probabilities... or N x 1
  #structure and settings have same list entries for calculations
  #tau will be J x K
  theta = as.matrix(theta)
  if (guess) {
    aa<-as.matrix(xi[,1:(ncol(xi)-2)])
    bb<-as.matrix(xi[,ncol(xi)-1])
    if (is.na(j)) {
      C<-t(as.matrix(xi[,ncol(xi)])%*%t(as.matrix(rep(1,nrow(theta)))))
    } else {
      C<-xi[j,ncol(xi)]
    }
    if (!is.na(tau)[1]) C<-array(rep(C,ncol(tau)),c(nrow(theta),ifelse(is.na(j),nrow(xi),1),ncol(tau)))
  } else {
    aa<-as.matrix(xi[,1:(ncol(xi)-1)])
    bb<-as.matrix(xi[,ncol(xi)])
  }
  if (is.na(tau)[1]) {
    if (is.na(j)) {
      AT<-t(aa%*%t(theta)) # t ( J x df  %*%  df x N )
      Bz<-t(bb%*%t(as.matrix(rep(1,nrow(theta)))))  # t ( J x df %*% df x N )
    } else {
      Bz<-bb[j]
      if (ncol(aa)>1) {
        AT<-(aa[j,]%*%t(theta))
      } else {
        AT<-aa[j]*as.vector(theta)
      }
    } 
  } else {
    if (is.na(j)) {
      AT<-array(rep(t(aa%*%t(theta)),ncol(tau)),c(nrow(theta),nrow(xi),ncol(tau))) #J x df  %*%  df x N
      B = list()
      for (i in 1:ncol(tau)) {
        B[[i]]<-t((tau[,i])%*%t(as.matrix(rep(1,nrow(theta)))))
      }
      Bz <- array(do.call(c,B),dim = c(nrow(theta),nrow(xi),ncol(tau)))
    } else {
      Bz<-matrix(rep(tau[j,],nrow(theta)),nrow(theta),ncol(tau))
      if (ncol(aa)>1) {
        AT<-matrix(rep(as.vector(aa[j,]%*%t(theta)),ncol(tau)),nrow(theta),ncol(tau))
      } else {
        AT<-matrix(rep(as.vector(aa[j]*theta),ncol(tau)),nrow(theta),ncol(tau))
      }
    }
  }
  if (!guess) {
    p<-pnorm(AT-Bz)
  } else {
    p<-C+(1-C)*pnorm(AT-Bz)
  }
  return (p)
}
