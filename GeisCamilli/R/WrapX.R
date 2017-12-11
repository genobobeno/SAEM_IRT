WrapX <- function(j,A,b,d,theta,refList) {
  n1cat<-ncol(d)
  # ker is m x (ncat-1) matrix of item kernels
  ker <- apply(theta%*%matrix(A[j,],Q,1)-b[j],1,
               function(x) x - c(d[j,]))
  Xi	= matrix(NA,n1cat,N)
  # Generate missing option propensities
  Xi[refList$MY[[j]]] <- rnorm(refList$missList[[j]]$miss)
  # Generate nonmissing option propensities
  r1 		<- matrix(refList$R[[j]][refList$MN[[j]]]  ,n1cat,refList$missList[[j]]$mcol)
  U  		<- matrix(runif(refList$missList[[j]]$nmis) ,n1cat,refList$missList[[j]]$mcol)
  P		  <- matrix(pnorm(-ker)[refList$MN[[j]]]      ,n1cat,refList$missList[[j]]$mcol)
  Xi[refList$MN[[j]]]  <-  qnorm( r1*U + P*(r1 + U - 2*r1*U) )
  ker + Xi
}