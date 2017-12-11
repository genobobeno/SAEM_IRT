WrapX <- function(j,A,b,d,theta) {
  # ker is m x (ncat-1) matrix of item kernels
  ker <- apply(theta%*%matrix(A[j,],Q,1)-b[j],1,
               function(x) x - c(d[j,]))
  Xi	= matrix(NA,n1cat,N)
  # Generate missing option propensities
  Xi[MY[[j]]] <- rnorm(missList[[j]]$miss)
  # Generate nonmissing option propensities
  r1 		<- matrix(R[[j]][MN[[j]]]  ,n1cat,missList[[j]]$mcol)
  U  		<- matrix(runif(missList[[j]]$nmis) ,n1cat,missList[[j]]$mcol)
  P		  <- matrix(pnorm(-ker)[MN[[j]]]      ,n1cat,missList[[j]]$mcol)
  Xi[MN[[j]]]  <-  qnorm( r1*U + P*(r1 + U - 2*r1*U) )
  ker + Xi
}