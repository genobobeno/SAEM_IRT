WrapX1 <- function(j,A,b,d,theta,settings=settings,R=R,MN=MN,missList=missList,MY=MY) {
  # ker is m x (ncat-1) matrix of item kernels
  if (ncol(as.matrix(A))>1) {
    # ker <- apply(theta%*%matrix(A[j,],settings$Adim,1)-b[j],1,
    #              function(x) x - c(d[j,]))
    ker <- matrix(theta%*%matrix(A[j,],settings$Adim,1)-b[j],nrow=length(d[j,]),ncol=nrow(theta),byrow=T)-
      matrix(d[j,],nrow=length(d[j,]),ncol=nrow(theta))
  } else {
    ker <- sapply(theta*A[j]-b[j],function(x) x - c(d[j,]))
  }
  #print(dim(ker))
  Xi	= matrix(NA,settings$ncat-1,nrow(as.matrix(theta)))
  # Generate missing option propensities
  Xi[MY[[j]]] <- rnorm(missList[[j]]$miss)
  # Generate nonmissing option propensities
  r1 		<- matrix(R[[j]][MN[[j]]]  ,settings$ncat-1,missList[[j]]$mcol)
  U  		<- matrix(runif(missList[[j]]$nmis) ,settings$ncat-1,missList[[j]]$mcol)
  P		  <- matrix(pnorm(-ker)[MN[[j]]]      ,settings$ncat-1,missList[[j]]$mcol)
  #P		  <- matrix(approxPnorm(-ker)[MN[[j]]]      ,settings$ncat-1,missList[[j]]$mcol)
  Xi[MN[[j]]]  <-  qnorm( r1*U + P*(r1 + U - 2*r1*U) )
  ker + Xi
}