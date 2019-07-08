WrapZ1 <- function(j,A,b,d,theta,settings=settings,indL=indL,indU=indU) {
  Zj = matrix(NA,nrow(as.matrix(theta)),1)
  if (ncol(A)>1) {
    eta  	<- theta%*%t(t(A[j,]))
  } else {
    eta  	<- theta%*%t(A[j,])
  }
  bd	<- matrix(b[j]+d[j,],nrow(as.matrix(theta)),settings$ncat-1,byrow=TRUE)
  hold <- sweep(bd,1,eta,"-")
  pp <- cbind(-Inf,hold,Inf)
  # Trim extreme values
  U  <- matrix(runif(nrow(as.matrix(theta)),0.000009,0.999991),nrow(as.matrix(theta)),1)
  # Item propensities for all examinees from truncated normal
  # yL <- matrix(y.a[,j]+1,n,1) 
  # indL <- cbind(nn,yL); 		indU <- cbind(nn,yL + 1)
  pL <- pnorm( pp[indL[[j]]] ) ; 	pU <- pnorm( pp[indU[[j]]] )
  #pL <- approxPnorm( pp[indL[[j]]] ) ; 	pU <- approxPnorm( pp[indU[[j]]] )
  Zj <- eta + qnorm( pL + U*(pU - pL)  ) 
  c(Zj)
}