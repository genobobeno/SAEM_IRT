WrapZ <- function(j,A,b,d,theta) {
  Zj = matrix(NA,N,1)
  eta  	<- theta%*%t(t(A[j,]))
  bd	<- matrix(b[j]+d[j,],N,n1cat,byrow=TRUE)
  hold <- sweep(bd,1,eta,"-")
  pp <- cbind(-Inf,hold,Inf)
  # Trim extreme values
  U  <- matrix(runif(N,0.000009,0.999991),N,1)

  # Item propensities for all examinees from truncated normal
  # yL <- matrix(y.a[,j]+1,n,1) 
  # indL <- cbind(nn,yL); 		indU <- cbind(nn,yL + 1)
  pL <- pnorm( pp[indL[[j]]] ) ; 	pU <- pnorm( pp[indU[[j]]] )
  
  Zj <- eta + qnorm( pL + U*(pU - pL)  ) 
  c(Zj)
}