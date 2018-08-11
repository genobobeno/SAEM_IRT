fastTruncNorm <-
function(j,Eta,hold) {
  #   Ez[which(w[,j]==1),j]<-rtruncnorm( length(Ez[which(w[,j]==1),j]), a=0,  b=Inf, mean=Eta[which(w[,j]==1),j], sd=1)
  #   Ez[which(w[,j]==0),j]<-rtruncnorm( length(Ez[which(w[,j]==0),j]), a=-Inf, b=0, mean=Eta[which(w[,j]==0),j], sd=1)
  pp <- cbind(-Inf,hold[,j],Inf)
  U  <- matrix(runif(N,0.000009,0.999991),N,1)
  pL <- pnorm( pp[indL[[j]]] ) ; pU <- pnorm( pp[indU[[j]]] ) 
  Eta[,j] + qnorm( pL + U*(pU - pL)  )  # Z[,j]
}
