TruncNormSlow <-
function(Eta,indL=indL,indU=indU) {
  #   Ez[which(w[,j]==1),j]<-rtruncnorm( length(Ez[which(w[,j]==1),j]), a=0,  b=Inf, mean=Eta[which(w[,j]==1),j], sd=1)
  #   Ez[which(w[,j]==0),j]<-rtruncnorm( length(Ez[which(w[,j]==0),j]), a=-Inf, b=0, mean=Eta[which(w[,j]==0),j], sd=1)
  Zj<-lapply(1:ncol(Eta),function(x) {  pp <- cbind(-Inf,-Eta[,x],Inf)
                                    U  <- matrix(runif(nrow(Eta),0.000009,0.999991),nrow(Eta),1)
                                    pL <- pnorm( pp[,indL[[x]]] ) 
                                    pU <- pnorm( pp[,indU[[x]]] ) 
                                    Eta[,x]+qnorm( pLU + U*(pU - pL)  ) 
  })
  do.call("cbind",Zj)
}
