fastSampZ <-
function(aa,bb,that,srp,w=NA,prl=settings$parallel) {
  Eta<-t(as.matrix(aa)%*%t(as.matrix(that))) 
  hold <- t(as.matrix(bb)%*%t(as.matrix(rep(1,nrow(srp))))) - Eta
  # Ez<-mat.or.vec(nrow(srp),ncol(srp))
  if (!is.na(w)[1]) { srp<-w } 
  #  Z = matrix(NA,N,J)
  # Item propensities for all examinees from truncated normal
  # yL <- matrix(y.a[,j]+1,n,1) 
  if (prl) {
    Z 		<- simplify2array(parSapply(cl,1:J,fastTruncNorm,Eta,hold,simplify=FALSE), 
                        higher=TRUE)
  } else {
    Z     <- TruncNormSlow(Eta=Eta,hold=hold,indL=indL,indU=indU)
  }
  
  Q
  
  
  return (Ez)
}
