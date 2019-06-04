SampZFast<-
function(aa,bb,that,indL,indU,srp,w=NA,prl=FALSE) {
  Eta<-t(as.matrix(aa)%*%t(as.matrix(that))) - t(as.matrix(bb)%*%t(as.matrix(rep(1,nrow(srp)))))
  #hold<- Eta - t(as.matrix(bb)%*%t(as.matrix(rep(1,nrow(srp)))))
  # Ez<-mat.or.vec(nrow(srp),ncol(srp))
  if (!is.na(w)[1]) { srp<-w } 
  #  Z = matrix(NA,N,J)
  # Item propensities for all examinees from truncated normal
  # yL <- matrix(y.a[,j]+1,n,1) 
  # Return Z
  if (prl) {
    simplify2array(parSapply(cl,1:ncol(srp),TruncNormFast,Eta,simplify=FALSE), 
                        higher=FALSE)
  } else {
    TruncNormSlow(Eta=Eta,indL=indL,indU=indU)
  }
}
