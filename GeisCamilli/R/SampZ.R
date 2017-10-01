SampZ <-
function(aa,bb,that,rp,w=NA,prl=settings$parallel) {
  Eta<-t(as.matrix(aa)%*%t(as.matrix(that))) - t(as.matrix(bb)%*%t(as.matrix(rep(1,nrow(rp)))))
  Ez<-mat.or.vec(nrow(rp),ncol(rp))
  if (!is.na(w)) {
    Ez[which(w==1)]<-TruncNorm(Eta[which(w==1)],"high",prl)
    Ez[which(w==0)]<-TruncNorm(Eta[which(w==0)],"low",prl)
  } else {
    Ez[which(rp==1)]<-TruncNorm(Eta[which(rp==1)],"high",prl)
    Ez[which(rp==0)]<-TruncNorm(Eta[which(rp==0)],"low",prl)
  }
  return (Ez)
}
