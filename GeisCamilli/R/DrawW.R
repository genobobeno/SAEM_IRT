DrawW <-
function(aa,bb,cc,tt,rp,prl=settings$parallel) {
  p<-pnorm(t(as.matrix(aa)%*%t(as.matrix(tt)))- 
           t(as.matrix(bb)%*%t(as.matrix(rep(1,nrow(rp))))))
  if (prl) {
    W<-simplify2array(parSapply(cl,1:ncol(rp),WBinomFast,rp,p,cc,simplify=FALSE), 
                   higher=FALSE)
  } else {
    W<-mat.or.vec(nrow(rp),ncol(rp))
    for (j in 1:ncol(rp)) {
      W[rp[,j]==1,j]<-rbinom(n=sum(rp[,j]),size=1,prob=p[rp[,j]==1,j]/(p[rp[,j]==1,j]+cc[j]*(1-p[rp[,j]==1,j])))
    }
  }
  W
}
