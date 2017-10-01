DrawW <-
function(aa,bb,cc,tt,rp) {
  AT<-t(as.matrix(aa)%*%t(as.matrix(tt))) 
  Bz<-t(as.matrix(bb)%*%t(as.matrix(rep(1,nrow(rp)))))
  p<-pnorm(AT-Bz)
  W<-mat.or.vec(nrow(rp),ncol(rp))
  for (i in 1:nrow(rp)) for (j in 1:ncol(rp)) {
    if (rp[i,j]) { W[i,j]<-rbinom(1,size=1,prob=p[i,j]/(p[i,j]+cc[j]*(1-p[i,j])))
    } else { W[i,j]<-0 } 
  }
  return(W)
}
