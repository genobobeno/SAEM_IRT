WBinomFast<-function(j,rp,p,cc) {
  W<-matrix(0,ncol=1,nrow=nrow(rp))
  W[rp[,j]==1,1]<-rbinom(n=sum(rp[,j]),size=1,prob=p[rp[,j]==1,j]/(p[rp[,j]==1,j]+cc[j]*(1-p[rp[,j]==1,j])))
  W
}