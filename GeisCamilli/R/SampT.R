SampT <-
function(aa,bb,zz,rp,prior,prl=FALSE,cores=settings$cores) {
  N<-nrow(rp)
  J<-ncol(rp)
  if (ncol(as.matrix(aa))>1){
    A2<-t(as.matrix(aa))%*%as.matrix(aa) # + diag(ncol(as.matrix(aa)))
    #A2<-prior$sigma+t(as.matrix(aa))%*% R %*%as.matrix(aa) <- If we have a correlation 
    V<-ginv(A2 + diag(ncol(as.matrix(aa))))
  } else {
    A2<-as.numeric(t(as.matrix(aa))%*%as.matrix(aa))
    V<-1/A2
  }
  ZB<-zz+t(as.matrix(bb)%*%t(as.matrix(rep(1,N)))) # N X J  
  if (ncol(as.matrix(aa))>1) {
    #### Solve for Theta using each item and sampled Z ####
    tHat<-t(V%*%t(as.matrix(aa))%*%t(ZB))  # N X Adim
  } else {
    AT<-t(as.matrix(aa)%*%t(as.matrix(rep(1,N)))) # N X J
    tHat<-rowSums((ZB)*AT)*V
  }
  that<-mat.or.vec(N,ncol(as.matrix(aa)))
  #### Resample Theta around that mean estimate ####
  if (ncol(as.matrix(aa))==1) {
    #for (n in 1:N) that[n]<-rnorm(1,mean=(tHat[n]/V+prior$tmu/prior$tsigma)/(1/V+1/prior$tsigma),sd=sqrt(1/(1/V+1/prior$tsigma)))
    Mean<-(tHat/V+prior$tmu/as.vector(prior$tsigma))/(1/V+1/as.vector(prior$tsigma))
    SD<-sqrt(1/(1/V+1/prior$tsigma))
    that<-rnorm(N,mean=Mean,sd=SD)
  } else {
    #for (n in 1:N) that[n,]<-mvrnorm(1,tHat[n,],Sigma=V)
    if (prl) {
      p.lst<-suppressWarnings(vsplit(1:N,f=1:cores))
      #clusterExport(cl,c("p.lst","tHat","V"))
      that<-do.call(rbind,parLapply(cl,1:cores,pMVNarma,pList=p.lst,thHat=tHat,VAR=V))
    } else {
      that<-mvrnormArma(N, tHat, V)
    }
  }
  that
}
