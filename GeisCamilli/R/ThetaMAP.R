ThetaMAP <-
function(aa,bb,cc,rp,settings) {
  xi<-as.matrix(cbind(aa,bb))
  if (settings$Adim==1) {
    tgrid<-seq(-3.3,3.3,length.out=761)
    p<-ProbOgive(xi=xi,theta=tgrid) # TQ x J
    THAT<-vector()
    for (n in 1:nrow(rp)) {
      y<-t(as.matrix(rp[n,])%*%t(as.matrix(rep(1,nrow(as.matrix(tgrid)))))) # TQ x J
      THAT<-c(THAT,tgrid[which.max(apply(y*log(p)+(1-y)*log(1-p),1,sum))])
    }
  } else {
    tgrid<-seq(-3.3,3.3,length.out=761)
    THAT<-mat.or.vec(nrow(rp),settings$Adim)
    for (n in 1:nrow(rp)) {
      xy<-rep(1,settings$Adim)
      xy0<-rep(0,settings$Adim)
      It<-0
      tsearch<-mat.or.vec(length(tgrid),settings$Adim)
      while(sum(xy!=xy0)>0 | It>500) {
        xy0<-xy
        for (q in 1:settings$Adim) {
          tsearch[,q]<-tgrid
          p<-ProbOgive(xi=xi,theta=tsearch) # TQ x J
          y<-t(as.matrix(rp[n,])%*%t(as.matrix(rep(1,nrow(as.matrix(tgrid)))))) # TQ x J
          xy[q]<-tgrid[which.max(apply(y*log(p)+(1-y)*log(1-p),1,sum))]
          tsearch[,q]<-xy[q]
          if (!prod(is.finite(xy))) xy<-xy0+rep(0.1,settings$Adim)
        }
        #if (n==1) print(tsearch[1:5,1:settings$Adim])
        It<-It+1
        #if (It>500) print(paste("MD MAP failed for Examinee:",n))
        #print(xy)
      }
      THAT[n,]<-xy
    }
  } 
  return(THAT)
}
