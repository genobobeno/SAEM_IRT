ThetaMAP <-
function(aa,bb,cc,rp,settings,TAU=NA) {
  xi<-as.matrix(cbind(aa,bb))
  if (settings$Adim==1) {
    tgrid<-seq(-3.3,3.3,length.out=761)
    p<-ProbOgive(xi=xi,theta=tgrid,tau=TAU) # TQ x J
    THAT<-vector()
    if (is.na(TAU)[1]) {
      for (n in 1:nrow(rp)) {
        y<-t(as.matrix(rp[n,])%*%t(as.matrix(rep(1,nrow(as.matrix(tgrid)))))) # TQ x J
        THAT<-c(THAT,tgrid[which.max(apply(y*log(p)+(1-y)*log(1-p),1,sum))])
      }
    } else {
      MNParray <- array(0,dim = c(length(tgrid),length(bb),ncol(TAU)+1))
      for (i in 1:length(tgrid)) for (j in 1:length(bb)) MNParray[i,j,] <- diff(-1*c(1,p[i,j,],0))
      #    RP[i,j]<-which(rmultinom(n = 1,size=1,prob=diff(-1*c(1,P[i,j,],0)))==1)
      for (n in 1:nrow(rp)) {
        y<-t(0+sapply(rp[n,],function(x) (0:ncol(TAU))==x))   # J X K
        yA<-array(0,dim=c(length(tgrid),length(bb),ncol(TAU)+1))
        for (tg in 1:length(tgrid))  yA[tg,,]<-y
        
         #y<-t(as.matrix(rp[n,])%*%t(as.matrix(rep(1,nrow(as.matrix(tgrid)))))) # TQ x J
        THAT<-c(THAT,tgrid[which.max(rowSums(yA*MNParray))])
      }
    }
  } else {
    tgrid<-seq(-3.3,3.3,length.out=761)
    THAT<-mat.or.vec(nrow(rp),settings$Adim)
    print(paste0("Iterating Over ",settings$Adim,"D Grid For Each of ",nrow(rp),
                 " examinees, Max: 500 iterations each"))
    if (is.na(TAU)[1]) {
      for (n in 1:nrow(rp)) {
        xy<-rep(1,settings$Adim)
        xy0<-rep(0,settings$Adim)
        It<-0
        tsearch<-mat.or.vec(length(tgrid),settings$Adim)
        while(sum(xy!=xy0)>0 | It>500) {
          xy0<-xy
          if (i%%10==1) cat(".")
          if (i%%100==1) cat(":")
          if (i%%500==1) cat("\n",i,"\t : ")
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
    } else {
      MNParray <- array(0,dim = c(length(tgrid),length(bb),ncol(TAU)+1))
      for (n in 1:nrow(rp)) {
        y<-t(0+sapply(rp[n,],function(x) (0:ncol(TAU))==x))   # J X K
        yA<-array(0,dim=c(length(tgrid),length(bb),ncol(TAU)+1))
        for (tg in 1:length(tgrid))  yA[tg,,]<-y
        xy<-rep(1,settings$Adim)
        xy0<-rep(0,settings$Adim)
        It<-0
        tsearch<-mat.or.vec(length(tgrid),settings$Adim)
        while(sum(xy!=xy0)>0 | It>500) {
          xy0<-xy
          for (q in 1:settings$Adim) {
            tsearch[,q]<-tgrid
            p<-ProbOgive(xi=xi,theta=tsearch,tau=TAU) # TQ x J
            for (i in 1:length(tgrid)) for (j in 1:length(bb)) MNParray[i,j,] <- diff(-1*c(1,p[i,j,],0))
            xy[q]<-tgrid[which.max(rowSums(yA*MNParray))]
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
  }
  THAT
}
