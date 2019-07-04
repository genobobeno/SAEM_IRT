GetLikelihood<-function(fit.data,LR=FALSE) {
  # load("RealData/FCI_NoGuess_A2.rda")
  # fit.data<-FitDATA
  if ("tau" %in% names(fit.data)) {
    TAU<-fit.data$tau
  } else {
    TAU<-NA
  }
  p<-ProbOgive(xi=fit.data$xi,theta=fit.data$That[,1:ncol(fit.data$A)],
               tau=TAU,guess = fit.data$settings$guess) # TQ x J
  if (!LR) {
    library(ks)
    if (ncol(fit.data$A)<7) {
      d.theta<-kde(x = fit.data$That[,1:ncol(fit.data$A)],
                   xmin=rep(-1*max(abs(range(fit.data$That[,1:ncol(fit.data$A)])))-0.1,ncol(fit.data$A)),
                   xmax=rep(max(abs(range(fit.data$That[,1:ncol(fit.data$A)])))+0.1,ncol(fit.data$A)))
      if (fit.data$settings$Adim==1) {
        t.grid<-d.theta$eval.points
      } else {
        t.grid<-expand.grid(d.theta$eval.points)
      }
      d.grid<-as.vector(d.theta$estimate)
      p.w<-ProbOgive(xi=fit.data$xi,theta=t.grid,
                   tau=TAU,guess = fit.data$settings$guess)
    }
  }
  if (is.na(TAU)[1]) {
#    for (n in 1:nrow(rp)) {
     y<-as.matrix(fit.data$RP) # N x J
     L<-rowSums(y*log(p+0.0000001)+(1-y)*log(1-p+0.0000001))
     L<-sum(L)
     if (fit.data$settings$Adim<7 & !LR) {
       w.aic <- 0
       p.waic <- 0
       for (i in 1:nrow(fit.data$That)) {
         y.w<-matrix(y[i,],nrow = nrow(p.w),ncol=ncol(y),byrow = TRUE)
         # LL<-y.w*log(p.w+0.00000001)+(1-y.w)*log(1-p.w+0.00000001)
         # w.aic<-w.aic+sum(colSums(matrix(d.grid,nrow=length(d.grid),ncol=ncol(y))*(LL)))
         # p.waic<-p.waic+var(rowSums(LL))
         w.LL<-((p.w+0.00000001)^(y.w))*((1-p.w+0.00000001)^(1-y.w))
         w.aic<-w.aic+sum(log(colSums(matrix(d.grid,nrow=length(d.grid),ncol=ncol(y))*(w.LL))))
         p.waic<-p.waic+var(rowSums(log(w.LL)))
       }
       waic=-2*(w.aic-p.waic)
     }
  } else {
    y<-array(data = rep(fit.data$RP,1+ncol(TAU)),dim = c(nrow(fit.data$RP),ncol(fit.data$RP),1+ncol(TAU)))
    for (k in 1:dim(y)[3]) {
      y[,,k]<-0+(y[,,k]==(k-1))
    }
    MNParray <- array(0,dim = c(nrow(fit.data$RP),ncol(fit.data$RP),1+ncol(TAU)))
    for (i in 1:nrow(fit.data$RP)) for (j in 1:ncol(fit.data$RP)) MNParray[i,j,] <- diff(-1*c(1,p[i,j,],0))
    L<-sum(rowSums(y*log(MNParray)))
  }
  
  dof<-ncol(fit.data$xi)*nrow(fit.data$xi)
  if ("tau" %in% names(fit.data)) dof<-dof+nrow(fit.data$tau)*ncol(fit.data$tau)
  aic.L<- (-2)*L+2*dof
  bic.L<- (-2)*L+log(nrow(fit.data$RP))*dof
  if (is.na(TAU)[1] & fit.data$settings$Adim<7 & !LR) {
    return(setNames(c(L,aic.L,bic.L,waic,dof),c("LogL","AIC","BIC","WAIC","DOF")))
  } else {
    return(setNames(c(L,aic.L,bic.L,dof),c("LogL","AIC","BIC","DOF")))
  }
}
