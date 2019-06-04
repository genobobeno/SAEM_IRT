GetLikelihood<-function(fit.data) {
  if ("tau" %in% names(fit.data)) {
    TAU<-fit.data$tau
  } else {
    TAU<-NA
  }
  p<-ProbOgive(xi=fit.data$xi,theta=fit.data$That[,1:ncol(fit.data$A)],
               tau=TAU,guess = fit.data$settings$guess) # TQ x J
  if (is.na(TAU)[1]) {
#    for (n in 1:nrow(rp)) {
     y<-as.matrix(fit.data$RP) # N x J
     L<-rowSums(y*log(p+0.0000001)+(1-y)*log(1-p+0.0000001))
     L<-sum(L)
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
  return(setNames(c(L,aic.L,bic.L,dof),c("LogL","AIC","BIC","DOF")))
}
