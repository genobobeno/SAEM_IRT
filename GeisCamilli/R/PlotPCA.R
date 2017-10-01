PlotPCA <-
function(gen,fit,estfile,simfile) {
  PCA<-as.matrix(cbind(fit$A,fit$A-2*fit$EmpSE$SEA,fit$A+2*fit$EmpSE$SEA))
  for (i in 0:2) PCA[,i*fit$settings$Adim+1:fit$settings$Adim]<-PCA[,i*fit$settings$Adim+1:fit$settings$Adim]%*%fit$AR$Th
  par(mfrow=c(2,fit$settings$Adim+1))
  TEST<-matrix(rep(1,nrow(gen$XI)*(fit$settings$Adim+1)),nrow(gen$XI),(fit$settings$Adim+1))
  for (i in 1:fit$settings$Adim) {
    plot(gen$XI[,i],PCA[,i],type="n",xlim=c(0,max(PCA)),ylim=range(PCA),main=paste("PCA Empirical Error Estimation Plot, A",i,sep=""))
    points(gen$XI[,i],PCA[,i],col="blue",pch=19)
    arrows(gen$XI[,i],PCA[,2*fit$settings$Adim+i],gen$XI[,i],PCA[,fit$settings$Adim+i],code=3,angle=90,length=0.07)
    abline(a=0,b=1)    
    TEST[,i]<-(PCA[,i]-gen$XI[,i])/((PCA[,i]-PCA[,fit$settings$Adim+i])/2)
  }
  plot(gen$XI[,ncol(gen$XI)],fit$B,type="n",xlim=range(gen$XI[,ncol(gen$XI)]),ylim=1.2*range(fit$B),main="PCA Empirical Error Estimation Plot, B")
  points(gen$XI[,ncol(gen$XI)],fit$B,col="blue",pch=19)
  arrows(gen$XI[,ncol(gen$XI)],fit$B-2*fit$EmpSE$SEB,gen$XI[,ncol(gen$XI)],fit$B+2*fit$EmpSE$SEB,code=3,angle=90,length=0.07)
  abline(a=0,b=1)
  TEST[,ncol(gen$XI)]<-(fit$B-gen$XI[,ncol(gen$XI)])/fit$EmpSE$SEB
  plot(density(as.vector(TEST)),main="Density of Z scores of Variables")
  plot(1:nrow(TEST),TEST[,1],type="n",main="Z scores of each item's fitted parameters",ylim=range(TEST))
  for (i in 1:(fit$settings$Adim+1)) points(1:nrow(TEST),TEST[,i],col=i+1)
  plot(1:nrow(TEST),rowSums(TEST^2),main="sum(Z^2) of each item",ylim=range(TEST^2))
}
