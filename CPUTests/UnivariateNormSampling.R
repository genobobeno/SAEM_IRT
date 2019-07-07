source("InitializeGeisCamilli.R")


iter = 25
powers<-4:8

PNorm<-function(x) {
#  v<-lt[[x]]
  rnorm(length(lt[[x]]),lt[[x]],1)
}
PMean<-function(x) {
  lt[[x]]
}
vsplit<-function(x,f) {
  setNames(lapply(f,function (y) x[ceiling((1:length(x))/(length(x)/max(f)))==y]),f)
}
enableJIT(1)

t0<-c();t1<-c();t2<-c()
for (i in powers) {
  settings<-list()
  cores<-settings$cores<-4
  N<-round(10^i) 
  Mean<-rnorm(N,sd=0.7)
  SD<-1
  cl <<- parallel::makePSOCKcluster(settings$cores)
  parallel::clusterSetRNGStream(cl, iseed = round(runif(settings$cores)*1001))
  # lt<-data.frame(m = Mean, l = ceiling((1:length(Mean))/(length(Mean)/cores)))
  lt<-suppressWarnings(vsplit(Mean,1:cores))
  clusterExport(cl,c("lt"))
  t0<-c(t0,system.time(
    for (i in 1:iter) that<-unlist(parLapply(cl,1:cores,PNorm))
  )["elapsed"])

  t2<-c(t2,system.time(
    for (i in 1:iter) {
      lt<-suppressWarnings(vsplit(Mean,1:cores))
      clusterExport(cl,c("lt"))
      that<-unlist(parLapply(cl,1:cores,PNorm))}
  )["elapsed"])
  
  parallel::stopCluster(cl)
  
  t1<-c(t1,system.time(
    for (i in 1:iter) that<-rnorm(N,mean=Mean,sd=SD)
  )["elapsed"])
  
}
# 1: CL.time =   ; N.time = 
### t0 is 4 CPUs rnorm without vsplit
### t2 is 4 CPUs rnorm WITH vsplit
### t1 is 1 CPU
par(mfrow=c(1,2))
plot(range(powers),log(c(0.01,1000)),main="Univariate Normal Sampling CPU Time (25X)", 
     xlab=expression(10^x),ylab="seconds",type="n",yaxt="n")
axis(2,at = log(c(0.01,0.1,1,10,100,1000)),labels = c(0.01,0.1,1,10,100,1000),las=2)
points(powers,log(t0),pch=0)
lines(powers,log(t0),lty=2)
points(powers,log(t1),pch=2)
lines(powers,log(t1),lty=3)
points(powers,log(t2),pch=15)
lines(powers,log(t2),lty=1)
legend("topleft",legend = c("1 CPU","4 CPUs","4 CPUs : total"),pch=c(2,0,15),lty=c(3,2,1))

plot(range(powers),range(c(t2-t1,t0-t1)),main="Time Difference (25X)", 
     xlab=expression(10^x),ylab="seconds",type="n",yaxt="n",ylim=c(-60,600))
axis(2,at = c(100*c(0:5)),labels = c(100*c(0:5)),las=2)
points(powers,t0-t1,pch=0)
lines(powers,t0-t1,lty=2)
points(powers,t2-t1,pch=15)
lines(powers,t2-t1,lty=1)
legend("topleft",legend = c(expression(Delta*t("4,1")),expression(Delta*t("4:total,1"))),pch=c(0,15),lty=c(2,1))



#### MV Norm
settings<-list()
cores<-settings$cores<-4
settings$Adim<-3
mt0<-c();mt1<-c();mt2<-c()
for (i in powers) {
  N<-round(10^i) 
  tHat<-matrix(mvrnorm(N,mu = rnorm(settings$Adim),Sigma = diag(settings$Adim)),
               nrow=N,ncol=settings$Adim,byrow=TRUE)
  V<-diag(settings$Adim)
  p.lst<-suppressWarnings(vsplit(1:N,f=1:cores))
  cl <<- parallel::makePSOCKcluster(settings$cores)
  clusterCall(cl, cppInit)
  parallel::clusterSetRNGStream(cl, iseed = round(runif(settings$cores)*1001))
  clusterExport(cl,"pMVNarma")
  mt0<-c(mt0,system.time(
    for (i in 1:iter) that<-do.call(rbind,parLapply(cl,1:cores,pMVNarma,pList=p.lst,thHat=tHat,VAR=V))
  )["elapsed"])
  dim(that)
  mt2<-c(mt2,system.time(
    for (i in 1:iter) {
      p.lst<-suppressWarnings(vsplit(1:N,f=1:cores))
      for (i in 1:iter) that<-do.call(rbind,parLapply(cl,1:cores,pMVNarma,pList=p.lst,thHat=tHat,VAR=V))
    }
  )["elapsed"])
  
  parallel::stopCluster(cl)
  mt1<-c(mt1,system.time(
    for (i in 1:iter) that<-mvrnormArma(N, tHat, V)
  )["elapsed"])
}

### t0 is 4 CPUs rnorm without vsplit
### t2 is 4 CPUs rnorm WITH vsplit
### t1 is 1 CPU
par(mfrow=c(1,2))
plot(range(powers[-5]),log(c(0.01,10000)),main="Multivariate Normal Sampling CPU Time (25X)", 
     xlab=expression(10^x),ylab="seconds",type="n",yaxt="n",xaxt="n")
axis(1,at = 4:7,labels = 4:7)
axis(2,at = log(c(0.01,0.1,1,10,100,1000,10000)),labels = c(0.01,0.1,1,10,100,1000,10000),las=2)
points(powers[-5],log(mt0),pch=0)
lines(powers[-5],log(mt0),lty=2)
points(powers[-5],log(mt1),pch=2)
lines(powers[-5],log(mt1),lty=3)
points(powers[-5],log(mt2),pch=15)
lines(powers[-5],log(mt2),lty=1)
legend("topleft",legend = c("1 CPU","4 CPUs","4 CPUs : total"),pch=c(2,0,15),lty=c(3,2,1))

plot(range(powers[-5]),range(c(mt2-mt1,mt0-mt1)),main="Time Difference (25X)", 
     xlab=expression(10^x),ylab="seconds",type="n",yaxt="n",xaxt="n",ylim=c(-60,6000))
axis(1,at = 4:7,labels = 4:7)
axis(2,at = c(1000*c(0:5)),labels = c(1000*c(0:5)),las=2)
points(powers[-5],mt0-mt1,pch=0)
lines(powers[-5],mt0-mt1,lty=2)
points(powers[-5],mt2-mt1,pch=15)
lines(powers[-5],mt2-mt1,lty=1)
legend("topleft",legend = c(expression(Delta*t("4,1")),expression(Delta*t("4:total,1"))),pch=c(0,15),lty=c(2,1))
