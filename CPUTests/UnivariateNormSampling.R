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
par(mfrow=c(1,2))
plot(range(powers),log(c(0.01,400)),main="Normal Sampling CPU Time (25X)", 
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
     xlab=expression(10^x),ylab="seconds",type="n",yaxt="n",ylim=c(-60,300))
axis(2,at = c(50*c(-1:5)),labels = c(50*c(-1:5)),las=2)
points(powers,t0-t1,pch=0)
lines(powers,t0-t1,lty=2)
points(powers,t2-t1,pch=15)
lines(powers,t2-t1,lty=1)
legend("topleft",legend = c(expression(Delta*t("4,1")),expression(Delta*t("4 : total,1"))),pch=c(0,15),lty=c(2,1))

#### MV Norm

settings<-list()
cores<-settings$cores<-4
N<-10001
tHat<-as.matrix(cbind(rnorm(N,sd=0.7),rnorm(N,mean=2,sd=0.7)))
V<-diag(2)
p.lst<-suppressWarnings(vsplit(1:N,f=1:cores))
cl <<- parallel::makePSOCKcluster(settings$cores)
parallel::clusterSetRNGStream(cl, iseed = round(runif(settings$cores)*1001))
clusterExport(cl,c("p.lst","tHat","V"))
clusterExport(cl,"mvrnormArma")
system.time(
  that<-do.call(rbind,parLapply(cl,1:cores,pMVNarma))
)
dim(that)
parallel::stopCluster(cl)

system.time(
  that<-mvrnormArma(N, tHat, V)
)

