install.packages("RMTstat")
library(RMTstat)
MP_Sample<-function(gamma,p=NA,n=NA,samples) {
if (is.na(p)) {p<-as.integer(gamma*n)} else {n<-as.integer(p/gamma)}
gminus = (1-sqrt(gamma))^2
gplus = (1+sqrt(gamma))^2
EV = seq(gminus,gplus,length.out = 1000)
dist<-(1/(2*pi*gamma*EV))*sqrt((EV-gminus)*(gplus-EV))
EVals=do.call("rbind",lapply(1:samples,function(y) {
rMat=do.call("rbind",lapply(1:n,function(x) rnorm(p)))  # The random matrix
eigen(cov(rMat))$values
}))
max.ev<-apply(EVals,1,max)
list(ev=EVals,x.mp=EV,y.mp=dist,max.ev=max.ev)
}
slices=10
mpd<-list()
mpda<-list()
for (s in 1:slices) {
cat (":",s)
mpda[[s]]<-MP_Sample(gamma=0.1,p=20*s,samples=500)
mpd[[s]]<-MP_Sample(gamma=0.01,p=20*s,samples=500)
}
par(mfrow=c(1,2))
cls<-rainbow(slices,start = 0.05,end = 0.74)
for (s in slices:1) {
if (s==slices) {
leg<-eval(parse(text=paste0("c(expression(lambda[1]^MP),",paste0("expression(italic(J)[",20*1:slices,"])",collapse=","),")")))
plot(density(mpd[[s]]$max.ev),col=cls[s],lty=1,xlab=expression(lambda[1]),xlim=c(0.925,1.1)*max(mpd[[s]]$x.mp[which(mpd[[s]]$y.mp>0)]),
main=expression("Asymptotics :" ~ italic(c) ~ "= 0.01 : Samples = 500")) #,freq = FALSE,col = 2,breaks = bks,main=expression(paste(gamma,"=",0.01)),xlab="largest eigenvalues")
abline(v=max(mpd[[s]]$x.mp[which(mpd[[s]]$y.mp>0)]),lwd=2,lty=2)
legend("topright",legend = leg,lty=c(2,rep(1,slices)),lwd=c(2,rep(1,slices)),col=c(1,cls))
} else {
lines(density(mpd[[s]]$max.ev),col=cls[s])
}
}
for (s in slices:1) {
if (s==slices) {
leg<-eval(parse(text=paste0("c(expression(lambda[1]^MP),",paste0("expression(italic(J)[",20*1:slices,"])",collapse=","),")")))
plot(density(mpda[[s]]$max.ev),col=cls[s],lty=1,xlab=expression(lambda[1]),xlim=c(0.8,1.15)*max(mpda[[s]]$x.mp[which(mpda[[s]]$y.mp>0)]),
main=expression("Asymptotics :" ~ italic(c) ~ "= 0.1 : Samples = 500")) #,freq = FALSE,col = 2,breaks = bks,main=expression(paste(gamma,"=",0.01)),xlab="largest eigenvalues")
abline(v=max(mpda[[s]]$x.mp[which(mpda[[s]]$y.mp>0)]),lwd=2,lty=2)
legend("topright",legend = leg,lty=c(2,rep(1,slices)),lwd=c(2,rep(1,slices)),col=c(1,cls))
} else {
lines(density(mpda[[s]]$max.ev),col=cls[s])
}
}
for (s in 1:slices) {
if (s==1) {
leg<-eval(parse(text=paste0("c(expression(lambda[1]^MP),",paste0("expression(italic(J)[",20*1:slices,"])",collapse=","),")")))
plot(density(mpd[[s]]$max.ev),col=cls[s],lty=1,xlab=expression(lambda[1]),xlim=c(0.925,1.1)*max(mpd[[s]]$x.mp[which(mpd[[s]]$y.mp>0)]),
main=expression("Asymptotics :" ~ italic(c) ~ "= 0.01 : Samples = 500")) #,freq = FALSE,col = 2,breaks = bks,main=expression(paste(gamma,"=",0.01)),xlab="largest eigenvalues")
abline(v=max(mpd[[s]]$x.mp[which(mpd[[s]]$y.mp>0)]),lwd=2,lty=2)
legend("topright",legend = leg,lty=c(2,rep(1,slices)),lwd=c(2,rep(1,slices)),col=c(1,cls))
} else {
lines(density(mpd[[s]]$max.ev),col=cls[s])
}
}
for (s in 1:slices) {
if (s==1) {
leg<-eval(parse(text=paste0("c(expression(lambda[1]^MP),",paste0("expression(italic(J)[",20*1:slices,"])",collapse=","),")")))
plot(density(mpda[[s]]$max.ev),col=cls[s],lty=1,xlab=expression(lambda[1]),xlim=c(0.8,1.15)*max(mpda[[s]]$x.mp[which(mpda[[s]]$y.mp>0)]),
main=expression("Asymptotics :" ~ italic(c) ~ "= 0.1 : Samples = 500")) #,freq = FALSE,col = 2,breaks = bks,main=expression(paste(gamma,"=",0.01)),xlab="largest eigenvalues")
abline(v=max(mpda[[s]]$x.mp[which(mpda[[s]]$y.mp>0)]),lwd=2,lty=2)
legend("topright",legend = leg,lty=c(2,rep(1,slices)),lwd=c(2,rep(1,slices)),col=c(1,cls))
} else {
lines(density(mpda[[s]]$max.ev),col=cls[s])
}
}
library(RMTstat)
TracyWidom<-function(ev,p,n) {
mu<-1/n*(sqrt(n-1/2)+sqrt(p-1/2))^2
sig<-sqrt(mu/n)*(1/sqrt(n-1/2)+1/sqrt(p-1/2))^(1/3)
(ev-mu)/sig
}
TWTransform<-function(ev,p,n) {
mu<-1/n*(sqrt(n-1/2)+sqrt(p-1/2))^2
sig<-sqrt(mu/n)*(1/sqrt(n-1/2)+1/sqrt(p-1/2))^(1/3)
(ev-mu)/sig
}
RMTstat::dtw(x = 0.2)
xtw<-seq(-3,2,length.out = 501)
plot(xtw,dtw(xtw,1),main="TW")
xtw<-seq(-6,3,length.out = 901)
plot(xtw,dtw(xtw,1),main="TW")
plot(xtw,dtw(xtw,1),type = "n",main="TW")
lines(xtw,dtw(xtw,1))
slices
lines(density(TWTransform(mpda[[slices]],20,200)),col=2)
lines(density(TWTransform(mpda[[slices]]$max.ev,20,200)),col=2)
lines(density(TWTransform(mpda[[1]]$max.ev,20,200)),col=2)
xtw<-seq(-6,3,length.out = 901)
plot(xtw,dtw(xtw,1),type = "n",main=expression("TW"[1] ~ " : J = 20"))
lines(density(TWTransform(mpda[[1]]$max.ev,20,200)),col=cls[1])
lines(xtw,dtw(xtw,1))
plot(xtw,dtw(xtw,1),type = "n",main=expression("TW"[1] ~ " : J = 200"))
lines(density(TWTransform(mpda[[10]]$max.ev,200,2000)),col=cls[10])
lines(xtw,dtw(xtw,1))
par(mfrow=c(1,2))
xtw<-seq(-6,3,length.out = 901)
plot(xtw,dtw(xtw,1),type = "n",main=expression("TW"[1] ~ " : J = 20"))
lines(density(TWTransform(mpda[[1]]$max.ev,20,200)),col=cls[1])
lines(xtw,dtw(xtw,1))
plot(xtw,dtw(xtw,1),type = "n",main=expression("TW"[1] ~ " : J = 200"))
lines(density(TWTransform(mpda[[10]]$max.ev,200,2000)),col=cls[10])
lines(xtw,dtw(xtw,1))
par(mfrow=c(1,2))
xtw<-seq(-6,3,length.out = 901)
plot(xtw,dtw(xtw,1),type = "n",main=expression("TW"[1] ~ " : J = 20"),ylim=c(0,0.35))
lines(density(TWTransform(mpda[[1]]$max.ev,20,200)),col=cls[1])
lines(xtw,dtw(xtw,1),lwd=2,lty=2)
plot(xtw,dtw(xtw,1),type = "n",main=expression("TW"[1] ~ " : J = 200"))
lines(density(TWTransform(mpda[[10]]$max.ev,200,2000)),col=cls[10])
lines(xtw,dtw(xtw,1),lwd=2,lty=2)
TWTransform<-function(ev,p,n,ptest=FALSE) {
mu<-1/n*(sqrt(n-1/2)+sqrt(p-1/2))^2
sig<-sqrt(mu/n)*(1/sqrt(n-1/2)+1/sqrt(p-1/2))^(1/3)
if (!ptest) {
(ev-mu)/sig
} else {
ptw(q = (ev-mu)/sig,beta = 1)
}
}
TWTransform(mpda[[1]]$max.ev,20,200,TRUE)
table(quantile(TWTransform(mpda[[1]]$max.ev,20,200,TRUE),p=0.1*0:10))
quantile(TWTransform(mpda[[1]]$max.ev,20,200,TRUE),p=0.1*0:10)
tws10<-MP_Sample(gamma=0.1,p=10,samples=500)
tws100<-MP_Sample(gamma=0.01,p=100,samples=500)
par(mfrow=c(1,2))
xtw<-seq(-6,3,length.out = 901)
plot(xtw,dtw(xtw,1),type = "n",main=expression("TW"[1] ~ " : " ~ italic(J) ~ " = 10, " ~ italic(N) ~" = 100"),ylim=c(0,0.35))
lines(density(TWTransform(tws10$max.ev,10,100)),col=2)
lines(xtw,dtw(xtw,1),lwd=2,lty=2)
legend("topright",legend = c(expression("TW"[1]),expression(lambda[1])),lty=c(2,1),lwd=c(2,1),col=c(1,2))
plot(xtw,dtw(xtw,1),type = "n",main=expression("TW"[1] ~ " : " ~ italic(J) ~ " = 100, " ~ italic(N) ~" = 10000"))
lines(density(TWTransform(tws100$max.ev,100,10000)),col=2)
lines(xtw,dtw(xtw,1),lwd=2,lty=2)
legend("topright",legend = c(expression("TW"[1]),expression(lambda[1])),lty=c(2,1),lwd=c(2,1),col=c(1,2))
par(mfrow=c(1,2))
xtw<-seq(-6,3,length.out = 901)
plot(xtw,dtw(xtw,1),type = "n",main=expression("TW"[1] ~ " : " ~ italic(J) ~ " = 10, " ~ italic(N) ~" = 100"),
ylim=c(0,0.35),xlab="s",ylab="density")
lines(density(TWTransform(tws10$max.ev,10,100)),col=2)
lines(xtw,dtw(xtw,1),lwd=2,lty=2)
legend("topright",legend = c(expression("TW"[1]),expression(lambda[1])),lty=c(2,1),lwd=c(2,1),col=c(1,2))
plot(xtw,dtw(xtw,1),type = "n",main=expression("TW"[1] ~ " : " ~ italic(J) ~ " = 100, " ~ italic(N) ~" = 10000"),
ylim=c(0,0.35),xlab="s",ylab="density")
lines(density(TWTransform(tws100$max.ev,100,10000)),col=2)
lines(xtw,dtw(xtw,1),lwd=2,lty=2)
legend("topright",legend = c(expression("TW"[1]),expression(lambda[1])),lty=c(2,1),lwd=c(2,1),col=c(1,2))
par(mfrow=c(1,2))
xtw<-seq(-6,3,length.out = 901)
plot(xtw,dtw(xtw,1),type = "n",main=expression("TW"[1] ~ " : " ~ italic(J) ~ " = 10, " ~ italic(N) ~" = 100 : 500 samples"),
ylim=c(0,0.35),xlab="s",ylab="density")
lines(density(TWTransform(tws10$max.ev,10,100)),col=2)
lines(xtw,dtw(xtw,1),lwd=2,lty=2)
legend("topright",legend = c(expression("TW"[1]),expression(lambda[1])),lty=c(2,1),lwd=c(2,1),col=c(1,2))
plot(xtw,dtw(xtw,1),type = "n",main=expression("TW"[1] ~ " : " ~ italic(J) ~ " = 100, " ~ italic(N) ~" = 10000 : 500 samples"),
ylim=c(0,0.35),xlab="s",ylab="density")
lines(density(TWTransform(tws100$max.ev,100,10000)),col=2)
lines(xtw,dtw(xtw,1),lwd=2,lty=2)
legend("topright",legend = c(expression("TW"[1]),expression(lambda[1])),lty=c(2,1),lwd=c(2,1),col=c(1,2))
tws10<-MP_Sample(gamma=0.1,p=10,samples=500)
tws100<-MP_Sample(gamma=0.01,p=100,samples=500)
par(mfrow=c(1,2))
xtw<-seq(-6,3,length.out = 901)
plot(xtw,dtw(xtw,1),type = "n",main=expression("TW"[1] ~ " : " ~ italic(J) ~ " = 10, " ~ italic(N) ~" = 100 : 500 samples"),
ylim=c(0,0.35),xlab="s",ylab="density")
lines(density(TWTransform(tws10$max.ev,10,100)),col=2)
lines(xtw,dtw(xtw,1),lwd=2,lty=2)
legend("topright",legend = c(expression("TW"[1]),expression(lambda[1])),lty=c(2,1),lwd=c(2,1),col=c(1,2))
plot(xtw,dtw(xtw,1),type = "n",main=expression("TW"[1] ~ " : " ~ italic(J) ~ " = 100, " ~ italic(N) ~" = 10000 : 500 samples"),
ylim=c(0,0.35),xlab="s",ylab="density")
lines(density(TWTransform(tws100$max.ev,100,10000)),col=2)
lines(xtw,dtw(xtw,1),lwd=2,lty=2)
legend("topright",legend = c(expression("TW"[1]),expression(lambda[1])),lty=c(2,1),lwd=c(2,1),col=c(1,2))
savehistory("C:/Dev/R_Code/SAEM_IRT/TWCode.Rhistory")
