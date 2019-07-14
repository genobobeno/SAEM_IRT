
#######   Gamma from 1:20
gLimit=20   # High Gamma limit
Samples=50    # Number of random matrices to sample for each condition
nZ = 50   # Number of Z's corresponding to number of items
for (i in 1:gLimit) {
  N = nZ*i # Number of Students taking an exam... 
  J = nZ # Number of items on the exam... J = Random Gaussian samples (analogous to all the Z's or X's)
  Gamma = N/J  # the only parameter for the Marcenko-Pastur distribution
  EVals=do.call("rbind",lapply(1:Samples,function(y) {    
    rMat=do.call("rbind",lapply(1:N,function(x) rnorm(J)))  # The random matrix
    eigen(cov(rMat))$values
  }))
  if(i==1){
    plot(density(EVals),xlim=c(-0.5,6),ylim=c(0,2),main=paste0("gamma = 1:",gLimit," ; ",Samples," samples"))
  } else {lines(density(EVals))}
}
###################   Gamma from 0.4 to 1
Nseq=20   # Number of 
Samples=100
nZ = 50
for (i in seq(0.4,1.0,length.out = Nseq)) {
  N = ceiling(nZ*i) # Number of Students taking an exam... 
  J = nZ # Number of items on the exam... J = Random Gaussian samples (analogous to all the Z's or X's)
  Gamma = N/J  # the only parameter for the Marcenko-Pastur distribution
  EVals=do.call("rbind",lapply(1:Samples,function(y) {    
    rMat=do.call("rbind",lapply(1:N,function(x) rnorm(J)))  # The random matrix
    eigen(cov(rMat))$values
  }))
  if(i==0.4){
    plot(density(EVals),xlim=c(-0.5,6),ylim=c(0,2),main=paste0("gamma = 0.4-1.0 ; ",Samples," samples"))
  } else {lines(density(EVals))}
}


gamma = J/N
if (plot) {
  par(mfrow=c(1,1))
  plot(EV,dist,xlim=c(-0.5,6),ylim=c(0,2*sqrt(J)),type="n",
       main=paste0("gamma = ",gamma,"; ",Samples," samples"))
  lines(EV,dist,lwd=2)}




mpd1<-MP_Sample(gamma=0.01,p=40,samples=500)
mpd2<-MP_Sample(gamma=0.1,p=40,samples=500)
mpd3<-MP_Sample(gamma=0.01,p=4,samples=500)
mpd4<-MP_Sample(gamma=0.1,p=4,samples=500)

par(mfrow=c(2,2))
bks<-seq(0,2.5,length.out = 100)
hist(mpd4$ev,freq = FALSE,col = 2,breaks = bks,main=expression(paste("500 samples : ",italic(J)," = 4, ",italic(N)," = 400")),xlab=expression(lambda),ylim=c(0,3.3))
lines(mpd4$x.mp,mpd4$y.mp,lwd=2,col=4)
legend("topright",legend = expression(italic(P)[MP] ~ " : " ~ italic(c) ~ "= 0.1"),lty=1,lwd=2,col=4)
hist(mpd3$ev,freq = FALSE,col = 2,breaks = bks,main=expression(paste("500 samples : ",italic(J)," = 4, ",italic(N)," = 4000")),xlab=expression(lambda))
lines(mpd3$x.mp,mpd3$y.mp,lwd=2,col=4)
legend("topright",legend = expression(italic(P)[MP] ~ " : " ~ italic(c) ~ "= 0.01"),lty=1,lwd=2,col=4)
hist(mpd2$ev,freq = FALSE,col = 2,breaks = bks,main=expression(paste("500 samples : ",italic(J)," = 40, ",italic(N)," = 4000")),xlab=expression(lambda),ylim=c(0,3.3))
lines(mpd2$x.mp,mpd2$y.mp,lwd=2,col=4)
legend("topright",legend = expression(italic(P)[MP] ~ " : " ~ italic(c) ~ "= 0.1"),lty=1,lwd=2,col=4)
hist(mpd1$ev,freq = FALSE,col = 2,breaks = bks,main=expression(paste("500 samples : ",italic(J)," = 40, ",italic(N)," = 40000")),xlab=expression(lambda))
lines(mpd1$x.mp,mpd1$y.mp,lwd=2,col=4)
legend("topright",legend = expression(italic(P)[MP] ~ " : " ~ italic(c) ~ "= 0.01"),lty=1,lwd=2,col=4)

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
#  lines(density(EVals))}

#========================

#######   Gamma from 1:20
library(RColorBrewer)

gLimit=20   # High Gamma limit
Samples=100    # Number of random matrices to sample for each condition
nZ = 50   # Number of Z's corresponding to number of items
cols<-rainbow(gLimit,start = 0,end = 0.75)
for (i in 1:gLimit) {
  N = nZ*i # Number of Students taking an exam... 
  J = nZ # Number of items on the exam... J = Random Gaussian samples (analogous to all the Z's or X's)
  Gamma = N/J  # the only parameter for the Marcenko-Pastur distribution
  EVals=do.call("rbind",lapply(1:Samples,function(y) {    
    rMat=do.call("rbind",lapply(1:N,function(x) rnorm(J)))  # The random matrix
    eigen(cov(rMat))$values
  }))
  if(i==1){
    eval(parse(text=paste0('plot(density(EVals),xlim=c(-0.5,6),ylim=c(0,2),col=cols[1],
                           main=expression(paste(gamma," = 1-",',gLimit,'," ; ",',Samples,'," samples")))')))
  } else {lines(density(EVals),col=cols[i])}
}
eval(parse(text=paste0('legend("topright",legend = c(expression(paste(gamma," = 1")),
                       expression(paste(gamma," = ',gLimit,'"))),
                       col=c(cols[1],cols[gLimit]),lty=1)')))
###################   Gamma from 0.4 to 1
Nseq=20   # Number of 
Samples=100
nZ = 50
gBounds<-c(0.4,1.0)
for (i in seq(gBounds[1],gBounds[2],length.out = Nseq)) {
  N = ceiling(nZ*i) # Number of Students taking an exam... 
  J = nZ # Number of items on the exam... J = Random Gaussian samples (analogous to all the Z's or X's)
  Gamma = N/J  # the only parameter for the Marcenko-Pastur distribution
  EVals=do.call("rbind",lapply(1:Samples,function(y) {    
    rMat=do.call("rbind",lapply(1:N,function(x) rnorm(J)))  # The random matrix
    eigen(cov(rMat))$values
  }))
  if(i==0.4){
    nC<-1
    eval(parse(text=paste0('plot(density(EVals),xlim=c(-0.5,6),ylim=c(0,2),col=cols[nC],
                           main=expression(paste(gamma," = ",',gBounds[1],',"-",',gBounds[2],'," ; ",',Samples,'," samples")))')))
  } else {
    nC<-nC+1
    lines(density(EVals),col=cols[nC])
  }
}
eval(parse(text=paste0('legend("topright",legend = c(expression(paste(gamma," = 0.4")),expression(paste(gamma," = 1.0"))),
                       col=c(cols[1],cols[gLimit]),lty=1)')))


#install.packages("RMTstat")

library(RMTstat)
# MP_Sample
# MP_Sample<-function(gamma,p=NA,n=NA,samples) {
# if (is.na(p)) {p<-as.integer(gamma*n)} else {n<-as.integer(p/gamma)}
# gminus = (1-sqrt(gamma))^2
# gplus = (1+sqrt(gamma))^2
# EV = seq(gminus,gplus,length.out = 1000)
# dist<-(1/(2*pi*gamma*EV))*sqrt((EV-gminus)*(gplus-EV))
# EVals=do.call("rbind",lapply(1:samples,function(y) {
# rMat=do.call("rbind",lapply(1:n,function(x) rnorm(p)))  # The random matrix
# eigen(cov(rMat))$values
# }))
# max.ev<-apply(EVals,1,max)
# list(ev=EVals,x.mp=EV,y.mp=dist,max.ev=max.ev)
# }


slices=3
m<-100
mpd<-list()
mpda<-list()
for (s in 1:slices) {
  cat (":",s)
  mpda[[s]]<-MP_Sample(gamma=0.1,p=m*s,samples=500)
  mpd[[s]]<-MP_Sample(gamma=0.01,p=m*s,samples=500)
}
par(mfrow=c(1,2))
cls<-c(2,3,4)
for (s in slices:1) {
  if (s==slices) {
    leg<-eval(parse(text=paste0("c(expression(lambda[1]^MP),",paste0("expression(italic(J)[",m*1:slices,"])",collapse=","),")")))
    plot(density(mpd[[s]]$max.ev),col=cls[s],lty=s,xlab=expression(lambda[1]),xlim=c(0.975,1.015)*max(mpd[[1]]$x.mp[which(mpd[[s]]$y.mp>0)]),
         main=expression("Asymptotics :" ~ italic(c) ~ "= 0.01 : Samples = 500"),lwd=2,yaxt="n") #,freq = FALSE,col = 2,breaks = bks,main=expression(paste(gamma,"=",0.01)),xlab="largest eigenvalues")
    abline(v=max(mpd[[s]]$x.mp[which(mpd[[s]]$y.mp>0)]),lwd=2,lty=5)
    legend("topright",legend = leg,lwd=2,col=c(1,cls),lty=c(5,1:slices))
  } else {
    lines(density(mpd[[s]]$max.ev),col=cls[s],lwd=2,lty=s)
  }
}
for (s in slices:1) {
  if (s==slices) {
    leg<-eval(parse(text=paste0("c(expression(lambda[1]^MP),",paste0("expression(italic(J)[",m*1:slices,"])",collapse=","),")")))
    plot(density(mpda[[s]]$max.ev),col=cls[s],lty=s,xlab=expression(lambda[1]),xlim=c(0.925,1.05)*max(mpda[[1]]$x.mp[which(mpda[[s]]$y.mp>0)]),
         main=expression("Asymptotics :" ~ italic(c) ~ "= 0.1 : Samples = 500"),lwd=2,yaxt="n") #,freq = FALSE,col = 2,breaks = bks,main=expression(paste(gamma,"=",0.01)),xlab="largest eigenvalues")
    abline(v=max(mpda[[s]]$x.mp[which(mpda[[s]]$y.mp>0)]),lwd=2,lty=5)
    legend("topright",legend = leg,lwd=2,col=c(1,cls),lty=c(2,1:slices))
  } else {
    lines(density(mpda[[s]]$max.ev),col=cls[s],lwd=2,lty=s)
  }
}
# for (s in 1:slices) {
#   if (s==1) {
#     leg<-eval(parse(text=paste0("c(expression(lambda[1]^MP),",paste0("expression(italic(J)[",m*1:slices,"])",collapse=","),")")))
#     plot(density(mpd[[s]]$max.ev),col=cls[s],lty=1,xlab=expression(lambda[1]),xlim=c(0.975,1.025)*max(mpd[[s]]$x.mp[which(mpd[[s]]$y.mp>0)]),
#          main=expression("Asymptotics :" ~ italic(c) ~ "= 0.01 : Samples = 500")) #,freq = FALSE,col = 2,breaks = bks,main=expression(paste(gamma,"=",0.01)),xlab="largest eigenvalues")
#     abline(v=max(mpd[[s]]$x.mp[which(mpd[[s]]$y.mp>0)]),lwd=2,lty=2)
#     legend("topright",legend = leg,lty=c(2,rep(1,slices)),lwd=c(2,rep(1,slices)),col=c(1,cls))
#   } else {
#     lines(density(mpd[[s]]$max.ev),col=cls[s])
#   }
# }
# for (s in 1:slices) {
#   if (s==1) {
#     leg<-eval(parse(text=paste0("c(expression(lambda[1]^MP),",paste0("expression(italic(J)[",m*1:slices,"])",collapse=","),")")))
#     plot(density(mpda[[s]]$max.ev),col=cls[s],lty=1,xlab=expression(lambda[1]),xlim=c(0.975,1.025)*max(mpda[[s]]$x.mp[which(mpda[[s]]$y.mp>0)]),
#          main=expression("Asymptotics :" ~ italic(c) ~ "= 0.1 : Samples = 500")) #,freq = FALSE,col = 2,breaks = bks,main=expression(paste(gamma,"=",0.01)),xlab="largest eigenvalues")
#     abline(v=max(mpda[[s]]$x.mp[which(mpda[[s]]$y.mp>0)]),lwd=2,lty=2)
#     legend("topright",legend = leg,lty=c(2,rep(1,slices)),lwd=c(2,rep(1,slices)),col=c(1,cls))
#   } else {
#     lines(density(mpda[[s]]$max.ev),col=cls[s])
#   }
# }


# TracyWidom<-function(ev,p,n) {
# mu<-1/n*(sqrt(n-1/2)+sqrt(p-1/2))^2
# sig<-sqrt(mu/n)*(1/sqrt(n-1/2)+1/sqrt(p-1/2))^(1/3)
# (ev-mu)/sig
# }
# TWTransform<-function(ev,p,n) {
# mu<-1/n*(sqrt(n-1/2)+sqrt(p-1/2))^2
# sig<-sqrt(mu/n)*(1/sqrt(n-1/2)+1/sqrt(p-1/2))^(1/3)
# (ev-mu)/sig
# }

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

