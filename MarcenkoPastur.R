MarcenkoPastur<-function(J,N,generate=FALSE,Samples=10,plot=FALSE) {
  # N=100000
  # J=100
  # generate=TRUE
  gamma = J/N
  gminus = (1-sqrt(gamma))^2
  gplus = (1+sqrt(gamma))^2
  EV = seq(gminus,gplus,length.out = 1000)
  dist<-(1/(2*pi*gamma*EV))*sqrt((EV-gminus)*(gplus-EV))
  HighEV<-max(dist)
  if (plot) {
  par(mfrow=c(1,1))
  plot(EV,dist,xlim=c(-0.5,6),ylim=c(0,2*sqrt(J)),type="n",
       main=paste0("gamma = ",gamma,"; ",Samples," samples"))
  lines(EV,dist,lwd=2)}
  if (generate) {
    j=8
    n<-ceiling(j/gPulls[i])
    EVals=do.call("rbind",lapply(1:Samples,function(y) {    
      rMat=do.call("rbind",lapply(1:j,function(x) rnorm(n)))  # The random matrix
      eigen(cov(rMat))$values
    }))
    HighEV<-max(HighEV,max(EVals))
    rm(EVals)
    gc()
    if (plot) {
      gPulls<-seq(gamma/2,2,length.out=5)
      for (i in 1:length(gPulls)) {
        n<-ceiling(j/gPulls[i])
        EVals=do.call("rbind",lapply(1:Samples,function(y) {    
          rMat=do.call("rbind",lapply(1:j,function(x) rnorm(n)))  # The random matrix
          eigen(cov(rMat))$values
        }))
        lines(density(EVals),col=i)
        rm(EVals)
        gc()
      }
    }
  }
  HighEV
}
  
  
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
