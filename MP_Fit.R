###Fitting to Marchenko-Pastur

basedir<-getwd()
source(paste0(basedir,"/","CreateSimulationStructure.R"))
d<-"S7"
basedir<-"C:/Dev/MyFiles/NextCloud/Documents/Rutgers/SAEM_IRT/"
simdir<-paste0(basedir,"/",gen.dir,"/",d)
fitdir<-paste0(basedir,"/",fit.dir,"/",d)
r=1
FitList<-readRDS(paste0(fitdir,"/",SFileString(sim.list[[d]],gen=FALSE,r = r),".rds"))
load(paste0(fitdir,"/",SFileString(sim.list[[d]],gen=FALSE,r = 1),".rda"))  

MPFit<-function(ZZ,RP,multidimensional=TRUE,polytomous=TRUE,Z=NA) {
  #ZZ = FitList$EZZ;multidimensional=TRUE;polytomous=TRUE;Z=NA;RP=FitList$RP
  if (!multidimensional) {
    if (polytomous) {
      S<-ZZ #-t(as.matrix(FitList$EZ))%*%as.matrix(FitList$EZ)
    } else {
      if (is.na(Z)) { print("Need the Z argument"); return(NULL)}
      S<-ZZ-as.matrix(Z)%*%t(as.matrix(Z))
    }
    gEV<-eigen(S, symmetric=TRUE)
    cEV<-eigen(cov2cor(S),symmetric = TRUE)
  } else {
    gEV<-eigen(ZZ,symmetric=TRUE)
    cEV<-eigen(cov2cor(ZZ),symmetric = TRUE)
  }
  ## x=0.01*1:300;gamma=0.2;sigma=1; x=MPdens$x
  MPFunc<-function(x,gamma,sigma) { 
    (1*(x>(1-sqrt(gamma))^2 & x<(1+sqrt(gamma))^2)/(2*pi*sigma*gamma*x)) * 
      sqrt((x>(1-sqrt(gamma))^2 & x<(1+sqrt(gamma))^2) *
             (x-sigma*((1-sqrt(gamma))^2))*(sigma*((1+sqrt(gamma))^2)-x))
  }
  
  MPdens<-KernSmooth::bkde(gEV$values,bandwidth = 0.1,range.x = c(-0.13,3.05))
  mi<-mat.or.vec(length(dg),length(ds))
  dg<-c(-0.00005,-0.00002,-0.00001,-0.000005,-0.000002,-0.000001,0,0.000001,0.000002,0.000005,0.00001,0.00002,0.00005)
  ds<-c(-0.05,-0.02,-0.01,-0.005,-0.002,-0.001,0,0.001,0.002,0.005,0.01,0.02,0.05)
  gammaT<-dim(RP)[2]/dim(RP)[1]; sigmaT<-1; i<-1
  while(i>0.001) {
    i0<-i1
    for (j in 1:length(dg)) { #j=1; xx=0.01
      mi[j,]<-sapply(ds,function(xx) sum((MPdens$y-ifelse(!is.nan(MPFunc(MPdens$x,gamma=gammaT+dg[j],sigma = sigmaT+xx)),
                                                                  MPFunc(MPdens$x,gamma=gammaT+dg[j],sigma = sigmaT+xx),
                                                                  0))^2))
    }
    mi[is.nan(mi)]<-10000
    change<-which(abs(mi)==min(abs(mi)),arr.ind=TRUE)
    gammaT<-gammaT+dg[change[1]];sigmaT<-sigmaT+ds[change[2]]    
    i1<-min(abs(mi))
    i<-abs(i1-i0)
    print(paste(gammaT,sigmaT,i))
  }

  print("Covariance(S), Fit to Marchenko-Pastur")
  print(paste("gamma =",gammaT,":: sigma =",sigmaT,":: Dimensions =",gammaT*dim(RP)[1]))
  
  plot(density(gEV$values),col=2) + lines(MPdens$x,MPFunc(MPdens$x,gamma=gammaT,sigma = sigmaT))
  
  MPdens<-KernSmooth::bkde(cEV$values,bandwidth = 0.1,range.x = c(-0.13,3.05))
  mi<-mat.or.vec(length(dg),length(ds))
  gammaT<-(1-sqrt(min(cEV$values)))^2; sigmaT<-1; i<-1
  while(i>0.001) {
    i0<-i1
    for (j in 1:length(dg)) { #j=1; xx=0.01
      mi[j,]<-sapply(ds,function(xx) sum((MPdens$y-ifelse(!is.nan(MPFunc(MPdens$x,gamma=gammaT+dg[j],sigma = sigmaT+xx)),
                                                          MPFunc(MPdens$x,gamma=gammaT+dg[j],sigma = sigmaT+xx),
                                                          0))^2))
    }
    mi[is.nan(mi)]<-10000
    change<-which(abs(mi)==min(abs(mi)),arr.ind=TRUE)
    gammaT<-gammaT+dg[change[1]];sigmaT<-sigmaT+ds[change[2]]    
    i1<-min(abs(mi))
    i<-abs(i1-i0)
    print(paste(gammaT,sigmaT,i))
    plot(density(cEV$values),col=2) + lines(MPdens$x,MPFunc(MPdens$x,gamma=gammaT,sigma = sigmaT))
  }
  print("Cor(S), Fit to Marchenko-Pastur")
  print(paste("gamma =",gammaT,":: sigma =",sigmaT,":: Dimensions =",gammaT*dim(RP)[1]))

  plot(density(cEV$values),col=2) + lines(MPdens$x,MPFunc(MPdens$x,gamma=gammaT,sigma = sigmaT))
}

par(mfrow=c(2,1))

MPFit(FitList$EZZ,RP = FitList$RP)
n=nrow(FitList$RP)
p=ncol(FitList$RP)
rMat=do.call("rbind",lapply(1:n,function(x) rnorm(p)))
dmp<-density(eigen(cov2cor(cov(rMat)))$values)
plot(density(eigen(cov2cor(cov(rMat)))$values))
lines(dmp$x,MPFunc(dmp$x,gamma = p/n,sigma = 1))
lines(density(eigen(cov(rMat))$values))
