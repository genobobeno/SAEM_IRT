
library(GPArotation)
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

TWTransform<-function(ev,p,n,ptest=FALSE) {
  #ev = 1.1081814; p=100; n=5000
  mu<-1/n*(sqrt(n-1/2)+sqrt(p-1/2))^2
  sig<-sqrt(mu/n)*(1/sqrt(n-1/2)+1/sqrt(p-1/2))^(1/3)
  if (!ptest) {
    (ev-mu)/sig
  } else {
    sapply((ev-mu)/sig,function(x) ptw(q = x,beta = 1))
  }
}

InverseTW<-function(prob,p,n) {
  q<-qtw(0.95,beta = 1)
  mu<-1/n*(sqrt(n-1/2)+sqrt(p-1/2))^2
  sig<-sqrt(mu/n)*(1/sqrt(n-1/2)+1/sqrt(p-1/2))^(1/3)
  q*sig+mu
}


TWTestAndPlot<-function(d,extra.dimensions=FALSE,E=10) {
  # extra.dimensions=TRUE; d="S1";r=1
  fd<-d
  source("CreateSimulationStructure.R")
  d<-fd
  if (extra.dimensions) {
    fit.dir.tw <- "TWFits"
    twd <- paste0(d,c("Plus1","Times5"))
    fitdir<-c(paste0(fit.dir,"/",d),paste0(fit.dir.tw,"/",twd))
  } else {
    fitdir<-paste0(fit.dir,"/",d)    
  }
  simdir<-paste0(gen.dir,"/",d)
  SimList<-readRDS(paste0(simdir,"/",SFileString(sim.list[[d]],gen=TRUE),"_1.rds"))
  FitList<-readRDS(paste0(fitdir[length(fitdir)],"/",SFileString(sim.list[[d]],gen=FALSE,r = 1),".rds"))
  if (d %in% c("S1","S2","S3")) {
    S<-FitList$EZZ-as.matrix(FitList$EZ)%*%t(as.matrix(FitList$EZ))
    gEV<-eigen(S, symmetric=TRUE)
    cEV<-eigen(cov2cor(S),symmetric = TRUE)
  } else {
    gEV<-eigen(FitList$EZZ,symmetric=TRUE)
    cEV<-eigen(cov2cor(FitList$EZZ),symmetric = TRUE)
  }
  if (d %in% c("S1","S2","S3")) {
    par(mfrow=c(1,2))
  } else {
    par(mfrow=c(1,2))
  }
  pCH<-matrix(c(0,2,5,15,17,18),nrow=3,ncol=2)
  Q<-1
  sTitle<-eval(parse(text=paste0("expression('Largest'~",E,"~'Eigenvalues')")))
  plot(Q:E,c(0.85,gEV$values[Q:(E-1)]),type="n",main=sTitle,xlab="number of factors",ylab="estimators",lty=1,pch=1)
  for (i in 1:length(fitdir)) {
    #for (r in 1:sim.list[[d]]$Reps) {
    FitList<-readRDS(paste0(fitdir[i],"/",SFileString(sim.list[[d]],gen=FALSE,r = 1),".rds"))
    if (d %in% c("S1","S2","S3")) {
      S<-FitList$EZZ-as.matrix(FitList$EZ)%*%t(as.matrix(FitList$EZ))
      gEV<-eigen(S, symmetric=TRUE)
      cEV<-eigen(cov2cor(S),symmetric = TRUE)
    } else {
      gEV<-eigen(FitList$EZZ,symmetric=TRUE)
      cEV<-eigen(cov2cor(FitList$EZZ),symmetric = TRUE)
    }
    gTW.TF<-TWTransform(gEV$values,p = sim.list[[d]]$J,n = sim.list[[d]]$N,ptest = T)
    print(paste0(fitdir[i],": TW transform and p-value of Q and Q+1 and Q+2"))
    print(TWTransform(gEV$values[sim.list[[d]]$Q+0:2],p = sim.list[[d]]$J,n = sim.list[[d]]$N,ptest = F))
    print(TWTransform(gEV$values[sim.list[[d]]$Q+0:2],p = sim.list[[d]]$J,n = sim.list[[d]]$N,ptest = T))
    cTW.TF<-TWTransform(cEV$values,p = sim.list[[d]]$J,n = sim.list[[d]]$N,ptest = T)
    #j.r<-sum(TW.TF>0.95)+1
    tw.999<-InverseTW(prob=0.999,p=sim.list[[d]]$J,n=sim.list[[d]]$N)
    lines(Q:E,gEV$values[Q:E],lty=1)
    points(Q:E,gEV$values[Q:E],pch=pCH[i,2],col=0)
    points(Q:E,gEV$values[Q:E],pch=ifelse(gEV$values[Q:E]>tw.999,pCH[i,2],pCH[i,1]),
           cex=ifelse(pCH[i,2]==18 & gEV$values[Q:E]>tw.999,1.5,1))
    lines(Q:E,cEV$values[Q:E],lty=4)
    points(Q:E,cEV$values[Q:E],pch=pCH[i,2],col=0)
    points(Q:E,cEV$values[Q:E],pch=ifelse(cEV$values[Q:E]>tw.999,pCH[i,2],pCH[i,1]),
           cex=ifelse(cEV$values[Q:E]>tw.999 & rep(pCH[i,2]==18,length(Q:E)),1.5,1))
  }
  abline(h=tw.999,lty=2)
  ltext =  c(expression(lambda["Q=1"]("S"[2])),expression(lambda["Q=1"]("Cor"("S"[2]))),
             expression(lambda["Q=2"]("S"[2])),expression(lambda["Q=2"]("Cor"("S"[2]))),
             expression(lambda["Q=5"]("S"[2])),expression(lambda["Q=5"]("Cor"("S"[2]))),
             expression(lambda["P=0.999"]("TW"[1])))
  strwidth(ltext)
  legend("topright",legend = ltext,text.width = max(strwidth(ltext))*1.05,
         lty=c(1,4,1,4,1,4,2),pch=c(0,0,2,2,5,5,NA))
  if (d %in% c("S1","S2","S3")) {
    Q<-ifelse(sim.list[[d]]$Q==1,2,1) 
    sTitle<-eval(parse(text=paste('expression("Eigenvalues ',Q,'thru',E,'")')))
    plot(Q:E,c(0.85,gEV$values[Q:(E-1)]),type="n",main=sTitle,
         xlab="number of factors",ylab="estimators",lty=1,pch=1)
    for (i in 1:length(fitdir)) {
      #for (r in 1:sim.list[[d]]$Reps) {
      FitList<-readRDS(paste0(fitdir[i],"/",SFileString(sim.list[[d]],gen=FALSE,r = 1),".rds"))
      if (d %in% c("S1","S2","S3")) {
        S<-FitList$EZZ-as.matrix(FitList$EZ)%*%t(as.matrix(FitList$EZ))
        gEV<-eigen(S, symmetric=TRUE)
        cEV<-eigen(cov2cor(S),symmetric = TRUE)
      } else {
        gEV<-eigen(FitList$EZZ,symmetric=TRUE)
        cEV<-eigen(cov2cor(FitList$EZZ),symmetric = TRUE)
      }
      gTW.TF<-TWTransform(gEV$values,p = sim.list[[d]]$J,n = sim.list[[d]]$N,ptest = T)
      cTW.TF<-TWTransform(cEV$values,p = sim.list[[d]]$J,n = sim.list[[d]]$N,ptest = T)
      tw.999<-InverseTW(prob=0.999,p=sim.list[[d]]$J,n=sim.list[[d]]$N)
      lines(Q:E,gEV$values[Q:E],lty=1)
      points(Q:E,gEV$values[Q:E],pch=pCH[i,2],col=0)
      points(Q:E,gEV$values[Q:E],pch=ifelse(gEV$values[Q:E]>tw.999,pCH[i,2],pCH[i,1]),
             cex=ifelse(pCH[i,2]==18 & gEV$values[Q:E]>tw.999,1.5,1))
      lines(Q:E,cEV$values[Q:E],lty=4)
      points(Q:E,cEV$values[Q:E],pch=pCH[i,2],col=0)
      points(Q:E,cEV$values[Q:E],pch=ifelse(cEV$values[Q:E]>tw.999,pCH[i,2],pCH[i,1]),
             cex=ifelse(cEV$values[Q:E]>tw.999 & rep(pCH[i,2]==18,length(Q:E)),1.5,1))
    }
    abline(h=tw.999,lty=2)
    
    ltext =  c(expression(lambda["Q=1"]("S"[2])),expression(lambda["Q=1"]("Cor"("S"[2]))),
              expression(lambda["Q=2"]("S"[2])),expression(lambda["Q=2"]("Cor"("S"[2]))),
              expression(lambda["Q=5"]("S"[2])),expression(lambda["Q=5"]("Cor"("S"[2]))),
              expression(lambda["P=0.999"]("TW"[1])))
    strwidth(ltext)
    legend("topright",legend = ltext,text.width = max(strwidth(ltext))*1.05,
           lty=c(1,4,1,4,1,4,2),pch=c(0,0,2,2,5,5,NA))
  }
  
  if (d %in% c("S1","S2","S3")) {
    par(mfrow=c(1,1))
  }
  FitList<-readRDS(paste0(fitdir[1],"/",SFileString(sim.list[[d]],gen=FALSE,r = 1),".rds"))
  if (d %in% c("S1","S2","S3")) {
    S<-FitList$EZZ-as.matrix(FitList$EZ)%*%t(as.matrix(FitList$EZ))
    gEV<-eigen(S, symmetric=TRUE)
    cEV<-eigen(cov2cor(S),symmetric = TRUE)
  } else {
    gEV<-eigen(FitList$EZZ,symmetric=TRUE)
    cEV<-eigen(cov2cor(FitList$EZZ),symmetric = TRUE)
  }
  Q<-1
  plot(Q:E,c(0.85,gEV$values[Q:(E-1)]/gEV$values[(Q+1):E]),type="n",main=expression("Adjacent Ratios of Eigenvalues"),
       xlab="number of factors",ylab="estimators",lty=1,pch=1)
  for (i in 1:length(fitdir)) {
    #for (r in 1:sim.list[[d]]$Reps) {
    FitList<-readRDS(paste0(fitdir[i],"/",SFileString(sim.list[[d]],gen=FALSE,r = 1),".rds"))
    if (d %in% c("S1","S2","S3")) {
      S<-FitList$EZZ-as.matrix(FitList$EZ)%*%t(as.matrix(FitList$EZ))
      gEV<-eigen(S, symmetric=TRUE)
      cEV<-eigen(cov2cor(S),symmetric = TRUE)
    } else {
      gEV<-eigen(FitList$EZZ,symmetric=TRUE)
      cEV<-eigen(cov2cor(FitList$EZZ),symmetric = TRUE)
    }
    print(fitdir[i])
    print(sapply(1:E,function(x) t.test(gEV$vectors[,x],mu = 0)$p.value))
    gTW.TF<-TWTransform(gEV$values,p = sim.list[[d]]$J,n = sim.list[[d]]$N,ptest = T)
    cTW.TF<-TWTransform(cEV$values,p = sim.list[[d]]$J,n = sim.list[[d]]$N,ptest = T)
    #j.r<-sum(TW.TF>0.95)+1
    gER<-gEV$values[1:(sim.list[[d]]$J-1)]/gEV$values[2:sim.list[[d]]$J]
    cER<-cEV$values[1:(sim.list[[d]]$J-1)]/cEV$values[2:sim.list[[d]]$J]
    lines(Q:E,gER[Q:E],lty=1)
    points(Q:E,gER[Q:E],pch=pCH[i,2],col=0)
    points(Q:E,gER[Q:E],pch=ifelse(gER[Q:E]==max(gER),pCH[i,2],pCH[i,1]),
           cex=ifelse(pCH[i,2]==18 & gER[Q:E]==max(gER),1.5,1))
    lines(Q:E,cER[Q:E],lty=4)
    points(Q:E,cER[Q:E],pch=pCH[i,2],col=0)
    points(Q:E,cER[Q:E],pch=ifelse(cER[Q:E]==max(cER),pCH[i,2],pCH[i,1]),
           cex=ifelse(cER[Q:E]==max(cER) & rep(pCH[i,2]==18,length(Q:E)),1.5,1))
  }
  ltext =  c(expression(ER["Q=1"]("S"[2])),expression(ER["Q=1"]("Cor"("S"[2]))),
             expression(ER["Q=2"]("S"[2])),expression(ER["Q=2"]("Cor"("S"[2]))),
             expression(ER["Q=5"]("S"[2])),expression(ER["Q=5"]("Cor"("S"[2]))))
  strwidth(ltext)
  legend("topright",legend = ltext,text.width = max(strwidth(ltext))*1.05,
         lty=c(1,4,1,4,1,4),pch=c(0,0,2,2,5,5))
}
