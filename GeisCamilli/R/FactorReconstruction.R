FactorReconstruction<-function(condition,repl=NA,all.reps=T,basedir="~/ParSAEM/SAEM_IRT/",ThetaFix=FALSE,...) {
  #condition="S3";repl=NA;all.reps=T;basedir=getwd()
  source(paste0(basedir,"/","CreateSimulationStructure.R"))
  library(Hmisc)
  d<-condition  
  simdir<-paste0(basedir,"/",gen.dir,"/",d)
  fitdir<-paste0(basedir,"/",fit.dir,"/",d)
  # if (!all.reps && !is.na(repl)) { #r=1
  #   SimList<-readRDS(paste0(simdir,"/",SFileString(sim.list[[d]],gen=TRUE),"_",r,".rds"))
  #   FitList<-readRDS(paste0(fitdir,"/",SFileString(sim.list[[d]],gen=FALSE,r = r),".rds"))
  #   par(mfrow=c(1,2))
  #   plot(SimList$gen.xi[,1:sim.list[[d]]$Q],FitList$xi[,1:sim.list[[d]]$Q],
  #        main=paste("Condition",d,": Replication",repl,": Slopes"),xlab=expression(italic(A)),
  #        ylab=expression(hat(italic(A))))
  #   plot(SimList$gen.xi[,sim.list[[d]]$Q+1],FitList$xi[,sim.list[[d]]$Q+1],
  #        main=paste("Condition",d,": Replication",repl,": Intercepts"),xlab=expression(italic(b)),
  #        ylab=expression(hat(italic(b))))
  # } else {
  gThat<-matrix(0, nrow=sim.list[[d]]$N*sim.list[[d]]$Reps,ncol=sim.list[[d]]$Q)
  That<-matrix(0, nrow=sim.list[[d]]$N*sim.list[[d]]$Reps,ncol=sim.list[[d]]$Q)
  dThat<-matrix(0, nrow=sim.list[[d]]$N*sim.list[[d]]$Reps,ncol=sim.list[[d]]$Q)
  FThat<-matrix(0, nrow=sim.list[[d]]$N*sim.list[[d]]$Reps,ncol=sim.list[[d]]$Q)
  dFThat<-matrix(0, nrow=sim.list[[d]]$N*sim.list[[d]]$Reps,ncol=sim.list[[d]]$Q)
  gThat<-as.data.frame(gThat)
  That<-as.data.frame(That)
  dThat<-as.data.frame(dThat)
  FThat<-as.data.frame(FThat)
  dFThat<-as.data.frame(dFThat)
  colnames(That)<-colnames(gThat)<-colnames(dThat)<-colnames(FThat)<-colnames(dFThat)<-paste0("T",1:sim.list[[d]]$Q)
  SimList<-readRDS(paste0(simdir,"/",SFileString(sim.list[[d]],gen=TRUE),"_1.rds"))
  FitList<-readRDS(paste0(fitdir,"/",SFileString(sim.list[[d]],gen=FALSE,r = 1),".rds"))
  gThat[1:sim.list[[d]]$N,]<-SimList$gen.theta
  if (sim.list[[d]]$Q>1) {
    vA<-colSums(FitList$xi[,1:sim.list[[d]]$Q])
    That[1:sim.list[[d]]$N,]<-FitList$Tmap[,1:sim.list[[d]]$Q]
    dThat[1:sim.list[[d]]$N,]<-FitList$Tmap[,1:sim.list[[d]]$Q]-SimList$gen.theta
  } else {
    vA<-sum(FitList$xi[,1])
    That[1:sim.list[[d]]$N,]<-FitList$Tmap
    dThat[1:sim.list[[d]]$N,]<-FitList$Tmap-SimList$gen.theta
  }
  if (sim.list[[d]]$Q>1) {
    FThat[1:sim.list[[d]]$N,]<-FitList$That[,1:sim.list[[d]]$Q]
    dFThat[1:sim.list[[d]]$N,]<-FitList$That[,1:sim.list[[d]]$Q]-SimList$gen.theta
  } else {
    FThat[1:sim.list[[d]]$N,]<-FitList$That[,1]
    dFThat[1:sim.list[[d]]$N,]<-FitList$That[,1]-SimList$gen.theta
  }
  for (i in 2:sim.list[[d]]$Reps) {
    sc2<-readRDS(paste0(simdir,"/",SFileString(sim.list[[d]],gen=TRUE),"_",i,".rds"))
    if (sum(abs(sc2$gen.xi[,1:sim.list[[d]]$Q] - SimList$gen.xi[,1:sim.list[[d]]$Q]))>0.00001) {
      print(paste("Check replication",i,"of condition",d,
                  "cause there is a discrepancy in generated slopes."))
    }
    FitList<-readRDS(paste0(fitdir,"/",SFileString(sim.list[[d]],gen=FALSE,r = i),".rds"))
    gThat[sim.list[[d]]$N*(i-1)+1:sim.list[[d]]$N,]<-sc2$gen.theta
    if (sim.list[[d]]$Q>1) {
      That[sim.list[[d]]$N*(i-1)+1:sim.list[[d]]$N,]<-FitList$Tmap[,1:sim.list[[d]]$Q]
      dThat[sim.list[[d]]$N*(i-1)+1:sim.list[[d]]$N,]<-FitList$Tmap[,1:sim.list[[d]]$Q]-sc2$gen.theta
    } else {
      That[sim.list[[d]]$N*(i-1)+1:sim.list[[d]]$N,]<-FitList$Tmap
      dThat[sim.list[[d]]$N*(i-1)+1:sim.list[[d]]$N,]<-FitList$Tmap-sc2$gen.theta
    }
    if (sim.list[[d]]$Q>1) {
      FThat[sim.list[[d]]$N*(i-1)+1:sim.list[[d]]$N,]<-FitList$That[,1:sim.list[[d]]$Q]
      dFThat[sim.list[[d]]$N*(i-1)+1:sim.list[[d]]$N,]<-FitList$That[,1:sim.list[[d]]$Q]-SimList$gen.theta
    } else {
      FThat[sim.list[[d]]$N*(i-1)+1:sim.list[[d]]$N,]<-FitList$That[,1]
      dFThat[sim.list[[d]]$N*(i-1)+1:sim.list[[d]]$N,]<-FitList$That[,1]-sc2$gen.theta
    }
  }
  
  if (sim.list[[d]]$Q==1) {
    par(mfrow=c(1,2+ThetaFix),mar=c(5,5,3,1),...)
  } else if (sim.list[[d]]$Q==3) {
    par(mfrow=c(1,2),mar=c(5,4,3,1),...)
  } else if (sim.list[[d]]$Q==5) {
    par(mfrow=c(5,1),mar=c(5,5,3,2),...)
  } else if (sim.list[[d]]$Q==10) {
    par(mfrow=c(5,2),mar=c(5,5,3,2),...)
  } else {
    print(paste("Check condition",d,
                "cause there is a discrepancy in expectations of parameter plots."))
  }
  for (q in 1:sim.list[[d]]$Q) {
    dThat$VTiles<-cut2(gThat[[paste0("T",q)]],g=100)
    xT<-as.vector(by(That[[paste0("T",q)]],INDICES = dThat$VTiles,mean))
    yT1<-as.vector(by(dThat[[paste0("T",q)]],INDICES = dThat$VTiles,mean))
    yT2<-as.vector(by(dThat[[paste0("T",q)]],INDICES = dThat$VTiles,sd))
    yFT1<-as.vector(by(dFThat[[paste0("T",q)]],INDICES = dThat$VTiles,mean))
    yFT2<-as.vector(by(dFThat[[paste0("T",q)]],INDICES = dThat$VTiles,sd))
    ptx<-eval(parse(text=paste0("expression(italic(Theta[",q,"]))")))
    pty<-eval(parse(text=paste0("expression(hat(italic(Theta))[",q,"])")))
    #       as.integer(5000000/sim.list[[d]]$N)," burn-in iterations")
    sTitle<-eval(parse(text=paste0('expression("Condition"~"',gsub("S","",d),
                                   ' : "~hat(italic(Theta))[',q,']~": Heatmap")')))
    ContourPlot(var1 = gThat[[paste0("T",q)]],var2 = That[[paste0("T",q)]],
                xlab=ptx,ylab=pty,main=sTitle)
    pty<-eval(parse(text=paste0("expression(hat(italic(Theta))[",q,"] - italic(Theta)[",q,"])")))
    sTitle<-eval(parse(text=paste0('expression("Percentiles of"~italic(Theta)[',q,']~": EAP")')))
    plot(xT,yT1,pch=16,xlim=0.8*range(gThat[,q]),ylim=1.05*range(c(yT1-1.96*yT2,yT1+1.96*yT2)),
         main=sTitle,xlab=ptx,ylab=pty)
    points(xT,yT1,pch=16)
    arrows(xT,yT1-1.96*yT2,xT,yT1+1.96*yT2,
           code=3,angle=90,length=0.04)
    print("Percentiles")
    print(cbind(xT,yT1,yT2))
    print(cbind(xT,yFT1,yFT2))
    print(t.test(yT1[40:60],mu = 0))
    
    if (ThetaFix) {
      sTitle<-eval(parse(text=paste0('expression("Percentiles of"~italic(Theta)[',q,']~": 100 Samples")')))
      plot(xT,yFT1,pch=16,xlim=0.8*range(gThat[,q]),ylim=1.05*range(c(yFT1-1.96*yFT2,yFT1+1.96*yFT2)),
           main=sTitle,xlab=ptx,ylab=pty)
      points(xT,yFT1,pch=16)
      arrows(xT,yFT1-1.96*yFT2,xT,yFT1+1.96*yFT2,
             code=3,angle=90,length=0.04)
    }
    
  }
  
  # 
  # cat("\n")
  # for (q in 1:ncol(SimList$gen.xi)) {
  #   ss<-summary(lm(bias[,q]~SimList$gen.xi[,q]))
  #   if (q<=sim.list[[d]]$Q) {
  #     first<-paste0("$A_",q,"$ & $")
  #   } else if (q==sim.list[[d]]$Q+1) {
  #     first<-paste0("$b$ & $")
  #   } else {
  #     first<-paste0("$g$ & $")
  #   }
  #   cat(paste0(first,signif(mean(bias[,q]),digits=3),ifelse(t.test(bias[,q])$p.value<0.001,"^{***}",
  #                             ifelse(t.test(bias[,q])$p.value<0.01,"^{**}",
  #                                    ifelse(t.test(bias[,q])$p.value<0.05,"^{*}",""))),"$ & $",
  #                signif(sd(bias[,q]),digits=3),"$ & $",
  #                signif(ss$coefficients[1,1],digits=3),
  #                ifelse(ss$coefficients[1,4]<0.001,"^{***}",
  #                       ifelse(ss$coefficients[1,4]<0.01,"^{**}",
  #                              ifelse(ss$coefficients[1,4]<0.05,"^{*}",""))),"$ & $",
  #                signif(ss$coefficients[1,2],digits=3),"$ & $",
  #                signif(ss$coefficients[2,1],digits=3),
  #                ifelse(ss$coefficients[2,4]<0.001,"^{***}",
  #                       ifelse(ss$coefficients[2,4]<0.01,"^{**}",
  #                              ifelse(ss$coefficients[2,4]<0.05,"^{*}",""))),"$ & $",
  #                signif(ss$coefficients[2,2],digits=3),"$ \\\\ \\\\hline \n"))
  # }
  # if (sim.list[[d]]$K>2) {
  #   for (k in 1:ncol(SimList$gen.tau)) {
  #     ss<-summary(lm(tau.bias[,k]~SimList$gen.tau[,k]))
  # 
  #     cat(paste0("$\\tau_",k,"$ & $",signif(mean(tau.bias[,k]),digits=3),ifelse(t.test(tau.bias[,q])$p.value<0.001,"^{***}",
  #                                            ifelse(t.test(tau.bias[,k])$p.value<0.01,"^{**}",
  #                                                   ifelse(t.test(tau.bias[,k])$p.value<0.05,"^{*}",""))),"$ & $",
  #                signif(sd(tau.bias[,k]),digits=3),"$ & $",
  #                signif(ss$coefficients[1,1],digits=3),
  #                ifelse(ss$coefficients[1,4]<0.001,"^{***}",
  #                       ifelse(ss$coefficients[1,4]<0.01,"^{**}",
  #                              ifelse(ss$coefficients[1,4]<0.05,"^{*}",""))),"$ & $",
  #                signif(ss$coefficients[1,2],digits=3),"$ & $",
  #                signif(ss$coefficients[2,1],digits=3),
  #                ifelse(ss$coefficients[2,4]<0.001,"^{***}",
  #                       ifelse(ss$coefficients[2,4]<0.01,"^{**}",
  #                              ifelse(ss$coefficients[2,4]<0.05,"^{*}",""))),"$ & $",
  #                signif(ss$coefficients[2,2],digits=3),"$ \\\\ \\hline \n"))
  #   }
  # }
}

