PlotEVRatios<-function(d,extra.dimensions=FALSE,E=10,tw=0.999,...) {
  # extra.dimensions=TRUE; d="S4";r=1
  fd<-d
  source("CreateSimulationStructure.R")
  d<-fd
  for (di in d) {
    if (extra.dimensions & di %in% c("S1","S2","S3")) {
      fit.dir.tw <- "TWFits"
      twd <- paste0(di,c("Plus1","Times5"))
      fitdir<-c(paste0(fit.dir,"/",di),paste0(fit.dir.tw,"/",twd))
    } else if (extra.dimensions & di %in% c("S4","S5","S6","S7","S8","S9")) {
      fit.dir.tw <- "TWFits"
      twd <- paste0(di,c("Minus2","Minus1","Plus1","Plus2"))
      fitdir<-c(paste0(fit.dir.tw,"/",twd[1:2]),
                paste0(fit.dir,"/",di),
                paste0(fit.dir.tw,"/",twd[3:4]))
    } else {
      fitdir<-paste0(fit.dir,"/",di)    
    }
    simdir<-paste0(gen.dir,"/",di)
    SimList<-readRDS(paste0(simdir,"/",SFileString(sim.list[[di]],gen=TRUE),"_1.rds"))
    FitList<-readRDS(paste0(fitdir[length(fitdir)],"/",SFileString(sim.list[[di]],gen=FALSE,r = 1),".rds"))
    if (di %in% c("S1","S2","S3")) {
      if (di=="S3") {
        S<-FitList$EZZ #-t(as.matrix(FitList$EZ))%*%as.matrix(FitList$EZ)
      } else {
        S<-FitList$EZZ-as.matrix(FitList$EZ)%*%t(as.matrix(FitList$EZ))
      }
      gEV<-eigen(S, symmetric=TRUE)
      cEV<-eigen(cov2cor(S),symmetric = TRUE)
    } else {
      gEV<-eigen(FitList$EZZ,symmetric=TRUE)
      cEV<-eigen(cov2cor(FitList$EZZ),symmetric = TRUE)
    }
    tw.999<-InverseTW(prob=tw,p=sim.list[[di]]$J,n=sim.list[[di]]$N)
    if (sum(d %in% c("S1","S2","S3"))>0 & di==d[1]) {
      #par(mfrow=c(length(d),2),...)
      pCH<-matrix(c(0,2,5,15,17,18),nrow=3,ncol=2)
      pBG<-c(0,3)
    } else if (sum(d %in% c("S4","S5","S6","S7","S8","S9"))>0 & di==d[1]) {
      #par(mfrow=c(length(d),1+ratios),...)
      pCH<-c(22,23,21,24,25)
      pCH.s<-c(1,1.2,1.2,1,1)
      pBG<-c(0,3)
    }
    if (length(d)==1 && sim.list[[d]]$Q<2){
      # if (d %in% c("S1","S2","S3")) {
      #   par(mfrow=c(1,1))
      # }
      FitList<-readRDS(paste0(fitdir[1],"/",SFileString(sim.list[[d]],gen=FALSE,r = 1),".rds"))
      # if (d %in% c("S1","S2","S3")) {
      if (d=="S3") {
        S<-FitList$EZZ #-t(as.matrix(FitList$EZ))%*%as.matrix(FitList$EZ)
      } else {
        S<-FitList$EZZ-as.matrix(FitList$EZ)%*%t(as.matrix(FitList$EZ))
      }
      gEV<-eigen(S, symmetric=TRUE)
      cEV<-eigen(cov2cor(S),symmetric = TRUE)
      # } else {
      #   gEV<-eigen(FitList$EZZ,symmetric=TRUE)
      #   cEV<-eigen(cov2cor(FitList$EZZ),symmetric = TRUE)
      # }
      Q<-1
      labNames <- c(paste0("Condition ",gsub('S','',di)," : Adjacent Eigenvalue Ratios"))
      sTitle<-bquote(.(labNames[1]))
      plot(Q:E,c(0.85,gEV$values[Q:(E-1)]/gEV$values[(Q+1):E]),type="n",
           main=sTitle,  #expression("Condition"~.(din)~": Adjacent Ratios of Eigenvalues"),
           xlab="number of factors",ylab="estimators",lty=1,pch=1)
      for (i in 1:length(fitdir)) {
        #for (r in 1:sim.list[[d]]$Reps) {
        FitList<-readRDS(paste0(fitdir[i],"/",SFileString(sim.list[[d]],gen=FALSE,r = 1),".rds"))
        if (d %in% c("S1","S2","S3")) {
          if (d=="S3") {
            S<-FitList$EZZ #-t(as.matrix(FitList$EZ))%*%as.matrix(FitList$EZ)
          } else {
            S<-FitList$EZZ-as.matrix(FitList$EZ)%*%t(as.matrix(FitList$EZ))
          }
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
    } else {
      FitList<-readRDS(paste0(fitdir[3],"/",SFileString(sim.list[[di]],gen=FALSE,r = 1),".rds"))
      gEV<-eigen(FitList$EZZ,symmetric=TRUE)
      cEV<-eigen(cov2cor(FitList$EZZ),symmetric = TRUE)
      Q<-1
      labNames <- c(paste0("Condition ",gsub('S','',di)," : Adjacent Eigenvalue Ratios"))
      sTitle<-bquote(.(labNames[1]))
      if (di %in% c("S4")) {
        plot(Q:E,c(0.85,gEV$values[Q:(E-1)]/gEV$values[(Q+1):E]),type="n",
             main=sTitle,#expression("Condition"~.(din)~": Adjacent Ratios of Eigenvalues"),
             xlab="number of factors",ylab="estimators",lty=1,pch=1,ylim=c(0.5,11.5))
      } else {
        plot(Q:E,c(0.85,gEV$values[Q:(E-1)]/gEV$values[(Q+1):E]),type="n",
             main=sTitle,#expression("Condition"~.(din)~": Adjacent Ratios of Eigenvalues"),
             xlab="number of factors",ylab="estimators",lty=1,pch=1)
      }
      for (i in 1:length(fitdir)) {
        #for (r in 1:sim.list[[d]]$Reps) {
        FitList<-readRDS(paste0(fitdir[i],"/",SFileString(sim.list[[di]],gen=FALSE,r = 1),".rds"))
        gEV<-eigen(FitList$EZZ,symmetric=TRUE)
        cEV<-eigen(cov2cor(FitList$EZZ),symmetric = TRUE)
        
        print(fitdir[i])
        print(sapply(1:E,function(x) t.test(gEV$vectors[,x],mu = 0)$p.value))
        gTW.TF<-TWTransform(gEV$values,p = sim.list[[di]]$J,n = sim.list[[di]]$N,ptest = T)
        cTW.TF<-TWTransform(cEV$values,p = sim.list[[di]]$J,n = sim.list[[di]]$N,ptest = T)
        #j.r<-sum(TW.TF>0.95)+1
        gER<-gEV$values[1:(sim.list[[di]]$J-1)]/gEV$values[2:sim.list[[di]]$J]
        cER<-cEV$values[1:(sim.list[[di]]$J-1)]/cEV$values[2:sim.list[[di]]$J]
        lines(Q:E,gER[Q:E],lty=1)
        points(Q:E,gER[Q:E],pch=pCH[i],bg=ifelse(gEV$values[Q:E]>tw.999,pBG[2],pBG[1]),col=1,cex=pCH.s[i])
        lines(Q:E,cER[Q:E],lty=4)
        points(Q:E,cER[Q:E],pch=pCH[i],bg=ifelse(cEV$values[Q:E]>tw.999,pBG[2],pBG[1]),col=1,cex=pCH.s[i])
      }
      if (di %in% c("S4","S5")) {
        ltext =  c(expression(lambda["Q=1"]("S"[2])),expression(lambda["Q=1"]("cor"("S"[2]))),
                   expression(lambda["Q=2"]("S"[2])),expression(lambda["Q=2"]("cor"("S"[2]))),
                   expression(lambda["Q=3"]("S"[2])),expression(lambda["Q=3"]("cor"("S"[2]))),
                   expression(lambda["Q=4"]("S"[2])),expression(lambda["Q=4"]("cor"("S"[2]))),
                   expression(lambda["Q=5"]("S"[2])),expression(lambda["Q=5"]("cor"("S"[2]))))
        strwidth(ltext)
        legend("topright",legend = ltext,text.width = max(strwidth(ltext))*1.02,
               lty=c(1,4,1,4,1,4,1,4,1,4),pch=c(rep(pCH,c(2,2,2,2,2))),pt.cex=c(rep(pCH.s,c(2,2,2,2,2))))
        
      } else if (di %in% c("S6","S7")) {
        ltext =  c(expression(lambda["Q=3"]("S"[2])),expression(lambda["Q=3"]("cor"("S"[2]))),
                   expression(lambda["Q=4"]("S"[2])),expression(lambda["Q=4"]("cor"("S"[2]))),
                   expression(lambda["Q=5"]("S"[2])),expression(lambda["Q=5"]("cor"("S"[2]))),
                   expression(lambda["Q=6"]("S"[2])),expression(lambda["Q=6"]("cor"("S"[2]))),
                   expression(lambda["Q=7"]("S"[2])),expression(lambda["Q=7"]("cor"("S"[2]))))
        strwidth(ltext)
        legend("topright",legend = ltext,text.width = max(strwidth(ltext))*1.02,
               lty=c(1,4,1,4,1,4,1,4,1,4),pch=c(rep(pCH,c(2,2,2,2,2))),pt.cex=c(rep(pCH.s,c(2,2,2,2,2))))
        
      } else if (di %in% c("S8","S9")) {
        ltext =  c(expression(lambda["Q=8"]("S"[2])),expression(lambda["Q=8"]("cor"("S"[2]))),
                   expression(lambda["Q=9"]("S"[2])),expression(lambda["Q=9"]("cor"("S"[2]))),
                   expression(lambda["Q=10"]("S"[2])),expression(lambda["Q=10"]("cor"("S"[2]))),
                   expression(lambda["Q=11"]("S"[2])),expression(lambda["Q=11"]("cor"("S"[2]))),
                   expression(lambda["Q=12"]("S"[2])),expression(lambda["Q=12"]("Cor"("S"[2]))))
        strwidth(ltext)
        legend("topright",legend = ltext,text.width = max(strwidth(ltext))*1.02,
               lty=c(1,4,1,4,1,4,1,4,1,4),pch=c(rep(pCH,c(2,2,2,2,2))),pt.cex=c(rep(pCH.s,c(2,2,2,2,2))))
      }
    }
  }
}
