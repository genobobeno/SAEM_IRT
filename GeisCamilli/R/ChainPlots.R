ChainPlots<-function(condition,repl=1,items=1:10,basedir="~/ParSAEM/SAEM_IRT/") {
  #condition<-"S1";repl=NA;items = c(53,37,5,11,38,63);basedir = "./"
  source(paste0(basedir,"CreateSimulationStructure.R"))
  d<-condition  
  simdir<-paste0(basedir,"/",gen.dir,"/",d)
  fitdir<-paste0(basedir,"/",fit.dir,"/",d)

  if (!is.na(repl)) {
    SimList<-readRDS(paste0(simdir,"/",SFileString(sim.list[[d]],gen=TRUE),"_",repl,".rds"))
    FitList<-readRDS(paste0(fitdir,"/",SFileString(sim.list[[d]],gen=FALSE,r = repl),".rds"))
    load(file=paste0(fitdir,"/",SFileString(sim.list[[d]],gen=FALSE,r = repl),".rda"))    
    bi<-FitList$settings$burnin
    if (sim.list[[d]]$Q==1) {
      par(mfrow=c(1,2),mar=c(5,4,4,2))
    } else {
      par(mfrow=c(1,ncol(SimList$gen.xi)),mar=c(5,4,4,2))
    }
    its<-dim(MCMCDATA$Aiter)[3]
    lcols<-rainbow(length(items),start=0.3,0.75)
    for (q in 1:(ncol(SimList$gen.xi)-1)) {
      sTitle<-eval(parse(text=paste0('expression("Condition"~"',gsub("S","",d),
                                     ' : "~italic(A[',q,'])~": MCMC Chain")')))
      plot(MCMCDATA$Aiter[1,q,],type="n",main=sTitle,
           xlab="Iterations",ylab="Slope",ylim=c(-0.1,0.1)+range(SimList$gen.xi[,q]))
      for (j in 1:length(items)) {
        lines(MCMCDATA$Aiter[items[j],q,],col=lcols[j])
        abline(h=SimList$gen.xi[items[j],q],col=lcols[j],lwd=1,lty=2)
        if (sum(MCMCDATA$Aiter[items[-j],q,its]-MCMCDATA$Aiter[items[j],q,its]<0.1 &
            MCMCDATA$Aiter[items[-j],q,its]-MCMCDATA$Aiter[items[j],q,its]>0)>0) {
          text(its-ceiling(its/12),MCMCDATA$Aiter[items[j],q,its]-0.03,paste0("j=",items[j]),cex=0.8)
        } else {          
          text(its-ceiling(its/12),MCMCDATA$Aiter[items[j],q,its]+0.03,paste0("j=",items[j]),cex=0.8)
        }
      }
      points(rep(1,length(items)),MCMCDATA$Aiter[items,q,1],pch=16)
      abline(v=c(ceiling(bi*0.8),bi),lty=2)
    }
    sTitle<-eval(parse(text=paste0('expression("Condition"~"',gsub("S","",d),
                                   ' : "~italic(b)~": MCMC Chain")')))
    plot(MCMCDATA$Biter[1,],type="n",main=sTitle,xlab="Iterations",ylab="Intercept",
         ylim=1.05*range(MCMCDATA$Biter[,1]))
    for (j in 1:length(items)) {
      lines(MCMCDATA$Biter[items[j],],col=lcols[j])
      abline(h=SimList$gen.xi[items[j],sim.list[[d]]$Q+1],col=lcols[j],lwd=1,lty=2)
      if (sum(MCMCDATA$Biter[items[-j],its]-MCMCDATA$Biter[items[j],its]<0.3 &
                  MCMCDATA$Biter[items[-j],its]-MCMCDATA$Biter[items[j],its]>0)>0) {
        text(its-ceiling(its/12),MCMCDATA$Biter[items[j],its]-0.15,paste0("j=",items[j]),cex=0.8)
      } else {          
        text(its-ceiling(its/12),MCMCDATA$Biter[items[j],its]+0.15,paste0("j=",items[j]),cex=0.8)
      }
    }
    points(rep(1,length(items)),MCMCDATA$Biter[items,1],pch=16)
    abline(v=c(ceiling(bi*0.8),bi),lty=2)
    if (sim.list[[d]]$Q==1) {
      par(mfrow=c(2,2))
      sTitle<-eval(parse(text=paste0('expression("MCMC Chain Paths :"~italic(t)~"= {1:',ceiling(bi*0.8),'}")')))
      xiMC<-c(1,0)
      plot(c(0.2,1.7),c(-5,5),type="n",xlab="A",ylab="b",main=sTitle,xlim=c(0.2,1.7),ylim=c(-5,5))
      for (j in 1:length(items)) {
        xiMCi<-cbind(MCMCDATA$Aiter[items[j],1,1:ceiling(bi*0.8)],MCMCDATA$Biter[items[j],1:ceiling(bi*0.8)])
        xiMC<-rbind(xiMC,xiMCi)
        lines(xiMCi[,1],xiMCi[,2],col=lcols[j])
        points(xiMCi[1,1],xiMCi[1,2])
        points(SimList$gen.xi[items[j],1],SimList$gen.xi[items[j],2],pch=24,bg="lightblue",col=1,lwd=2,cex=1.2)
        text(SimList$gen.xi[items[j],1]-0.04,SimList$gen.xi[items[j],2]-0.3,paste0("j=",items[j]),cex=0.8)
      }
      sTitle<-eval(parse(text=paste0('expression("MCMC Chain Heatmap :"~italic(t)~"= {1:',ceiling(bi*0.8),'}")')))
      ContourPlot(var1 = xiMC[,1],var2=xiMC[,2],xlab="A",ylab="b",main=sTitle,xlim=c(0.2,1.7),ylim=c(-5,5))
      points(SimList$gen.xi[items,1],SimList$gen.xi[items,2],pch=24,bg="lightblue",col=1,lwd=2,cex=1.2)
      text(SimList$gen.xi[items,1]-0.04,SimList$gen.xi[items,2]-0.3,paste0("j=",items),cex=0.8)
      
      sTitle<-eval(parse(text=paste0('expression("MCMC Chain Paths :"~italic(t)~"= {',ceiling(bi*0.8),':',bi,'}")')))
      xiMC<-c(1,0)
      plot(c(0.2,1.7),c(-5,5),type="n",xlab="A",ylab="b",main=sTitle,xlim=c(0.2,1.7),ylim=c(-5,5))
      for (j in 1:length(items)) {
        xiMCi<-cbind(MCMCDATA$Aiter[items[j],1,ceiling(bi*0.8):bi],MCMCDATA$Biter[items[j],ceiling(bi*0.8):bi])
        xiMC<-rbind(xiMC,xiMCi)
        lines(xiMCi[,1],xiMCi[,2],col=lcols[j])
        points(SimList$gen.xi[items[j],1],SimList$gen.xi[items[j],2],pch=24,bg="lightblue",col=1,lwd=2,cex=1.2)
        text(SimList$gen.xi[items[j],1]-0.04,SimList$gen.xi[items[j],2]-0.3,paste0("j=",items[j]),cex=0.8)
      }
      sTitle<-eval(parse(text=paste0('expression("MCMC Chain Heatmap :"~italic(t)~"= {',ceiling(bi*0.8),':',bi,'}")')))
      ContourPlot(var1 = xiMC[,1],var2=xiMC[,2],xlab="A",ylab="b",main=sTitle,xlim=c(0.25,1.7),ylim=c(-5,4))
      points(SimList$gen.xi[items,1],SimList$gen.xi[items,2],pch=24,bg="lightblue",col=1,lwd=2,cex=1.2)
      text(SimList$gen.xi[items,1]-0.04,SimList$gen.xi[items,2]-0.3,paste0("j=",items),cex=0.8)
      
    }
  } else {
    SimList<-readRDS(paste0(simdir,"/",SFileString(sim.list[[d]],gen=TRUE),"_1.rds"))
    FitList<-readRDS(paste0(fitdir,"/",SFileString(sim.list[[d]],gen=FALSE,r = 1),".rds"))
    bi<-FitList$settings$burnin
    Aiter<-list()
    Biter<-list()
    for (r in 1:sim.list[[d]]$Reps) {
      load(file=paste0(fitdir,"/",SFileString(sim.list[[d]],gen=FALSE,r = r),".rda"))
      Ai<-list()
      Bi<-list()
      for (j in 1:length(items)) {
        if (r==1) {
          Aiter[[paste0("J",items[j])]]<-MCMCDATA$Aiter[items[j],1,]
          Biter[[paste0("J",items[j])]]<-MCMCDATA$Biter[items[j],]
        } else {
          Ai[[paste0("J",items[j])]]<-MCMCDATA$Aiter[items[j],1,]
          Bi[[paste0("J",items[j])]]<-MCMCDATA$Biter[items[j],]
        }
      }
      if (r==1) {
        Aiter<-as.data.frame(Aiter)
        Biter<-as.data.frame(Biter)
        Biter$Run<-Aiter$Run<-c(0,rep(1,ceiling(bi*0.8)-1),
                                rep(2,ceiling(bi*0.2)),
                                3,rep(4,(nrow(Aiter)-bi-1)))
      } else {
        Ai<-as.data.frame(Ai)
        Bi<-as.data.frame(Bi)
        Bi$Run<-Ai$Run<-c(0,rep(1,ceiling(bi*0.8)-1),
                          rep(2,ceiling(bi*0.2)),
                          3,rep(4,(nrow(Ai)-bi-1)))
        Aiter<-rbind(Aiter,Ai)
        Biter<-rbind(Biter,Bi)
      }
    }
    if (sim.list[[d]]$Q==1) {
      par(mfrow=c(2,2))
      lcols<-rainbow(length(items),start=0.3,0.75)
      sTitle<-eval(parse(text=paste0('expression("MCMC Chain Paths :"~italic(t)~"= {1:',ceiling(bi*0.8),'}")')))
      plot(c(0.2,1.7),c(-5,5),type="n",xlab="A",ylab="b",main=sTitle,xlim=c(0.2,1.7),ylim=c(-5,5))
      Aiter$chunk<-cumsum(Aiter$Run==0)
      Biter$chunk<-cumsum(Biter$Run==0)
      for (j in 1:length(items)) {
        linecol<-ifelse(Aiter$Run[Aiter$Run<2]==0,"#FFFFFFFF",lcols[j])
        for (ch in 1:max(Aiter$chunk)) {
          lines(Aiter[Aiter$Run<2 & Aiter$chunk==ch,paste0("J",items[j])],
                Biter[Biter$Run<2 & Biter$chunk==ch,paste0("J",items[j])],
                col=lcols[j])
        }
        points(Aiter[Aiter$Run==0,paste0("J",items[j])],Biter[Biter$Run==0,paste0("J",items[j])],col=1)
        points(SimList$gen.xi[items[j],1],SimList$gen.xi[items[j],2],pch=24,bg="lightblue",col=1,lwd=2,cex=1.2)
        text(SimList$gen.xi[items[j],1]-0.04,SimList$gen.xi[items[j],2]-0.3,paste0("j=",items[j]),cex=0.8)
      }
      sTitle<-eval(parse(text=paste0('expression("MCMC Chain Heatmap :"~italic(t)~"= {1:',ceiling(bi*0.8),'}")')))
      ContourPlot(var1 = unlist(Aiter[Aiter$Run<2,grepl("J",colnames(Aiter))]),
                  var2=unlist(Biter[Biter$Run<2,grepl("J",colnames(Biter))]),
                  xlab="A",ylab="b",main=sTitle,xlim=c(0.2,1.7),ylim=c(-5,5))
      points(SimList$gen.xi[items,1],SimList$gen.xi[items,2],pch=24,bg="lightblue",col=1,lwd=2,cex=1.2)
      text(SimList$gen.xi[items,1]-0.04,SimList$gen.xi[items,2]-0.3,paste0("j=",items),cex=0.8)
      
      sTitle<-eval(parse(text=paste0('expression("MCMC Chain Paths :"~italic(t)~"= {',ceiling(bi*0.8),':',bi,'}")')))
      plot(c(0.2,1.7),c(-5,5),type="n",xlab="A",ylab="b",main=sTitle,xlim=c(0.2,1.7),ylim=c(-5,5))
      for (j in 1:length(items)) {
        lines(Aiter[Aiter$Run>1&Aiter$Run<4,paste0("J",items[j])],Biter[Biter$Run>1&Biter$Run<4,paste0("J",items[j])],
              col=ifelse(Biter$Run==3,"#FFFFFFFF",lcols[j]))
        points(SimList$gen.xi[items[j],1],SimList$gen.xi[items[j],2],pch=24,bg="lightblue",col=1,lwd=2,cex=1.2)
        text(SimList$gen.xi[items[j],1]-0.04,SimList$gen.xi[items[j],2]-0.3,paste0("j=",items[j]),cex=0.8)
      }
      sTitle<-eval(parse(text=paste0('expression("MCMC Chain Heatmap :"~italic(t)~"= {',ceiling(bi*0.8),':',bi,'}")')))
      ContourPlot(var1 = unlist(Aiter[Aiter$Run>1&Aiter$Run<4,grepl("J",colnames(Aiter))]),
                  var2=unlist(Biter[Biter$Run>1&Biter$Run<4,grepl("J",colnames(Biter))]),
                  xlab="A",ylab="b",main=sTitle,xlim=c(0.2,1.7),ylim=c(-5,5))
      points(SimList$gen.xi[items,1],SimList$gen.xi[items,2],pch=24,bg="lightblue",col=1,lwd=2,cex=1.2)
      text(SimList$gen.xi[items,1]-0.04,SimList$gen.xi[items,2]-0.3,paste0("j=",items),cex=0.8)
      
    }
  }
}