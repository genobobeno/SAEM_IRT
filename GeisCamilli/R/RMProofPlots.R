RMProofPlots<-function(condition,items=NA,basedir="~/ParSAEM/SAEM_IRT/") {
  source(paste0(basedir,"CreateSimulationStructure.R"))
  d<-condition  
  simdir<-paste0(basedir,"/",gen.dir,"/",d)
  fitdir<-paste0(basedir,"/",fit.dir,"/",d)
  SimList<-readRDS(paste0(simdir,"/",SFileString(sim.list[[d]],gen=TRUE),"_1.rds"))
  FitList<-readRDS(paste0(fitdir,"/",SFileString(sim.list[[d]],gen=FALSE,r = 1),".rds"))
  par(mfrow=c(1,1))
  plot(c(0.4,1.7),c(-5,3),type="n",xlab="A",ylab="b",main=expression("Item Parameter Space"),xlim=c(0.3,1.55),ylim=c(-3.3,3.0))
  points(SimList$gen.xi[,1],SimList$gen.xi[,2],pch=16,col=1)
  text(SimList$gen.xi[,1]+0.03,SimList$gen.xi[,2]+0.15,1:sim.list[[d]]$J,cex=0.8)
  if (is.na(items)[1]) {
    items<-QueryUser("Pick Your Favorite Six items, separated by commas:",choices = NA,type = "character",defaultchoice = "1:6")
    items<-eval(parse(text=paste0("c(",items,")")))
  }
  points(SimList$gen.xi[items,1],SimList$gen.xi[items,2],pch=0,col=1,cex=2)
  legend("topright",legend = c("{A,b}","Chosen"),pch=c(1,0))
  # items<-c(66,37,5,11,38,63);basedir<-"./";d<-"S1"
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
      Biter$Run<-Aiter$Run<-c(0,rep(1,ceiling(bi*0.8)-2),2,
                              rep(3,ceiling(bi*0.2)-1),4,
                              5,rep(6,(nrow(Aiter)-bi-2)),7)
    } else {
      Ai<-as.data.frame(Ai)
      Bi<-as.data.frame(Bi)
      Bi$Run<-Ai$Run<-c(0,rep(1,ceiling(bi*0.8)-2),2,
                        rep(3,ceiling(bi*0.2)-1),4,
                        5,rep(6,(nrow(Ai)-bi-2)),7)
      Aiter<-rbind(Aiter,Ai)
      Biter<-rbind(Biter,Bi)
    }
  }
  if (sim.list[[d]]$Q==1) {
    par(mfrow=c(3,2))
    Aiter$chunk<-cumsum(Aiter$Run==0)
    Biter$chunk<-cumsum(Biter$Run==0)
    for (j in 1:length(items)) {
      TF<-Aiter$Run %in% c(2,4)
      sTitle<-eval(parse(text=paste0('expression("Item ',items[j],' : Annealing Drift")')))
      plot(range(Aiter[TF,paste0("J",items[j])]),range(Aiter[TF,paste0("J",items[j])]),
           type="n",xlab="A",ylab="b",main=sTitle,
           xlim=c(0.98,1.02)*range(Aiter[TF,paste0("J",items[j])]),ylim=range(Biter[TF,paste0("J",items[j])]))
      for (ch in 1:max(Aiter$chunk)) {
        points(Aiter[Aiter$Run==2 & Aiter$chunk==ch,paste0("J",items[j])],
               Biter[Biter$Run==2 & Biter$chunk==ch,paste0("J",items[j])],col=1)
        arrows(Aiter[Aiter$Run==2 & Aiter$chunk==ch,paste0("J",items[j])],
               Biter[Biter$Run==2 & Biter$chunk==ch,paste0("J",items[j])],
               Aiter[Aiter$Run==4 & Aiter$chunk==ch,paste0("J",items[j])],
               Biter[Biter$Run==4 & Biter$chunk==ch,paste0("J",items[j])],code=2,angle = 20,length = 0.08)
        points(Aiter[Aiter$Run==4 & Aiter$chunk==ch,paste0("J",items[j])],
               Biter[Biter$Run==4 & Biter$chunk==ch,paste0("J",items[j])],pch=16)
      }
      points(SimList$gen.xi[items[j],1],SimList$gen.xi[items[j],2], pch = 24, cex=2, col="blue", bg="red", lwd=2)
      print(paste("Annealing Covariance :: J =",items[j]))
      print(cov(cbind(Aiter[Aiter$Run==2,paste0("J",items[j])],Biter[Biter$Run==2,paste0("J",items[j])])))
      print(cov(cbind(Aiter[Aiter$Run==4,paste0("J",items[j])],Biter[Biter$Run==4,paste0("J",items[j])])))
      cv1<-cov(cbind(Aiter[Aiter$Run==2,paste0("J",items[j])],Biter[Biter$Run==2,paste0("J",items[j])]))
      cv2<-cov(cbind(Aiter[Aiter$Run==4,paste0("J",items[j])],Biter[Biter$Run==4,paste0("J",items[j])]))
      cat(paste(items[j],"& $\\left[ \\begin{array}{cc}",signif(cv1[1],digits=3),"&",signif(cv1[3],digits=3),
                  "\\\\",signif(cv1[2],digits=3),"&",signif(cv1[4],digits=3),"\\end{array} \\right]$ & ",
                  "$\\left[ \\begin{array}{cc}",signif(cv2[1],digits=3),"&",signif(cv2[3],digits=3),
                  "\\\\",signif(cv2[2],digits=3),"&",signif(cv2[4],digits=3),"\\end{array} \\right]$ \\\\ \\hline \n"))
    }
    for (j in 1:length(items)) {
      TF<-Aiter$Run %in% c(5,7)
      sTitle<-eval(parse(text=paste0('expression("Item ',items[j],' : RM Drift")')))
      plot(range(Aiter[TF,paste0("J",items[j])]),range(Aiter[TF,paste0("J",items[j])]),
           type="n",xlab="A",ylab="b",main=sTitle,
           xlim=c(0.98,1.02)*range(Aiter[TF,paste0("J",items[j])]),ylim=range(Biter[TF,paste0("J",items[j])]))
      for (ch in 1:max(Aiter$chunk)) {
        points(Aiter[Aiter$Run==5 & Aiter$chunk==ch,paste0("J",items[j])],
               Biter[Biter$Run==5 & Biter$chunk==ch,paste0("J",items[j])],col=1)
        arrows(Aiter[Aiter$Run==5 & Aiter$chunk==ch,paste0("J",items[j])],
               Biter[Biter$Run==5 & Biter$chunk==ch,paste0("J",items[j])],
               Aiter[Aiter$Run==7 & Aiter$chunk==ch,paste0("J",items[j])],
               Biter[Biter$Run==7 & Biter$chunk==ch,paste0("J",items[j])],code=2,angle = 30,length = 0.08)
        points(Aiter[Aiter$Run==7 & Aiter$chunk==ch,paste0("J",items[j])],
               Biter[Biter$Run==7 & Biter$chunk==ch,paste0("J",items[j])],pch=16)
      }
      points(SimList$gen.xi[items[j],1],SimList$gen.xi[items[j],2], pch = 24, cex=2, col="blue", bg="red", lwd=2)
      print(paste("RM Covariance :: J =",items[j]))
      print(cov(cbind(Aiter[Aiter$Run==5,paste0("J",items[j])],Biter[Biter$Run==5,paste0("J",items[j])])))
      print(cov(cbind(Aiter[Aiter$Run==7,paste0("J",items[j])],Biter[Biter$Run==7,paste0("J",items[j])])))
      cv1<-cov(cbind(Aiter[Aiter$Run==5,paste0("J",items[j])],Biter[Biter$Run==5,paste0("J",items[j])]))
      cv2<-cov(cbind(Aiter[Aiter$Run==7,paste0("J",items[j])],Biter[Biter$Run==7,paste0("J",items[j])]))
      cat(paste(items[j],"& $\\left[ \\begin{array}{cc}",signif(cv1[1],digits=3),"&",signif(cv1[3],digits=3),
                  "\\\\",signif(cv1[2],digits=3),"&",signif(cv1[4],digits=3),"\\end{array} \\right]$ & ",
                  "$\\left[ \\begin{array}{cc}",signif(cv2[1],digits=3),"&",signif(cv2[3],digits=3),
                  "\\\\",signif(cv2[2],digits=3),"&",signif(cv2[4],digits=3),"\\end{array} \\right]$ \\\\ \\hline \n"))
      
    }
  }
}