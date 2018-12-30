QuickDiagnostic<-function(simCondition=1,replication=1) {

  gen.dir<-"GeneratedFiles"
  fit.dir<-"ConvergedModelFits"
  sim.list<-readRDS(paste0(gen.dir,"/Simulations.rds"))  
    
  d <- paste0("S",simCondition)
  r <- replication
  
  simdir<-paste0(gen.dir,"/",d)
  SimList<-readRDS(paste0(simdir,"/",SFileString(sim.list[[d]],gen=TRUE),"_",r,".rds"))
  print(SimList$gen.list$gen.xi)
  
  fitdir<-paste0(fit.dir,"/",d)
  FitList<-readRDS(paste0(fitdir,"/",SFileString(sim.list[[d]],gen=FALSE,r = r),".rds"))
  
  if (d=="S2") {
    par(mfrow=c(2,3))
    plot(FitList$XI[,2],FitList$xi[,2],pch=16,col=4,main="blue=B, red=A, green=C")
    points(FitList$XI[,1],FitList$xi[,1],pch=16,col=2)
    points(FitList$XI[,3],FitList$xi[,3],pch=16,col=3)
    plot(FitList$XI[,3],FitList$xi[,3],pch=16,col=3,main="Guessing with linear fit")
    Clm<-lm(FitList$xi[,3]~FitList$XI[,3])
    abline(Clm)
    text(0.15,0.15,paste0(paste0(coef(Clm),collapse=" + ")," * x"))
    plot(density(FitList$xi[,1] - FitList$XI[,1]),main="A",col=2)
    plot(density(FitList$xi[,2] - FitList$XI[,2]),main="B",col=4)
    plot(density(FitList$xi[,3] - FitList$XI[,3]),main="C",col=3)
    plot(density(FitList$THETA-FitList$That[,1]),main="Theta")
  } else {
    AA<-ncol(FitList$XI)
    KK<-ifelse(with(FitList,exists("TAU")),ncol(FitList$TAU),0)
    #  sim.list$S3
    par(mfrow=c(5,ceiling((AA+KK+2)/5)),mar=c(3,2,0.5,0.5))
    for (i in 1:(AA-1)) {
      Ah<-sign(sum(FitList$xi[,i]))*FitList$xi[,i]
      plot(FitList$XI[,i],Ah,pch=16,col=1+i,main=paste0("Slope ",i))
      Clm<-lm(Ah~FitList$XI[,i])
      abline(Clm)
      text(0.2+(i==1)*0.8,0.6,paste0(paste0(format(coef(Clm),digits = 3),collapse="+"),"*x"))
    }
    Ah<-sign(sum(FitList$xi[,1]))*FitList$xi[,1]
    plot(density(Ah - FitList$XI[,1]),main="A",col=2)
    if (AA>2) {
      for (i in 2:(AA-1)) {
        Ah<-sign(sum(FitList$xi[,i]))*FitList$xi[,i]
        lines(density(Ah - FitList$XI[,i]),col=1+i) 
      }
    }
    plot(FitList$XI[,AA],FitList$xi[,AA],pch=16,col=4,main="Intercepts")
    Clm<-lm(FitList$xi[,AA]~FitList$XI[,AA])
    abline(Clm)
    text(0.5,-0.5,paste0(paste0(format(coef(Clm),digits = 3),collapse="+"),"*x"))
    plot(density(FitList$xi[,AA] - FitList$XI[,AA]),main="Intercepts",col=4)
    if (with(FitList,exists("TAU"))) {
      for (i in 1:KK) {
        plot(FitList$TAU[,i],FitList$tau[,i],pch=16,col=1+i,main=paste0("Offset ",i))
        Clm<-lm(FitList$tau[,i]~FitList$TAU[,i])
        abline(Clm)
        text(0+(i-1)*0.5,-1+(i-1)*0.5,paste0(paste0(format(coef(Clm),digits = 3),collapse="+"),"*x"))
      }
      plot(density(FitList$tau[,1] - FitList$TAU[,1]),main="Offset",col=2)
      for (i in 2:KK) {
        lines(density(FitList$tau[,i] - FitList$TAU[,i]),col=1+i) 
      }
    }
  }
  
  delta<-FitList$XI[,1] - FitList$xi[,1]
  tests<-cbind(FitList$XI[,1],rankY(FitList$XI[,1]),delta,rankY(-abs(delta)))
  colnames(tests)<-c("A1","A1rank","dA1","dA1rank")
  if (ncol(FitList$XI)>2 && d!="S2") {
    for (i in 2:(ncol(FitList$XI)-1)) {
      delta<-FitList$XI[,i] - FitList$xi[,i]
      tests<-cbind(tests,FitList$XI[,i],rankY(FitList$XI[,i]),delta,rankY(-abs(delta)))
      colnames(tests)[ncol(tests)-3:0]<-sprintf(fmt = c("A%d","A%drank","dA%d","dA%drank"),i)
    } 
  }
  if (d!="S2") {
    delta<-FitList$XI[,ncol(FitList$XI)] - FitList$xi[,ncol(FitList$XI)]
    tests<-cbind(tests,FitList$XI[,ncol(FitList$XI)],rankY(FitList$XI[,ncol(FitList$XI)]),delta,rankY(-abs(delta)))
    colnames(tests)[ncol(tests)-3:0]<-sprintf(fmt = c("B%d","B%drank","dB%d","dB%drank"),i)
  } else {
    delta<-FitList$XI[,ncol(FitList$XI)-1] - FitList$xi[,ncol(FitList$XI)-1]
    tests<-cbind(tests,FitList$XI[,ncol(FitList$XI)-1],rankY(FitList$XI[,ncol(FitList$XI)-1]),delta,rankY(-abs(delta)))
    colnames(tests)[ncol(tests)-3:0]<-sprintf(fmt = c("B%d","B%drank","dB%d","dB%drank"),i)
    delta<-FitList$XI[,ncol(FitList$XI)] - FitList$xi[,ncol(FitList$XI)]
    tests<-cbind(tests,FitList$XI[,ncol(FitList$XI)],rankY(FitList$XI[,ncol(FitList$XI)]),delta,rankY(-abs(delta)))
    colnames(tests)[ncol(tests)-3:0]<-sprintf(fmt = c("C%d","C%drank","dC%d","dC%drank"),i)
    
  }
  print(tests)

}