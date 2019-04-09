


BenchmarkPlots<-function(condition,repl=NA,all.reps=T,basedir="~/ParSAEM/SAEM_IRT/") {
  #condition="S1";repl=NA;all.reps=T;basedir="~/ParSAEM/SAEM_IRT/"
  source(paste0(basedir,"CreateSimulationStructure.R"))
  d<-condition  
  simdir<-paste0(basedir,"/",gen.dir,"/",d)
  fitdir<-paste0(basedir,"/",fit.dir,"/",d)
  if (!all.reps && !is.na(repl)) {
    SimList<-readRDS(paste0(simdir,"/",SFileString(sim.list[[d]],gen=TRUE),"_",r,".rds"))
    FitList<-readRDS(paste0(fitdir,"/",SFileString(sim.list[[d]],gen=FALSE,r = r),".rds"))
    par(mfrow=c(1,2))
    plot(SimList$gen.xi[,1:sim.list[[d]]$Q],FitList$xi[,1:sim.list[[d]]$Q],
         main=paste("Condition",d,": Replication",repl,": Slopes"),xlab=expression(italic(A)),
         ylab=expression(hat(italic(A))))
    plot(SimList$gen.xi[,sim.list[[d]]$Q+1],FitList$xi[,sim.list[[d]]$Q+1],
         main=paste("Condition",d,": Replication",repl,": Intercepts"),xlab=expression(italic(b)),
         ylab=expression(hat(italic(b))))
  } else {
    XIhat<-array(0, dim=c(sim.list[[d]]$J,1+sim.list[[d]]$Q+sim.list[[d]]$Guessing,sim.list[[d]]$Reps))
    dXIhat<-array(0, dim=c(sim.list[[d]]$J,1+sim.list[[d]]$Q+sim.list[[d]]$Guessing,sim.list[[d]]$Reps))
    That<-array(0, dim=c(sim.list[[d]]$N,sim.list[[d]]$Q,sim.list[[d]]$Reps))
    dThat<-array(0, dim=c(sim.list[[d]]$N,sim.list[[d]]$Q,sim.list[[d]]$Reps))
    SimList<-readRDS(paste0(simdir,"/",SFileString(sim.list[[d]],gen=TRUE),"_1.rds"))
    FitList<-readRDS(paste0(fitdir,"/",SFileString(sim.list[[d]],gen=FALSE,r = 1),".rds"))
    if (sim.list[[d]]$Q>1) {
      vA<-colSums(FitList$xi[,1:sim.list[[d]]$Q])
    } else {
      vA<-sum(FitList$xi[,1])
    }
    XIhat[,,1]<-FitList$xi*matrix(c(sign(vA),rep(1,1+sim.list[[d]]$Guessing)),nrow = sim.list[[d]]$J,
                                  ncol=sim.list[[d]]$Q+1+sim.list[[d]]$Guessing,byrow=T)
    dXIhat[,,1]<-FitList$xi*matrix(c(sign(vA),
                                     rep(1,1+sim.list[[d]]$Guessing)),nrow = sim.list[[d]]$J,
                                   ncol=sim.list[[d]]$Q+1+sim.list[[d]]$Guessing,byrow=T)-SimList$gen.xi
    That[,,1]<-FitList$That[,1:sim.list[[d]]$Q]
    dThat[,,1]<-FitList$That[,1:sim.list[[d]]$Q]-SimList$gen.theta
    if (!is.na(sim.list[[d]]$K) & sim.list[[d]]$K>2) {
      Tauhat<-array(0, dim=c(sim.list[[d]]$J,sim.list[[d]]$K-1,sim.list[[d]]$Reps))
      dTauhat<-array(0, dim=c(sim.list[[d]]$J,sim.list[[d]]$K-1,sim.list[[d]]$Reps))
      Tauhat[,,1]<-FitList$tau
      dTauhat[,,1]<-FitList$tau-SimList$gen.tau
    }
    for (i in 2:sim.list[[d]]$Reps) {
      sc2<-readRDS(paste0(simdir,"/",SFileString(sim.list[[d]],gen=TRUE),"_",i,".rds"))
      if (sum(abs(sc2$gen.xi[,1:sim.list[[d]]$Q] - SimList$gen.xi[,1:sim.list[[d]]$Q]))>0.00001) {
        print(paste("Check replication",i,"of condition",d,
                    "cause there is a discrepancy in generated slopes."))
      }
      FitList<-readRDS(paste0(fitdir,"/",SFileString(sim.list[[d]],gen=FALSE,r = i),".rds"))
      if (sim.list[[d]]$Q>1) {
        vA<-colSums(FitList$xi[,1:sim.list[[d]]$Q])
      } else {
        vA<-sum(FitList$xi[,1])
      }
      XIhat[,,i]<-FitList$xi*matrix(c(sign(vA),
                                      rep(1,1+sim.list[[d]]$Guessing)),nrow = sim.list[[d]]$J,
                                    ncol=sim.list[[d]]$Q+1+sim.list[[d]]$Guessing,byrow=T)
      dXIhat[,,i]<-FitList$xi*matrix(c(sign(vA),
                                       rep(1,1+sim.list[[d]]$Guessing)),nrow = sim.list[[d]]$J,
                                     ncol=sim.list[[d]]$Q+1+sim.list[[d]]$Guessing,byrow=T)-SimList$gen.xi
      if (!is.na(sim.list[[d]]$K) & sim.list[[d]]$K>2) {
        Tauhat[,,i]<-FitList$tau
        dTauhat[,,i]<-FitList$tau-SimList$gen.tau
      }
    }
    mXI<-apply(XIhat,c(1,2),mean)
    bias<-mXI-SimList$gen.xi
    sXI<-apply(XIhat-array(rep(mXI,sim.list[[d]]$Reps),
                           dim=c(sim.list[[d]]$J,1+sim.list[[d]]$Q+sim.list[[d]]$Guessing,
                                 sim.list[[d]]$Reps)),
               c(1,2),sd)
    if (sim.list[[d]]$Q==1 & !sim.list[[d]]$Guessing & 
        (is.na(sim.list[[d]]$K) | sim.list[[d]]$K!=4)) {
      par(mfrow=c(2,1),mar=c(5,5,3,2))
    } else if (sim.list[[d]]$Q==1 & sim.list[[d]]$Guessing) {
      par(mfrow=c(3,1),mar=c(5,5,3,2))
    } else if (sim.list[[d]]$Q==1 & sim.list[[d]]$K==4) {
      par(mfrow=c(5,1),mar=c(5,5,3,2))
    } else if (sim.list[[d]]$Q==3) {
      par(mfrow=c(3,1),mar=c(5,5,3,2))
    } else if (sim.list[[d]]$Q==5) {
      par(mfrow=c(5,1),mar=c(5,5,3,2))
    } else if (sim.list[[d]]$Q==10) {
      par(mfrow=c(5,2),mar=c(5,5,3,2))
    } else {
      print(paste("Check condition",d,
                  "cause there is a discrepancy in expectations of parameter plots."))
    }
    for (q in 1:sim.list[[d]]$Q) {
      pty<-eval(parse(text=paste0("expression(hat(italic(A))[",q,"] - italic(A)[",q,"])")))
      ptx<-eval(parse(text=paste0("expression(hat(italic(A[",q,"])))")))
      sTitle<-eval(parse(text=paste0('expression("Condition"~"',gsub("S","",d),
                                     ' : "~italic(A)[',q,']~": Bias & Error")')))
      #       as.integer(5000000/sim.list[[d]]$N)," burn-in iterations")
      plot(SimList$gen.xi[,q],bias[,q],pch=19,cex=1.5,xlim=c(-0.1,1.05*max(SimList$gen.xi[,q])),
           ylim=1.05*range(c(bias[,q]-2*sXI[,q],bias[,q]+2*sXI[,q])),
           main=sTitle,xlab=ptx,ylab=pty)
      arrows(SimList$gen.xi[,q],bias[,q]-2*sXI[,q],SimList$gen.xi[,q],bias[,q]+2*sXI[,q],
             code=3,angle=90,length=0.07)
      abline(lm(bias[,q]~SimList$gen.xi[,q]),col=2,lty=2)
      print(summary(lm(bias[,q]~SimList$gen.xi[,q])))
    }
    if (!is.na(sim.list[[d]]$K) & sim.list[[d]]$K>2 & sim.list[[d]]$Q>1) {par(mfrow=c(4,1),mar=c(5,5,3,2))}
    pty<-expression(hat(italic(b)) - italic(b))
    ptx<-expression(hat(italic(b)))
    sTitle<-eval(parse(text=paste0('expression("Condition"~"',gsub("S","",d),
                                   ' : "~italic(b)~": Bias & Error")')))
    q<-1+sim.list[[d]]$Q
    plot(SimList$gen.xi[,q],bias[,q],pch=19,cex=1.5,xlim=c(-1.03,1.03)*max(abs(SimList$gen.xi[,q])),
         ylim=1.05*range(c(bias[,q]-2*sXI[,q],bias[,q]+2*sXI[,q])),
         main=sTitle,xlab=ptx,ylab=pty)
    arrows(SimList$gen.xi[,q],bias[,q]-2*sXI[,q],SimList$gen.xi[,q],bias[,q]+2*sXI[,q],code=3,angle=90,length=0.07)
    abline(lm(bias[,q]~SimList$gen.xi[,q]),col=2,lty=2)
    print(summary(lm(bias[,q]~SimList$gen.xi[,q])))
    if (!is.na(sim.list[[d]]$K) & sim.list[[d]]$K>2) {
      mTau<-apply(Tauhat,c(1,2),mean)
      tau.bias<-mTau-SimList$gen.tau
      sTau<-apply(Tauhat-array(rep(mTau,sim.list[[d]]$Reps),
                             dim=c(sim.list[[d]]$J,sim.list[[d]]$K-1,sim.list[[d]]$Reps)),
                 c(1,2),sd)
      for (k in 1:(sim.list[[d]]$K-1)) {
        pty<-eval(parse(text=paste0("expression(hat(italic(tau))[",k,"] - italic(tau)[",k,"])")))
        ptx<-eval(parse(text=paste0("expression(hat(italic(tau[",k,"])))")))
        sTitle<-eval(parse(text=paste0('expression("Condition"~"',gsub("S","",d),
                                       ' : "~italic(tau)[',k,']~": Bias & Error")')))
        #       as.integer(5000000/sim.list[[d]]$N)," burn-in iterations")
        plot(SimList$gen.tau[,k],tau.bias[,k],pch=19,cex=1.5,xlim=1.05*range(SimList$gen.tau[,k]),
             ylim=1.05*range(c(tau.bias-2*sTau,tau.bias+2*sTau)),
             main=sTitle,xlab=ptx,ylab=pty)
        arrows(SimList$gen.tau[,k],tau.bias[,k]-2*sTau[,k],SimList$gen.tau[,k],tau.bias[,k]+2*sTau[,k],code=3,angle=90,length=0.07)
        abline(lm(tau.bias[,k]~SimList$gen.tau[,k]),col=2,lty=2)
        print(summary(lm(tau.bias[,k]~SimList$gen.tau[,k])))      }
    }
    if (sim.list[[d]]$Guessing) {
      pty<-expression(hat(italic(c)) - italic(c))
      ptx<-expression(hat(italic(c)))
      sTitle<-eval(parse(text=paste0('expression("Condition"~"',gsub("S","",d),
                                     ' : "~italic(c)~": Bias & Error")')))
      q<-2+sim.list[[d]]$Q
      plot(SimList$gen.xi[,q],bias[,q],pch=19,cex=1.5,xlim=c(-0.01,1.05*max(SimList$gen.xi[,q])),
           ylim=1.05*range(c(bias[,q]-2*sXI[,q],bias[,q]+2*sXI[,q])),
           main=sTitle,xlab=ptx,ylab=pty)
      arrows(SimList$gen.xi[,q],bias[,q]-2*sXI[,q],SimList$gen.xi[,q],bias[,q]+2*sXI[,q],
             code=3,angle=90,length=0.07)
      abline(lm(bias[,q]~SimList$gen.xi[,q]),col=2,lty=2)
      print(summary(lm(bias[,q]~SimList$gen.xi[,q])))
    }
    if (sim.list[[d]]$Q==1) {
      par(mfrow=c(1,2),mar=c(5,5,4,4))
      library(plot3D)
      grid.size<-4
      xC<-seq(0,1.05*max(SimList$gen.xi[,1]),length.out = grid.size)
      dx<-mean(diff(xC))/2
      yC<-seq(-1.05*max(abs(SimList$gen.xi[,2])),1.05*max(abs(SimList$gen.xi[,2])),length.out = grid.size)
      dy<-mean(diff(yC))/2
      zzA<-matrix(NA,nrow = grid.size,ncol = grid.size)
      zzb<-matrix(NA,nrow = grid.size,ncol = grid.size)
      ErrorArray<-abs(XIhat-array(rep(mXI,sim.list[[d]]$Reps),
                                  dim=c(sim.list[[d]]$J,1+sim.list[[d]]$Q+sim.list[[d]]$Guessing,
                                        sim.list[[d]]$Reps)))
      for (xx in 1:length(xC)) {
        for (yy in 1:length(yC)) {
          if (sum(SimList$gen.xi[,1]>(xC[xx]-dx) & SimList$gen.xi[,1]<(xC[xx]+dx) &
                  SimList$gen.xi[,2]>(yC[yy]-dy) & SimList$gen.xi[,2]<(yC[yy]+dy))>0) {
            zzA[xx,yy]<-mean(ErrorArray[SimList$gen.xi[,1]>(xC[xx]-dx) & SimList$gen.xi[,1]<(xC[xx]+dx) &
                                          SimList$gen.xi[,2]>(yC[yy]-dy) & SimList$gen.xi[,2]<(yC[yy]+dy),1,],
                             na.rm=TRUE)
            zzb[xx,yy]<-mean(ErrorArray[SimList$gen.xi[,1]>(xC[xx]-dx) & SimList$gen.xi[,1]<(xC[xx]+dx) &
                                          SimList$gen.xi[,2]>(yC[yy]-dy) & SimList$gen.xi[,2]<(yC[yy]+dy),2,],
                             na.rm=TRUE)
          }
        }          
      }
      contour2D(z = zzA, x = xC, y=yC, xlab=expression(italic(A)), ylab=expression(italic(b)),
                main=expression("Contours of Mean Slope Errors"~ hat(italic(A))-italic(A)))
      contour2D(z = zzb, x = xC, y=yC, xlab=expression(italic(A)), ylab=expression(italic(b)),
                main=expression("Contours of Mean Intercept Errors"~ hat(italic(b))-italic(b)))
    }
  }
}




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
# apply(FitList$XI - FitList$xi,2,function(x) cbind(order(-abs(x)),sort(x,decreasing = T)))
# 
# for ( i in 2:ncol(FitList$XI)) lines(density(FitList$XI[,i] - FitList$xi[,i]),col=i)
# if (exists("FitList$TAU") && !is.na(FitList$TAU[1])) {
#   print(FitList$TAU - FitList$tau)
#   plot(density(as.vector(FitList$TAU) - as.vector(FitList$tau)))
# }
# 
# 
# plot(FitList$Iterations$Ait[1,1,],type="n")
# lines(FitList$Iterations$Ait[1,1,])
# acf(FitList$Iterations$Bit[1,])
# eigen(FitList$EZZ)
