PlotsCCF <-
function(Aiter,Biter,settings,ItCC,ItAC) {
  J<-nrow(Biter)
  It<-ncol(Biter)
  start<-1 #(It-40)
  par(mfrow=c(floor(J/5),floor(J/5)),mar=c(2,2,0.5,0.5))
  TotCC=0
  TotAC=rep(1,10)
  for (j in 1:J) {
    for (jj in 1:J) {
      if (j!=jj) {
        A1cf<-ccf(Aiter[j,1,start:It],Aiter[jj,1,start:It],lag.max = 10,cex=0.8,ylim=c(-1,1))
        text(-5,-0.5,paste("A1:",j,"-A1:",jj,sep=""))
        TotCC=TotCC+max(abs(A1cf$acf))
      } else {
        A1cf<-ccf(Aiter[j,1,start:It],Aiter[jj,2,start:It],lag.max = 10,cex=0.8,ylim=c(-1,1))
        text(-5,-0.5,paste("A1:",j,"-A2:",jj,sep=""))
        TotCC=TotCC+max(abs(A1cf$acf))
        A1cf<-acf(Aiter[jj,1,1:It],plot = FALSE,cex=0.8,ylim=c(-1,1),type="correlation")
        TotAC=TotAC*abs(A1cf$acf[2:11])
      }
    }
  }
  for (j in 1:J) {
    for (jj in 1:J) {
      if (j!=jj) {
        A2cf<-ccf(Aiter[j,2,start:It],Aiter[jj,2,start:It],lag.max = 10,cex=0.8,ylim=c(-1,1))
        text(-5,-0.5,paste("A2:",j,"-A2:",jj,sep=""))
        TotCC=TotCC+max(abs(A2cf$acf))              
      } else {
        A2cf<-ccf(Aiter[j,2,start:It],Aiter[jj,1,start:It],lag.max = 10,cex=0.8,ylim=c(-1,1))
        text(-5,-0.5,paste("A2:",j,"-A1:",jj,sep=""))
        TotCC=TotCC+max(abs(A2cf$acf))
        A2cf<-acf(Aiter[jj,2,1:It],plot = FALSE,cex=0.8,ylim=c(-1,1),type="correlation")
        TotAC=TotAC*abs(A2cf$acf[2:11])
      }
    }
  }
  for (j in 1:J) {
    for (jj in 1:J) {
      A1cf<-ccf(Aiter[j,1,start:It],Biter[jj,start:It],lag.max = 10,cex=0.8,ylim=c(-1,1))
      text(-5,-0.5,paste("A1:",j,"-B:",jj,sep=""))
      TotCC=TotCC+max(abs(A1cf$acf))
    }
  }
  for (j in 1:J) {
    for (jj in 1:J) {
      A2cf<-ccf(Aiter[j,2,start:It],Biter[jj,start:It],lag.max = 10,cex=0.8,ylim=c(-1,1))
      text(-5,-0.5,paste("A2:",j,"-B:",jj,sep=""))
      TotCC=TotCC+max(abs(A2cf$acf))              
    }
  }
  for (j in 1:J) {
    for (jj in 1:J) {
      if (j!=jj) {
        Bcf<-ccf(Biter[j,start:It],Biter[jj,start:It],lag.max = 10,cex=0.8,ylim=c(-1,1))
        text(-5,-0.5,paste("B:",j,"-B:",jj,sep=""))
        TotCC=TotCC+max(abs(Bcf$acf))
      } else {
        Bcf<-acf(Biter[jj,start:It],cex=0.8,ylim=c(-1,1),type="correlation")
        text(-5,-0.5,paste("ACF B:",jj,sep=""))
        TotAC=TotAC*abs(Bcf$acf[2:11])
      }
    }
  }
  ItCC<<-c(ItCC,TotCC)
  ItAC<<-c(ItAC,sum(TotAC))
  if (It>50) {
    par(mfrow=c(2,1))
    plot(ItCC,main="Total Cross Correlation starting from Iteration 30")
    plot(ItAC,main="Total Auto Correlation starting from Iteration 30")
  }
}
