TWSplitSCData<-function(fit.data.list,extra.dimensions=TRUE,E=10,tw=0.999,sTitle="Real Data",...) {
  # extra.dimensions=TRUE; d="S4";r=1
  maxGE<-maxCE<-0
  for (i in 1:length(fit.data.list)) {
    FitList<-fit.data.list[[i]]
    if (!is.na(FitList$settings$ncat) && FitList$settings$ncat>2) {
      S<-FitList$EZZ #-t(as.matrix(FitList$EZ))%*%as.matrix(FitList$EZ)
    } else {
      S<-FitList$EZZ-as.matrix(FitList$EZ)%*%t(as.matrix(FitList$EZ))
    }
    gEV<-eigen(S, symmetric=TRUE)
    cEV<-eigen(cov2cor(S),symmetric = TRUE)
    maxGE<-max(c(gEV$values[-1],maxGE))
    maxCE<-max(c(cEV$values[-1],maxCE))
  }
  yGLIM<-c(0.6,1.05*maxGE)
  yCLIM<-c(0.6,1.05*maxCE)
  f.exp<-c()
  pCH<-c(22,23,21,24,25)
  pCH.s<-1.3*c(1,1.2,1.2,1,1)
  pBG<-c(0,3)

  for (i in 1:length(fit.data.list)) {
    if (!is.na(fit.data.list[[i]]$settings$ncat) && fit.data.list[[i]]$settings$ncat>2) {
      S<-fit.data.list[[i]]$EZZ #-t(as.matrix(FitList$EZ))%*%as.matrix(FitList$EZ)
    } else {
      S<-fit.data.list[[i]]$EZZ-as.matrix(fit.data.list[[i]]$EZ)%*%t(as.matrix(fit.data.list[[i]]$EZ))
    }
    gEV<-eigen(S, symmetric=TRUE)
    cEV<-eigen(cov2cor(S),symmetric = TRUE)
    
    print(paste0("Percentage of S Variance By Eigenvalue (",fit.data.list[[i]]$settings$Adim,") : "))
    print(signif(gEV$values[1:13]/sum(gEV$values),digits=3))
    print(paste0("Percentage of cor(S) Variance By Eigenvalue (",fit.data.list[[i]]$settings$Adim,") : "))
    print(signif(cEV$values[1:13]/sum(cEV$values),digits=3))
    
    gTW.TF<-TWTransform(gEV$values,p = ncol(fit.data.list[[i]]$RP),n = nrow(fit.data.list[[i]]$RP),ptest = T)
    #print(paste0(fitdir[i],": TW transform and p-value of Q and Q+1 and Q+2"))
    #cTW.TF<-TWTransform(cEV$values,p = ncol(fit.data.list[[i]]$RP),n = nrow(fit.data.list[[i]]$RP),ptest = T)
    #j.r<-sum(TW.TF>0.95)+1
    tw.999<-InverseTW(prob=tw,p=ncol(fit.data.list[[i]]$RP),n=nrow(fit.data.list[[i]]$RP))
    if (i==1) {
      Q<-2
      labNames <- c("Largest Eigenvalues of",paste0("From ",sTitle))
      pTitle<-bquote(.(labNames[1])~S[2]~.(labNames[2]))
      plot(Q:E,c(0.7,gEV$values[Q:(E-1)]),ylim=yGLIM,
           type="n",main=pTitle,xlab="number of factors",ylab="estimators",lty=1,pch=1)
    }
    gTW.TF<-TWTransform(gEV$values,p = ncol(fit.data.list[[i]]$RP),n = nrow(fit.data.list[[i]]$RP),ptest = T)
    tw.999<-InverseTW(prob=tw,p=ncol(fit.data.list[[i]]$RP),n=nrow(fit.data.list[[i]]$RP))
    print(paste("Tracy Widom value of significance:",tw.999))
    print(paste(i,"in list: TWTransform no ptest - 1 thru",E))
    print(TWTransform(gEV$values[1:E],p = ncol(fit.data.list[[i]]$RP),n = nrow(fit.data.list[[i]]$RP),ptest = F))
    print(paste(i,"in list: TWTransform with ptest - 1 thru",E))
    print(TWTransform(gEV$values[1:E],p = ncol(fit.data.list[[i]]$RP),n = nrow(fit.data.list[[i]]$RP),ptest = T))
    lines(Q:E,gEV$values[Q:E],lty=1)
    points(Q:E,gEV$values[Q:E],pch=pCH[i],col=0)
    points(Q:E,gEV$values[Q:E],pch=pCH[i],cex=pCH.s[i],
           bg=ifelse(gEV$values[Q:E]>tw.999,pBG[2],pBG[1]))
    f.exp<-eval(parse(text=paste0('c(f.exp,expression(lambda["Q=',
                                  fit.data.list[[i]]$settings$Adim,'"]("S"[2])))')))
  }
  abline(h=tw.999,lty=2)
  tw.exp<-eval(parse(text=paste0('expression(lambda["P=',apaformat(tw,digits=round(abs(log(1-tw,base = 10)))),'"]("TW"[1]))')))
  ltext =  c(f.exp,tw.exp)
  strwidth(ltext)
  legend("topright",legend = ltext,text.width = max(strwidth(ltext))*1.02,
         lty=c(rep(1,length(f.exp)),2),
         pch=c(pCH[1:length(f.exp)],NA),pt.cex=c(pCH.s[1:length(f.exp)],NA))
  
  f.exp<-c()
  for (i in 1:length(fit.data.list)) {
    if (!is.na(fit.data.list[[i]]$settings$ncat) && fit.data.list[[i]]$settings$ncat>2) {
      S<-fit.data.list[[i]]$EZZ #-t(as.matrix(FitList$EZ))%*%as.matrix(FitList$EZ)
    } else {
      S<-fit.data.list[[i]]$EZZ-as.matrix(fit.data.list[[i]]$EZ)%*%t(as.matrix(fit.data.list[[i]]$EZ))
    }
    cEV<-eigen(cov2cor(S),symmetric = TRUE)
    cTW.TF<-TWTransform(cEV$values,p = ncol(fit.data.list[[i]]$RP),n = nrow(fit.data.list[[i]]$RP),ptest = T)
    #j.r<-sum(TW.TF>0.95)+1
    tw.999<-InverseTW(prob=tw,p=ncol(fit.data.list[[i]]$RP),n=nrow(fit.data.list[[i]]$RP))
    if (i==1) {
      Q<-2
      # sTitle<-eval(parse(text=paste0("expression('",sTitle," : Largest'~",E,"~'Eigenvalues of'~S[2])")))
      # plot(Q:E,c(0.7,gEV$values[Q:(E-1)]),type="n",main=sTitle,xlab="number of factors",ylab="estimators",lty=1,pch=1)
      labNames <- c("Largest Eigenvalues of cor(",paste0(") From ",sTitle))
      pTitle<-bquote(.(labNames[1])*S[2]*.(labNames[2]))
      plot(Q:E,c(0.7,cEV$values[Q:(E-1)]),ylim=yCLIM,
           type="n",main=pTitle,xlab="number of factors",ylab="estimators",lty=1,pch=1)
    }
    cTW.TF<-TWTransform(cEV$values,p = ncol(fit.data.list[[i]]$RP),n = nrow(fit.data.list[[i]]$RP),ptest = T)
    tw.999<-InverseTW(prob=tw,p=ncol(fit.data.list[[i]]$RP),n=nrow(fit.data.list[[i]]$RP))
    print(TWTransform(cEV$values[Q:E],p = ncol(fit.data.list[[i]]$RP),n = nrow(fit.data.list[[i]]$RP),ptest = F))
    print(TWTransform(cEV$values[Q:E],p = ncol(fit.data.list[[i]]$RP),n = nrow(fit.data.list[[i]]$RP),ptest = T))
    lines(Q:E,cEV$values[Q:E],lty=4)
    points(Q:E,cEV$values[Q:E],pch=pCH[[i]],col=0)
    points(Q:E,cEV$values[Q:E],pch=pCH[[i]],cex=pCH.s[i],
           bg=ifelse(cEV$values[Q:E]>tw.999,pBG[2],pBG[1]))
    f.exp<-eval(parse(text=paste0('c(f.exp,expression(lambda["Q=',
                                  fit.data.list[[i]]$settings$Adim,'"]("cor"("S"[2]))))')))
  }
  abline(h=tw.999,lty=2)
  
  tw.exp<-eval(parse(text=paste0('expression(lambda["P=',apaformat(tw,digits=round(abs(log(1-tw,base = 10)))),'"]("TW"[1]))')))
  ltext =  c(f.exp,tw.exp)
  strwidth(ltext)
  legend("topright",legend = ltext,text.width = max(strwidth(ltext))*1.02,
         lty=c(rep(4,length(f.exp)),2),
         pch=c(pCH[1:length(f.exp)],NA),pt.cex=c(pCH.s[1:length(f.exp)],NA))
  
}    

