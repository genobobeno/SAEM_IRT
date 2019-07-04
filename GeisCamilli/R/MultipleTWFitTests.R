MultipleTWFitTests<-function(fit.data.list,ratios=TRUE,tw=0.999,title="",plot.layout=TRUE,...) {
  if (length(fit.data.list)>5) (print("Probably too many Plotted objects"))
  Adim.max<-which.max(sapply(fit.data.list,function(x) x$settings$Adim))
  E=fit.data.list[[Adim.max]]$settings$Adim+3
  f.er.exp<-f.exp<-c()
  if (plot.layout) {
    par(mfrow=c(1,1+ratios),mar=c(5,4,3,2),...)
  } else {
    par(mar=c(5,4,3,2),...)
  }
  pCH<-c(22,23,21,24,25)[1:length(fit.data.list)]
  pCH.s<-c(1.1,1.2,1.2,1.1,1.1)[1:length(fit.data.list)]
  pBG<-c(0,3)
  max.gEV<-0
  max.gER<-0
  for (f.data in 1:length(fit.data.list)) {
    if (!is.na(fit.data.list[[f.data]]$settings$ncat) && fit.data.list[[f.data]]$settings$ncat>2) {
      S<-fit.data.list[[f.data]]$EZZ #-t(as.matrix(FitList$EZ))%*%as.matrix(FitList$EZ)
    } else {
      S<-fit.data.list[[f.data]]$EZZ-as.matrix(fit.data.list[[f.data]]$EZ)%*%t(as.matrix(fit.data.list[[f.data]]$EZ))
    }
    gEV<-eigen(S, symmetric=TRUE)
    gER<-gEV$values[1:(ncol(fit.data.list[[f.data]]$RP)-1)]/gEV$values[2:ncol(fit.data.list[[f.data]]$RP)]
    max.gEV<-max(gEV$values,max.gEV)
    max.gER<-max(gER,max.gER)
  }
  for (f.data in 1:length(fit.data.list)) {
    if (!is.na(fit.data.list[[f.data]]$settings$ncat) && fit.data.list[[f.data]]$settings$ncat>2) {
      S<-fit.data.list[[f.data]]$EZZ #-t(as.matrix(FitList$EZ))%*%as.matrix(FitList$EZ)
    } else {
      S<-fit.data.list[[f.data]]$EZZ-as.matrix(fit.data.list[[f.data]]$EZ)%*%t(as.matrix(fit.data.list[[f.data]]$EZ))
    }
    
    gEV<-eigen(S, symmetric=TRUE)
    cEV<-eigen(cov2cor(S),symmetric = TRUE)
    print(paste0("Percentage of S Variance By Eigenvalue (",fit.data.list[[f.data]]$settings$Adim,") : "))
    print(signif(gEV$values[1:13]/sum(gEV$values),digits=3))
    print(paste0("Percentage of cor(S) Variance By Eigenvalue (",fit.data.list[[f.data]]$settings$Adim,") : "))
    print(signif(cEV$values[1:13]/sum(cEV$values),digits=3))
    # if (sum(d %in% c("S1","S2","S3"))>0 & di==d[1]) {
    #     par(mfrow=c(length(d),2),...)
    #     pCH<-matrix(c(0,2,5,15,17,18),nrow=3,ncol=2)
    #   } else if (sum(d %in% c("S4","S5","S6","S7","S8","S9"))>0 & di==d[1]) {
    Q<-1
    #plot(Q:E,c(0.7,gEV$values[Q:(E-1)]),type="n",main=sTitle,xlab="number of factors",ylab="estimators",lty=1,pch=1)
    if (f.data==1) {
      sTitle<-paste0("Largest Eigenvalues of ",title," Fits")
      plot(Q:E,log(c(0.7,gEV$values[Q:(E-2)],max.gEV)),type="n",main=sTitle,xlab="number of factors",ylab="log(estimators)",lty=1,pch=1)
    }
    gTW.TF<-TWTransform(gEV$values,p = ncol(fit.data.list[[f.data]]$RP),n = nrow(fit.data.list[[f.data]]$RP),ptest = T)
    cTW.TF<-TWTransform(cEV$values,p = ncol(fit.data.list[[f.data]]$RP),n = nrow(fit.data.list[[f.data]]$RP),ptest = T)
    tw.999<-InverseTW(prob=tw,p=ncol(fit.data.list[[f.data]]$RP),n=nrow(fit.data.list[[f.data]]$RP))
    print(TWTransform(gEV$values[Q:E],p = ncol(fit.data.list[[f.data]]$RP),n = nrow(fit.data.list[[f.data]]$RP),ptest = F))
    print(TWTransform(gEV$values[Q:E],p = ncol(fit.data.list[[f.data]]$RP),n = nrow(fit.data.list[[f.data]]$RP),ptest = T))
    lines(Q:E,log(gEV$values[Q:E]),lty=1)
    points(Q:E,log(gEV$values[Q:E]),pch=pCH[[f.data]],col=0)
    points(Q:E,log(gEV$values[Q:E]),pch=pCH[[f.data]],cex=pCH.s[f.data],
           bg=ifelse(gEV$values[Q:E]>tw.999,pBG[2],pBG[1]))
    lines(Q:E,log(cEV$values[Q:E]),lty=4)
    points(Q:E,log(cEV$values[Q:E]),pch=pCH[[f.data]],col=0)
    points(Q:E,log(cEV$values[Q:E]),pch=pCH[[f.data]],cex=pCH.s[f.data],
           bg=ifelse(cEV$values[Q:E]>tw.999,pBG[2],pBG[1]))
    f.exp<-eval(parse(text=paste0('c(f.exp,expression(lambda["Q=',
                                  fit.data.list[[f.data]]$settings$Adim,'"]("S"[2])),expression(lambda["Q=',
                                  fit.data.list[[f.data]]$settings$Adim,'"]("Cor"("S"[2]))))')))
  }
  abline(h=log(tw.999),lty=2)
  
  tw.exp<-eval(parse(text=paste0('expression(lambda["P=',tw,'"]("TW"[1]))')))
  ltext =  c(f.exp,tw.exp)
  strwidth(ltext)
  legend("topright",legend = ltext,text.width = max(strwidth(ltext))*1.02,
         lty=c(rep(c(1,4),length(f.exp)/2),2),
         pch=c(rep(pCH,rep(2,length(pCH))),NA),pt.cex=c(rep(pCH.s,rep(2,length(pCH.s))),NA))
  if (ratios) {
    plot(Q:E,c(0.85,rep(1,E-Q-1),max.gER),type="n",
         main="Adjacent Ratios of Eigenvalues",
         xlab="number of factors",ylab="estimators",lty=1,pch=1)
    for (f.data in 1:length(fit.data.list)) {
      if (!is.na(fit.data.list[[f.data]]$settings$ncat) && fit.data.list[[f.data]]$settings$ncat>2) {
        S<-fit.data.list[[f.data]]$EZZ #-t(as.matrix(FitList$EZ))%*%as.matrix(FitList$EZ)
      } else {
        S<-fit.data.list[[f.data]]$EZZ-as.matrix(fit.data.list[[f.data]]$EZ)%*%t(as.matrix(fit.data.list[[f.data]]$EZ))
      }
      
      gEV<-eigen(S, symmetric=TRUE)
      cEV<-eigen(cov2cor(S),symmetric = TRUE)
      
      gER<-gEV$values[1:(ncol(fit.data.list[[f.data]]$RP)-1)]/gEV$values[2:ncol(fit.data.list[[f.data]]$RP)]
      cER<-cEV$values[1:(ncol(fit.data.list[[f.data]]$RP)-1)]/cEV$values[2:ncol(fit.data.list[[f.data]]$RP)]
      
      lines(Q:E,gER[Q:E],lty=1)
      points(Q:E,gER[Q:E],pch=pCH[f.data],cex=pCH.s[[f.data]],
             bg=ifelse(gEV$values[Q:E]>tw.999,pBG[2],pBG[1]))
      lines(Q:E,cER[Q:E],lty=4)
      points(Q:E,cER[Q:E],pch=pCH[f.data],cex=pCH.s[[f.data]],
             bg=ifelse(cEV$values[Q:E]>tw.999,pBG[2],pBG[1]))
      f.er.exp<-eval(parse(text=paste0('c(f.er.exp,expression("ER"["Q=',
                                       fit.data.list[[f.data]]$settings$Adim,'"]("S"[2])),expression("ER"["Q=',
                                       fit.data.list[[f.data]]$settings$Adim,'"]("Cor"("S"[2]))))')))
      
    }
    ltext =  c(f.er.exp)
    strwidth(ltext)
    legend("topright",legend = ltext,text.width = max(strwidth(ltext))*1.02,
           lty=c(rep(c(1,4),length(f.er.exp)/2),2),
           pch=c(rep(pCH,rep(2,length(pCH)))),pt.cex=c(rep(pCH.s,rep(2,length(pCH.s)))))
  }
}
