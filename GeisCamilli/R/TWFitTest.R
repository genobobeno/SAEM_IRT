TWFitTest<-function(fit.data,ratios=TRUE,tw=0.999,...){
  E=fit.data$settings$Adim+5
  if (!is.na(fit.data$settings$ncat) && fit.data$settings$ncat>2) {
    S<-fit.data$EZZ #-t(as.matrix(FitList$EZ))%*%as.matrix(FitList$EZ)
  } else {
    S<-fit.data$EZZ-as.matrix(fit.data$EZ)%*%t(as.matrix(fit.data$EZ))
  }
  gEV<-eigen(S, symmetric=TRUE)
  cEV<-eigen(cov2cor(S),symmetric = TRUE)
  # if (sum(d %in% c("S1","S2","S3"))>0 & di==d[1]) {
  #     par(mfrow=c(length(d),2),...)
  #     pCH<-matrix(c(0,2,5,15,17,18),nrow=3,ncol=2)
  #   } else if (sum(d %in% c("S4","S5","S6","S7","S8","S9"))>0 & di==d[1]) {
  par(mfrow=c(1,1+ratios),...)
  pCH<-c(24,22)
  pCH.s<-1.2
  pBG<-c(0,3)
  #   }
  Q<-1
  sTitle<-paste("Largest Eigenvalues of Fit\nAssumed",fit.data$settings$Adim,"Dimensions")
  #plot(Q:E,c(0.7,gEV$values[Q:(E-1)]),type="n",main=sTitle,xlab="number of factors",ylab="estimators",lty=1,pch=1)
  plot(Q:E,log(c(0.7,gEV$values[Q:(E-1)])),type="n",main=sTitle,xlab="number of factors",ylab="log(estimators)",lty=1,pch=1)
  gTW.TF<-TWTransform(gEV$values,p = ncol(fit.data$RP),n = nrow(fit.data$RP),ptest = T)
  cTW.TF<-TWTransform(cEV$values,p = ncol(fit.data$RP),n = nrow(fit.data$RP),ptest = T)
  tw.999<-InverseTW(prob=tw,p=ncol(fit.data$RP),n=nrow(fit.data$RP))
  print(TWTransform(gEV$values[Q:E],p = ncol(fit.data$RP),n = nrow(fit.data$RP),ptest = F))
  print(TWTransform(gEV$values[Q:E],p = ncol(fit.data$RP),n = nrow(fit.data$RP),ptest = T))
  lines(Q:E,log(gEV$values[Q:E]),lty=1)
  points(Q:E,log(gEV$values[Q:E]),pch=pCH[1],col=0)
  points(Q:E,log(gEV$values[Q:E]),pch=pCH[1],
         cex=pCH.s,bg=ifelse(gEV$values[Q:E]>tw.999,pBG[2],pBG[1]))
  lines(Q:E,log(cEV$values[Q:E]),lty=4)
  points(Q:E,log(cEV$values[Q:E]),pch=pCH[2],col=0)
  points(Q:E,log(cEV$values[Q:E]),pch=pCH[2],cex=pCH.s,
         bg=ifelse(cEV$values[Q:E]>tw.999,pBG[2],pBG[1]))
  abline(h=log(tw.999),lty=2)
  tw.exp<-eval(parse(text=paste0('expression(lambda["P=',tw,'"]("TW"[1]))')))
  ltext =  c(expression(lambda("S"[2])),expression(lambda("Cor"("S"[2]))),tw.exp)
  strwidth(ltext)
  legend("topright",legend = ltext,text.width = max(strwidth(ltext))*1.02,
         lty=c(1,4,2),pch=c(24,22,NA),pt.cex=c(pCH.s,pCH.s,NA))
  if (ratios) {
    plot(Q:E,c(0.85,gEV$values[Q:(E-1)]/gEV$values[(Q+1):E]),type="n",
         main="Adjacent Ratios of Eigenvalues",
         xlab="number of factors",ylab="estimators",lty=1,pch=1)
    gER<-gEV$values[1:(ncol(fit.data$RP)-1)]/gEV$values[2:ncol(fit.data$RP)]
    cER<-cEV$values[1:(ncol(fit.data$RP)-1)]/cEV$values[2:ncol(fit.data$RP)]
    lines(Q:E,gER[Q:E],lty=1)
    points(Q:E,gER[Q:E],pch=pCH[1],bg=ifelse(gEV$values[Q:E]>tw.999,pBG[2],pBG[1]),col=1,cex=pCH.s)
    lines(Q:E,cER[Q:E],lty=4)
    points(Q:E,cER[Q:E],pch=pCH[2],bg=ifelse(cEV$values[Q:E]>tw.999,pBG[2],pBG[1]),col=1,cex=pCH.s)
    ltext =  c(expression(lambda("S"[2])),expression(lambda("Cor"("S"[2]))),
               expression(lambda["P=0.999"]("TW"[1])))
    strwidth(ltext)
    legend("topright",legend = ltext,text.width = max(strwidth(ltext))*1.02,
           lty=c(1,4,2),pch=c(pCH,NA),pt.cex=c(rep(pCH.s,2),NA))
  }
}