### Suzanne's Tools


### Create category of URM [1 if yes, 0 if no, NA if other]
TestX
CutCategories()
urm<-c("1","3","4","5")
nurm<-c("6","8")
IDF$URM<-NA
for (i in 1:nrow(IDF)) {
  if (IDF$ETHNIC_CD[i] %in% urm) IDF$URM[i]<-1
  if (IDF$ETHNIC_CD[i] %in% nurm) IDF$URM[i]<-0
}


### Create SES category: 1 in lowest group, 2 in middle, 3 in highest
IDF$SES<-findInterval(IDF$MHI,unique(quantile(IDF$MHI,seq(0,1.0,length.out=4),na.rm=TRUE)),all.inside=TRUE)




### Item Plots 
CompareDemo<-function(score=NA,item=NA,class=NA,year=NA,test=1) {
  DF<-IDF
  ColCut<-c("GENDER_CD","URM","SES")
  if (!is.na(class)) {
    print(paste("Removing students not in class:",class))
    Cut<-paste("LetterGrade",class,sep="")
    DF<-DF[which(!is.na(DF[,Cut])),]
  }
  if (!is.na(year)) {
    print(paste("Removing students not graduating in year:",year))
    Cut<-paste("LetterGrade",class,sep="")
    DF<-DF[which(DF[,"Year"]==as.character(year)),]
  }
  if (!is.na(score)) {
    print(paste("Getting rid of rows without values for:",score))
    Cut<-score
    DF<-DF[which(!is.na(DF[,score])),c(ColCut,score)]
    pvar<-score
    title<-paste("Test",score)
  } else {
    if (is.na(item)) print("You must call this function with an item or score column")
    pvar<-paste("C",item,sep="")
    DF<-DF[which(!is.na(DF[,pvar])),c(ColCut,pvar)]
    title<-paste("Item",item,":",TestX[item,"Item"])
  }

  CutList=list(All=1:nrow(DF),MCut=which(DF$GENDER_CD=="M"),FCut=which(DF$GENDER_CD=="F"),SES1=which(DF$SES==1),SES2=which(DF$SES==2),SES3=which(DF$SES==3),URM=which(DF$URM==1),NURM=which(DF$URM==0))
  length(CutList)

  if (!is.na(score)) {
    Marg<-rep(0,length(CutList))
    ErrM<-rep(0,length(CutList))
    nn<-rep(0,length(CutList))
    names(Marg)<-c("ALL","Male","Fem","SES1","SES2","SES3","URM","NURM")
    for(i in 1:length(CutList)) {
      nn[i]<-length(which(!is.na(DF[CutList[[i]],pvar])))
      Marg[i]<-mean(DF[CutList[[i]],pvar],na.rm=TRUE)
      ErrM[i]<-sd(DF[CutList[[i]],pvar],na.rm=TRUE)/sqrt(nn[i])
      #ErrM[i]<-sqrt(Marg[i]*(1-Marg[i])/length(which(!is.na(DF[CutList[[i]],pvar]))))
    }
    par(mfrow=c(1,1),mar=c(4,4,3,2))
    centers<-barplot(Marg,ylim=c(0,1.05),ylab="Percentage",las=2,main=paste("Marginal Averages of",title))
    arrows(centers, Marg-ErrM*2, centers, Marg+ErrM*2, lwd=2, angle=90,code=3,length=0.1)
    text(centers,Marg+ErrM*2+0.05,labels=nn,cex=0.7)
    print(data.frame(N=nn,AVG=format(Marg,digits=3),SD=format(ErrM,digits=3)))
    
    par(mfrow=c(2,2),mar=c(4,4,3,2))
    centers<-barplot(Marg,ylim=c(0,1.05),ylab="Percentage",las=2,main=paste("Marginal Averages of",title))
    arrows(centers, Marg-ErrM*2, centers, Marg+ErrM*2, lwd=2, angle=90,code=3,length=0.03)
    text(centers,Marg+ErrM*2+0.05,labels=nn,cex=0.7)
    print(data.frame(N=nn,AVG=format(Marg,digits=3),SD=format(ErrM,digits=3)))
  
    MSES<-rep(0,3)
    FSES<-rep(0,3)
    MSErr<-rep(0,3)
    FSErr<-rep(0,3)
    MSn<-rep(0,3)
    FSn<-rep(0,3)
    for (i in 1:3) {
      MSn[i]<-length(which(!is.na(DF[CutList$MCut[!is.na(match(CutList$MCut,CutList[[3+i]]))],pvar])))
      FSn[i]<-length(which(!is.na(DF[CutList$FCut[!is.na(match(CutList$FCut,CutList[[3+i]]))],pvar])))
      MSES[i]<-mean(DF[CutList$MCut[!is.na(match(CutList$MCut,CutList[[3+i]]))],pvar],na.rm=TRUE)
      FSES[i]<-mean(DF[CutList$FCut[!is.na(match(CutList$FCut,CutList[[3+i]]))],pvar],na.rm=TRUE)
      MSErr[i]<-sd(DF[CutList$MCut[!is.na(match(CutList$MCut,CutList[[3+i]]))],pvar],na.rm=TRUE)/sqrt(MSn[i])
      FSErr[i]<-sd(DF[CutList$FCut[!is.na(match(CutList$FCut,CutList[[3+i]]))],pvar],na.rm=TRUE)/sqrt(FSn[i])
      #MSErr[i]<-sqrt(MSES[i]*(1-MSES[i])/length(which(!is.na(DF[CutList$MCut[!is.na(match(CutList$MCut,CutList[[3+i]]))],pvar]))))
      #FSErr[i]<-sqrt(FSES[i]*(1-FSES[i])/length(which(!is.na(DF[CutList$FCut[!is.na(match(CutList$FCut,CutList[[3+i]]))],pvar]))))
    }
    plot(1:3,MSES,ylim=c(0,1),ylab="Precentage",xlab=NA,main=paste("Gender vs. SES",title),pch=16,col=4,xaxt="n")
    axis(1,at=1:3,label=c("SES1","SES2","SES3"),las=2)
    points(1:3,FSES,pch=16,col=2)
    lines(1:3,MSES,col=4)
    lines(1:3,FSES,col=2)
    arrows(1:3, MSES-MSErr*2, 1:3, MSES+MSErr*2, lwd=2, angle=90,code=3,length=0.03,col=4)
    arrows(1:3, FSES-FSErr*2, 1:3, FSES+FSErr*2, lwd=2, angle=90,code=3,length=0.03,col=2)
    print(data.frame(N=c(MSn,FSn),AVG=format(c(MSES,FSES),digits=3),SD=format(c(MSErr,FSErr),digits=3)))
    
    USES<-rep(0,3)
    NSES<-rep(0,3)
    USErr<-rep(0,3)
    NSErr<-rep(0,3)
    USn<-rep(0,3)
    NSn<-rep(0,3)
    for (i in 1:3) {
      USn[i]<-length(which(!is.na(DF[CutList$URM[!is.na(match(CutList$URM,CutList[[3+i]]))],pvar])))
      NSn[i]<-length(which(!is.na(DF[CutList$NURM[!is.na(match(CutList$NURM,CutList[[3+i]]))],pvar])))
      USES[i]<-mean(DF[CutList$URM[!is.na(match(CutList$URM,CutList[[3+i]]))],pvar],na.rm=TRUE)
      NSES[i]<-mean(DF[CutList$NURM[!is.na(match(CutList$NURM,CutList[[3+i]]))],pvar],na.rm=TRUE)
      USErr[i]<-sd(DF[CutList$URM[!is.na(match(CutList$URM,CutList[[3+i]]))],pvar],na.rm=TRUE)/sqrt(USn[i])
      NSErr[i]<-sd(DF[CutList$NURM[!is.na(match(CutList$NURM,CutList[[3+i]]))],pvar],na.rm=TRUE)/sqrt(NSn[i])
      #USErr[i]<-sqrt(USES[i]*(1-USES[i])/length(which(!is.na(DF[CutList$URM[!is.na(match(CutList$URM,CutList[[3+i]]))],pvar]))))
      #NSErr[i]<-sqrt(NSES[i]*(1-NSES[i])/length(which(!is.na(DF[CutList$NURM[!is.na(match(CutList$NURM,CutList[[3+i]]))],pvar]))))
    }
    plot(1:3,USES,ylim=c(0,1),ylab="Precentage",main=paste("N/URM vs. SES",title),pch=16,col="brown",xaxt="n",xlab=NA)
    axis(1,at=1:3,label=c("SES1","SES2","SES3"),las=2)
    points(1:3,NSES,pch=16,col=3)
    lines(1:3,USES,col="brown")
    lines(1:3,NSES,col=3)
    arrows(1:3, USES-USErr*2, 1:3, USES+USErr*2, lwd=2, angle=90,code=3,length=0.03,col="brown")
    arrows(1:3, NSES-NSErr*2, 1:3, NSES+NSErr*2, lwd=2, angle=90,code=3,length=0.03,col=3)
    print(data.frame(N=c(USn,NSn),AVG=format(c(USES,NSES),digits=3),SD=format(c(USErr,NSErr),digits=3)))
  
    UMF<-rep(0,2)
    NMF<-rep(0,2)
    UMFErr<-rep(0,2)
    NMFErr<-rep(0,2)
    UMFn<-rep(0,2)
    NMFn<-rep(0,2)
    for (i in 1:2) {
      UMFn[i]<-length(which(!is.na(DF[CutList$URM[!is.na(match(CutList$URM,CutList[[1+i]]))],pvar])))
      NMFn[i]<-length(which(!is.na(DF[CutList$NURM[!is.na(match(CutList$NURM,CutList[[1+i]]))],pvar])))
      UMF[i]<-mean(DF[CutList$URM[!is.na(match(CutList$URM,CutList[[1+i]]))],pvar],na.rm=TRUE)
      NMF[i]<-mean(DF[CutList$NURM[!is.na(match(CutList$NURM,CutList[[1+i]]))],pvar],na.rm=TRUE)
      UMFErr[i]<-sd(DF[CutList$URM[!is.na(match(CutList$URM,CutList[[1+i]]))],pvar],na.rm=TRUE)/sqrt(UMFn[i])
      NMFErr[i]<-sd(DF[CutList$NURM[!is.na(match(CutList$NURM,CutList[[1+i]]))],pvar],na.rm=TRUE)/sqrt(NMFn[i])
      #<-sqrt(UMF[i]*(1-UMF[i])/length(which(!is.na(DF[CutList$URM[!is.na(match(CutList$URM,CutList[[1+i]]))],pvar]))))
      #<-sqrt(NMF[i]*(1-NMF[i])/length(which(!is.na(DF[CutList$NURM[!is.na(match(CutList$NURM,CutList[[1+i]]))],pvar]))))
    }
    plot(1:2,UMF,ylim=c(0,1),ylab="Precentage",main=paste("N/URM vs. Gender",title),pch=16,col="brown",xaxt="n",xlab=NA)
    axis(1,at=1:2,label=c("Male","Female"),las=2)
    points(1:2,NMF,pch=16,col=3)
    lines(1:2,UMF,col="brown")
    lines(1:2,NMF,col=3)
    arrows(1:2, UMF-UMFErr*2, 1:2, UMF+UMFErr*2, lwd=2, angle=90,code=3,length=0.03,col="brown")
    arrows(1:2, NMF-NMFErr*2, 1:2, NMF+NMFErr*2, lwd=2, angle=90,code=3,length=0.03,col=3)
    print(data.frame(N=c(UMFn,NMFn),AVG=format(c(UMF,NMF),digits=3),SD=format(c(UMFErr,NMFErr),digits=3)))
    
  } else {
    Marg<-rep(0,length(CutList))
    ErrM<-matrix(rep(0,2*length(CutList)),nrow=length(CutList),ncol=2)
    nn<-rep(0,length(CutList))
    names(Marg)<-c("ALL","Male","Fem","SES1","SES2","SES3","URM","NURM")
    for(i in 1:length(CutList)) {
      nn[i]<-length(which(!is.na(DF[CutList[[i]],pvar])))
      Marg[i]<-mean(DF[CutList[[i]],pvar],na.rm=TRUE)
      ErrM[i,1]<-binom.confint(Marg[i]*nn[i], nn[i], conf.level = 0.95, methods = "wilson")$lower
      ErrM[i,2]<-binom.confint(Marg[i]*nn[i], nn[i], conf.level = 0.95, methods = "wilson")$upper
      #sqrt(Marg[i]*(1-Marg[i])/length(which(!is.na(DF[CutList[[i]],pvar]))))
    }
    par(mfrow=c(1,1),mar=c(4,4,3,2))
    centers<-barplot(Marg,ylim=c(0,1.05),ylab="Percentage",las=2,main=paste("Marginal Averages of",title))
    arrows(centers, ErrM[,1], centers, ErrM[,2], lwd=2, angle=90,code=3,length=0.1)
    text(centers,ErrM[,2]+0.05,labels=nn,cex=0.7)
    print(data.frame(N=nn,AVG=format(Marg,digits=3),lower=format(ErrM[,1],digits=3),upper=format(ErrM[,2],digits=3)))
    
    par(mfrow=c(2,2),mar=c(4,4,3,2))
    centers<-barplot(Marg,ylim=c(0,1.05),ylab="Percentage",las=2,main=paste("Marginal Averages of",title))
    arrows(centers, ErrM[,1], centers, ErrM[,2], lwd=2, angle=90,code=3,length=0.03)
    text(centers,ErrM[,2]+0.05,labels=nn,cex=0.7)
    
    MSES<-rep(0,3)
    FSES<-rep(0,3)
    MSn<-rep(0,3)
    FSn<-rep(0,3)
    MSErr<-matrix(rep(0,2*3),nrow=3,ncol=2)
    FSErr<-matrix(rep(0,2*3),nrow=3,ncol=2)
    for (i in 1:3) {
      MSn[i]<-length(which(!is.na(DF[CutList$MCut[!is.na(match(CutList$MCut,CutList[[3+i]]))],pvar])))
      FSn[i]<-length(which(!is.na(DF[CutList$FCut[!is.na(match(CutList$FCut,CutList[[3+i]]))],pvar])))
      MSES[i]<-mean(DF[CutList$MCut[!is.na(match(CutList$MCut,CutList[[3+i]]))],pvar],na.rm=TRUE)
      FSES[i]<-mean(DF[CutList$FCut[!is.na(match(CutList$FCut,CutList[[3+i]]))],pvar],na.rm=TRUE)
      MSErr[i]<-sqrt(MSES[i]*(1-MSES[i])/length(which(!is.na(DF[CutList$MCut[!is.na(match(CutList$MCut,CutList[[3+i]]))],pvar]))))
      FSErr[i]<-sqrt(FSES[i]*(1-FSES[i])/length(which(!is.na(DF[CutList$FCut[!is.na(match(CutList$FCut,CutList[[3+i]]))],pvar]))))
      MSErr[i,1]<-binom.confint(MSn[i]*MSES[i], MSn[i], conf.level = 0.95, methods = "wilson")$lower
      MSErr[i,2]<-binom.confint(MSn[i]*MSES[i], MSn[i], conf.level = 0.95, methods = "wilson")$upper
      FSErr[i,1]<-binom.confint(FSn[i]*FSES[i], FSn[i], conf.level = 0.95, methods = "wilson")$lower
      FSErr[i,2]<-binom.confint(FSn[i]*FSES[i], FSn[i], conf.level = 0.95, methods = "wilson")$upper
    }
    plot(1:3,MSES,ylim=c(0,1),ylab="Precentage",xlab=NA,main=paste("Gender vs. SES",title),pch=16,col=4,xaxt="n")
    axis(1,at=1:3,label=c("SES1","SES2","SES3"),las=2)
    points(1:3,FSES,pch=16,col=2)
    lines(1:3,MSES,col=4)
    lines(1:3,FSES,col=2)
    arrows(1:3, MSErr[,1], 1:3, MSErr[,2], lwd=2, angle=90,code=3,length=0.03,col=4)
    arrows(1:3, FSErr[,1], 1:3, FSErr[,2], lwd=2, angle=90,code=3,length=0.03,col=2)
    print(data.frame(N=c(MSn,FSn),AVG=format(c(MSES,FSES),digits=3),lower=format(c(MSErr[,1],FSErr[,1]),digits=3),upper=format(c(MSErr[,2],FSErr[,2]),digits=3)))
    
    USES<-rep(0,3)
    NSES<-rep(0,3)
    USErr<-matrix(rep(0,2*3),nrow=3,ncol=2)
    NSErr<-matrix(rep(0,2*3),nrow=3,ncol=2)
    USn<-rep(0,3)
    NSn<-rep(0,3)
    for (i in 1:3) {
      USn[i]<-length(which(!is.na(DF[CutList$URM[!is.na(match(CutList$URM,CutList[[3+i]]))],pvar])))
      NSn[i]<-length(which(!is.na(DF[CutList$NURM[!is.na(match(CutList$NURM,CutList[[3+i]]))],pvar])))
      USES[i]<-mean(DF[CutList$URM[!is.na(match(CutList$URM,CutList[[3+i]]))],pvar],na.rm=TRUE)
      NSES[i]<-mean(DF[CutList$NURM[!is.na(match(CutList$NURM,CutList[[3+i]]))],pvar],na.rm=TRUE)
      USErr[i,1]<-binom.confint(USn[i]*USES[i], USn[i], conf.level = 0.95, methods = "wilson")$lower
      USErr[i,2]<-binom.confint(USn[i]*USES[i], USn[i], conf.level = 0.95, methods = "wilson")$upper
      NSErr[i,1]<-binom.confint(NSn[i]*NSES[i], NSn[i], conf.level = 0.95, methods = "wilson")$lower
      NSErr[i,2]<-binom.confint(NSn[i]*NSES[i], NSn[i], conf.level = 0.95, methods = "wilson")$upper
    }
    plot(1:3,USES,ylim=c(0,1),ylab="Precentage",main=paste("N/URM vs. SES",title),pch=16,col="brown",xaxt="n",xlab=NA)
    axis(1,at=1:3,label=c("SES1","SES2","SES3"),las=2)
    points(1:3,NSES,pch=16,col=3)
    lines(1:3,USES,col="brown")
    lines(1:3,NSES,col=3)
    arrows(1:3, USErr[,1], 1:3, USErr[,2], lwd=2, angle=90,code=3,length=0.03,col="brown")
    arrows(1:3, NSErr[,1], 1:3, NSErr[,2], lwd=2, angle=90,code=3,length=0.03,col=3)
    print(data.frame(N=c(USn,NSn),AVG=format(c(USES,NSES),digits=3),lower=format(c(USErr[,1],NSErr[,1]),digits=3),upper=format(c(USErr[,2],NSErr[,2]),digits=3)))
    
    UMF<-rep(0,2)
    NMF<-rep(0,2)
    UMFErr<-matrix(rep(0,2*2),nrow=2,ncol=2)
    NMFErr<-matrix(rep(0,2*2),nrow=2,ncol=2)
    UMFn<-rep(0,2)
    NMFn<-rep(0,2)
    for (i in 1:2) {
      UMFn[i]<-length(which(!is.na(DF[CutList$URM[!is.na(match(CutList$URM,CutList[[1+i]]))],pvar])))
      NMFn[i]<-length(which(!is.na(DF[CutList$NURM[!is.na(match(CutList$NURM,CutList[[1+i]]))],pvar])))
      UMF[i]<-mean(DF[CutList$URM[!is.na(match(CutList$URM,CutList[[1+i]]))],pvar],na.rm=TRUE)
      NMF[i]<-mean(DF[CutList$NURM[!is.na(match(CutList$NURM,CutList[[1+i]]))],pvar],na.rm=TRUE)
      UMFErr[i,1]<-binom.confint(UMFn[i]*UMF[i], UMFn[i], conf.level = 0.95, methods = "wilson")$lower
      UMFErr[i,2]<-binom.confint(UMFn[i]*UMF[i], UMFn[i], conf.level = 0.95, methods = "wilson")$upper
      NMFErr[i,1]<-binom.confint(NMFn[i]*NMF[i], NMFn[i], conf.level = 0.95, methods = "wilson")$lower
      NMFErr[i,2]<-binom.confint(NMFn[i]*NMF[i], NMFn[i], conf.level = 0.95, methods = "wilson")$upper
    }
    plot(1:2,UMF,ylim=c(0,1),ylab="Precentage",main=paste("N/URM vs. Gender",title),pch=16,col="brown",xaxt="n",xlab=NA)
    axis(1,at=1:2,label=c("Male","Female"),las=2)
    points(1:2,NMF,pch=16,col=3)
    lines(1:2,UMF,col="brown")
    lines(1:2,NMF,col=3)
    arrows(1:2, UMFErr[,1], 1:2, UMFErr[,2], lwd=2, angle=90,code=3,length=0.03,col="brown")
    arrows(1:2, NMFErr[,1], 1:2, NMFErr[,2], lwd=2, angle=90,code=3,length=0.03,col=3)
    print(data.frame(N=c(UMFn,NMFn),AVG=format(c(UMF,NMF),digits=3),lower=format(c(UMFErr[,1],NMFErr[,1]),digits=3),upper=format(c(UMFErr[,2],NMFErr[,2]),digits=3)))    
  }
}




ItemNumbers<-function() {
  print(TestX[which(TestX$Category=="MATH_REASON"),c("ItemIndex","Item")] )
}