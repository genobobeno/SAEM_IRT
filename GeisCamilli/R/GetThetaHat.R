GetThetaHat <-  function(aa,bb,cc,rp,tHat,zHat,w,prior,setting,R=NA,TAU=NA,refList=NA) {
  #print("Theta Estimation GO!")
  if (!is.na(R) & "Th" %in% names(R) & setting$Adim>1) {
    for (i in 1:ncol(R$loadings)) {
      if (sum(R$loadings[,i])<0) R$loadings[,i]<- (-1)*R$loadings[,i]
    }
    aa<-R$loadings
  }
  ZHAT<-zHat  # N x J
  if (setting$Adim>1) {
    THAT<-array(tHat,dim=c(nrow(rp),setting$Adim,1))
  } else {
    THAT<-tHat
  }                      
  #  if (tolower(settings$esttheta)=="mcmc") {
  cat("\n *** Starting Theta Estimation *** \n")
  cat(paste(setting$nesttheta,"Draws \n"))
  if (is.na(TAU)[1]) {
    for (i in 1:setting$nesttheta) {
      if (i%%500==1) cat(" \n",i,"\t")
      if (i%%10==1) cat(".")
      if (i%%100==1) cat(":")
      tHat<-SampT(aa=aa,bb=bb,zz=zHat,rp=rp,prior=prior)
      #zHat<-SampZ(aa=aa,bb=bb,that=tHat,rp=rp,w=w)
      zHat<-SampZFast(aa=aa,bb=bb,that=tHat,srp=rp,w=w)
      ZHAT<-cbind(ZHAT,zHat)
      if (setting$Adim>1) {
        THAT<-abind(THAT,tHat,along=3)
      } else {
        THAT<-cbind(THAT,tHat)
      }
    }
  } else {
    ATA 		<- t(aa)%*%aa #*4
    BTB_INV	<- solve(diag(setting$Adim) + ATA)
    for (i in 1:setting$nesttheta) {
      X2		<- simplify2array(parSapply(cl,1:length(bb),WrapX,simplify=FALSE,A=aa,b=bb,
                                      d=TAU-matrix(rep(bb,ncol(TAU)),length(bb),ncol(TAU)),theta=tHat), higher=TRUE)	
      X3		<- t(apply(X2,c(1,3),mean))
      if (!setting$dbltrunc) {
        zHat		<- colMeans(X2)
      } else {
        X1 		<- parSapply(cl,1:length(bb),WrapZ,simplify=FALSE,A=aa,b=bb,
                          d=TAU-matrix(rep(bb,ncol(TAU)),length(bb),ncol(TAU)),theta=tHat)
        zHat		<- simplify2array(X1, higher=TRUE)	
      }
      if (setting$Adim==1) {
        tHat	<- as.matrix(parSapply(cl,1:nrow(zHat),WrapT,A=aa,Z = zHat,BTB_INV=BTB_INV,b=bb,dbltrunc=setting$dbltrunc))
      } else {
        tHat	<- t(parSapply(cl,1:nrow(zHat),WrapTmv,A=aa,Z = zHat,BTB_INV=BTB_INV,b=bb,dbltrunc=setting$dbltrunc))
      }

      if (setting$Adim>1) {
        THAT<-abind(THAT,tHat,along=3)
      } else {
        THAT<-cbind(THAT,tHat)
      }
    }
  }
  print("Finished MCMC theta Estimate")
  ifelse(setting$Adim>1,THAT<-abind(THAT[,,-1],tHat,along=3),THAT<-cbind(THAT[,-1],tHat))
  if (setting$Adim>1) {
    THETA<-apply(THAT,c(1,2),mean)
    THETA<-cbind(THETA,THETA-apply(THAT,c(1,2),sd),THETA+apply(THAT,c(1,2),sd))
    colnames(THETA)<-paste(rep(c("Theta","Theta-SE","Theta+SE"),c(setting$Adim,setting$Adim,setting$Adim)),
                           rep(c(1:setting$Adim),rep(3,setting$Adim)))
  } else {
    THETA<-rowMeans(THAT)
    THETA<-cbind(THETA,THETA-apply(THAT,1,sd),THETA+apply(THAT,1,sd))
    colnames(THETA)<-c("Theta","Theta-SE","Theta+SE")
  }
  if (setting$thetamap) {
    print("Running Theta MAP Estimate")
    TMAP<-ThetaMAP(aa=aa,bb=bb,cc=cc,rp=rp,settings=setting,TAU=TAU,That=THETA[,1:setting$Adim])
    TMAP<-as.matrix(TMAP)
    ifelse(setting$Adim==1,colnames(TMAP)<-"TMAP",colnames(TMAP)<-paste("TMAP",1:setting$Adim,sep=""))
    if (!is.na(R)) {
      aa=R$loadings
      TRMAP<-ThetaMAP(aa=aa,bb=bb,cc=cc,rp=rp,settings=setting,TAU=TAU,That=THETA[,1:setting$Adim])
      TRMAP<-as.matrix(TRMAP)
      ifelse(setting$Adim==1,colnames(TRMAP)<-"TRMAP",colnames(TRMAP)<-paste("TRMAP",1:setting$Adim,sep=""))
    } else {
      TRMAP=NA
    }
  } else {
    TMAP<-NA
    TRMAP=NA
  }
  return(list(THETA=THETA,TMAP=TMAP,TRMAP=TRMAP))
}
