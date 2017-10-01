GetThetaHat <-
function(aa,bb,cc,rp,tHat,zHat,w,prior,setting,R=NA) {
  #print("Theta Estimation GO!")
  #if (!is.na(R)) aa=R$loadings
  ZHAT<-zHat  # N x J
  if (setting$Adim>1) {
    THAT<-array(tHat,dim=c(nrow(rp),setting$Adim,1))
  } else {
    THAT<-tHat
  }                      
  #  if (tolower(settings$esttheta)=="mcmc") {
  for (i in 1:setting$nesttheta) {
    tHat<-SampT(aa=aa,bb=bb,zz=zHat,rp=rp,prior=prior)
    zHat<-SampZ(aa=aa,bb=bb,that=tHat,rp=rp,w=w)
    ZHAT<-cbind(ZHAT,zHat)
    ifelse(setting$Adim>1,THAT<-abind(THAT,tHat,along=3),THAT<-cbind(THAT,tHat))
  }
  tHat<-SampT(aa=aa,bb=bb,zz=zHat,rp=rp,prior=prior)
  ifelse(setting$Adim>1,THAT<-abind(THAT[,,-1],tHat,along=3),THAT<-cbind(THAT[,-1],tHat))
  if (setting$Adim>1) {
    THETA<-apply(THAT,c(1,2),mean)
    THETA<-cbind(THETA,THETA-apply(THAT,c(1,2),sd),THETA+apply(THAT,c(1,2),sd))
    colnames(THETA)<-paste(c("Theta","Theta-SE","Theta+SE"),rep(c(1:setting$Adim),rep(3,setting$Adim)))
  } else {
    THETA<-rowMeans(THAT)
    THETA<-cbind(THETA,THETA-apply(THAT,1,sd),THETA+apply(THAT,1,sd))
    colnames(THETA)<-c("Theta","Theta-SE","Theta+SE")
  }
  if (setting$thetamap) {
    TMAP<-ThetaMAP(aa=aa,bb=bb,cc=cc,rp=rp,settings=setting)
    TMAP<-as.matrix(TMAP)
    ifelse(setting$Adim==1,colnames(TMAP)<-"TMAP",colnames(TMAP)<-paste("TMAP",1:setting$Adim,sep=""))
    if (!is.na(R)) {
      aa=R$loadings
      TRMAP<-ThetaMAP(aa=aa,bb=bb,cc=cc,rp=rp,settings=setting)
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
