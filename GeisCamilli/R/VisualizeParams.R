VisualizeParams <- function(struct=NA,samples=1000) {
  struct<-CheckParams(struct)
  par(mfrow=c(2,3),mar=c(5,4,3,1))
  #1
  A<-as.matrix(GenerateA(samples,struct$Adim,tolower(struct$Adist),struct$Aparams))[,1]
  plot(density(A),main=paste0("Distribution of First Dimension Slopes\nSamples = ",samples),ylab="Distribution",xlab="A")
  #2  
  if (struct$Adim>1) {
    A<-GenerateA(samples,struct$Adim,tolower(struct$Adist),struct$Aparams)[,2]
    plot(density(A),main=paste0("Distribution of Higher Dimension slopes\nSamples = ",samples),ylab="Distribution",xlab="A")
  } else {
    plot(c(1,0),c(0,1),type="n",main="Only 1 dimension",xlab="",ylab="")
    text(0.5,0.5,"The model only requires one dimension, \nthus no distribution of higher dimensions.")
  }
  #3 & 4
  if (!is.na(struct$ncat) && struct$ncat!=2) {
    tau = GenerateTau(samples,ncat=struct$ncat,taudist=tolower(struct$taudist),tauparams=struct$tauparams)
    b<-as.vector(rowMeans(tau))
    plot(density(b),main=paste0("Distribution of intercepts\n Samples = ",samples),xlab="b",ylab="Distribution")
    plot(density(as.vector(tau)),main=paste0("Distribution of category intercepts\n Samples = ",samples),xlab="tau",ylab="Distribution")
  } else {
    tau<-NA
    b<-GenerateB(samples,tolower(struct$bdist),struct$bparams)
    plot(density(b),main=paste0("Distribution of intercepts\n Samples = ",samples),xlab="b",ylab="Distribution")
    plot(c(1,0),c(0,1),type="n",main="Dichotomous",xlab="",ylab="")
    text(0.5,0.5,"The model is dichotomous, \nthus no distribution of category intercepts.")
  }
  xi<-cbind(A,b)
  #5
  if (struct$guess) {
    c<-GenerateC(samples,tolower(struct$cdist),struct$cparams)
    xi<-cbind(xi,c)
    plot(density(c),main=paste0("Distribution of guessing parameters\n Samples = ",samples),xlab="c",ylab="Distribution")
  }  else {
    plot(c(1,0),c(0,1),type="n",main="Guessing",xlab="",ylab="")
    text(0.5,0.5,"The model does not include guessing, \nthus no distribution of guessing parameters.")
  }
  #6
  t<-GenerateTheta(samples,struct$tmu,struct$tsigma)
  plot(density(as.vector(t)),main=paste0("Distribution of latent factors\n Samples = ",samples),xlab="THETA",ylab="Distribution")
  # rp<-GenerateRP(xi=xi,theta=t,struct=struct,tau=tau)
}