OgiveICC <-
function(xi,guess=settings$guess,j=NA) { # gives back a N X J frame of probabilities... or N x 1
  #structure and settings have same list entries for calculations
  stopifnot((ncol(xi)<=4 & guess)|(ncol(xi)<4 & !guess))
  if (guess) {
    aa<-as.matrix(xi[,1:(ncol(xi)-2)])
    bb<-as.matrix(xi[,ncol(xi)-1])
    cc<-as.matrix(xi[,ncol(xi)])
    ifelse(is.na(j),C<-t(as.matrix(cc)%*%t(as.matrix(rep(1,nrow(as.matrix(theta)))))),C<-cc[j])
  } else {
    aa<-as.matrix(xi[,1:(ncol(xi)-1)])
    bb<-as.matrix(xi[,ncol(xi)])
  }
  if (ncol(aa)==2) {
    theta1 = -35:35*0.1
    theta2 = matrix(rep(-35:35*0.1,rep(71,71)),71,71)
    if (is.na(j)) {
      par(mfrow=c(5,4))
      for (j in 1:length(bb)) {
        if (j==1) par(mar=c(0,5,4,0))
        if (j>1 & j<4) par(mar=c(0,0,4,0))
        if (j==4) par(mar=c(0,0,4,3))
        if (j%in%c(5,9,13)) par(mar=c(0,5,0,0))
        if (j%in%c(6,7,10,11,14,15)) par(mar=c(0,0,0,0))
        if (j%in%c(8,12,16)) par(mar=c(0,0,0,3))
        if (j==17) par(mar=c(5,5,0,0))
        if (j>17 & j<20) par(mar=c(5,0,0,0))
        if (j==20) par(mar=c(5,0,0,3))
        for (i in 1:71) {
          theta = matrix(c(theta1,theta2[,i]),71,2)
          Bz<-bb[j]
          ifelse(ncol(aa)>1,AT<-(aa[j,]%*%t(as.matrix(theta))),AT<-aa[j]*theta)
          ifelse(!guess,p<-pnorm(AT-Bz),p<-C+(1-C)*pnorm(AT-Bz))
          if (i==1) plot(theta1,p,type="n",main=NA,ylab="P",xlab="Theta",xlim=range(theta),ylim=c(0,1))
          lines(theta1,p)
          if (i == 71) text(-2,0.9,paste("Item",j))
        }
      }
    } else {
      par(mfrow=c(1,1),mar=c(5,5,4,3))
      for (i in 1:71) {
        theta = matrix(c(theta1,theta2[,i]),71,2)
        Bz<-bb[j]
        ifelse(ncol(aa)>1,AT<-(aa[j,]%*%t(as.matrix(theta))),AT<-aa[j]*theta)
        ifelse(!guess,p<-pnorm(AT-Bz),p<-C+(1-C)*pnorm(AT-Bz))
        if (i==1) plot(theta1,p,type="n",main=paste("ICC for Item",j),ylab="P",xlab="Theta",xlim=range(theta),ylim=c(0,1))
        lines(theta1,p)
      } 
      ifelse(!guess,p<-pnorm(AT-Bz),p<-C+(1-C)*pnorm(AT-Bz))
    }
  } else {
    theta = -35:35*0.1
    if (is.na(j)) {
      AT<-t(aa%*%t(as.matrix(theta))) #J x df  %*%  df x N
      Bz<-t(bb%*%t(as.matrix(rep(1,nrow(as.matrix(theta))))))
      ifelse(!guess,p<-pnorm(AT-Bz),p<-C+(1-C)*pnorm(AT-Bz))
      par(mfrow=c(5,4))
      for (j in 1:length(bb)) {
        if (j==1) par(mar=c(0,5,4,0))
        if (j>1 & j<4) par(mar=c(0,0,4,0))
        if (j==4) par(mar=c(0,0,4,3))
        if (j%in%c(5,9,13)) par(mar=c(0,5,0,0))
        if (j%in%c(6,7,10,11,14,15)) par(mar=c(0,0,0,0))
        if (j%in%c(8,12,16)) par(mar=c(0,0,0,3))
        if (j==17) par(mar=c(5,5,0,0))
        if (j>17 & j<20) par(mar=c(5,0,0,0))
        if (j==20) par(mar=c(5,0,0,3))
        plot(theta,p[,j],type="n",main=NA,ylab="P",xlab="Theta",xlim=range(theta),ylim=c(0,1))
        lines(theta,p[,j])
        text(-2,0.9,paste("Item",j))
      }
    } else {
      par(mfrow=c(1,1),mar=c(5,5,4,3))
      Bz<-bb[j]
      ifelse(ncol(aa)>1,AT<-(aa[j,]%*%t(as.matrix(theta))),AT<-aa[j]*theta)
      par(mfrow=c(1,1),mar=c(5,5,4,3))
      ifelse(!guess,p<-pnorm(AT-Bz),p<-C+(1-C)*pnorm(AT-Bz))
      plot(theta,p,type="n",main=paste("ICC for Item",j),ylab="P",xlab="Theta",xlim=range(theta),ylim=c(0,1))
      lines(theta,p) 
    }  
  }
}
