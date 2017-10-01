CompareFits <-
function(Gen=NA,Fit1,Fit2,Fit3=NA) {
  Err1<-abs(Fit1$xiError)^0.5
  Err2<-abs(Fit2$xiError)^0.5
  if (!is.na(Fit3)) Err3<-abs(Fit3$xiError)^0.5
  UBnds<-apply(rbind(Fit1$xi,Fit2$xi),2,max)  
  LBnds<-apply(rbind(Fit1$xi,Fit2$xi),2,min)  
  plot(Fit1$xi[,1],Fit1$xi[,2],pch=1,cex=1.2,col="red",main="B vs. A", xlim=c(LBnds[1]-0.2,UBnds[1]+0.2),ylim=c(LBnds[2]-0.2,UBnds[2]+0.2))
  arrows(Fit1$xi[,1]-2*Err1[,1],Fit1$xi[,2],Fit1$xi[,1]+2*Err1[,1],Fit1$xi[,2],code=3,col="red",angle=90,length=0.07)
  arrows(Fit1$xi[,1],Fit1$xi[,2]-2*Err1[,2],Fit1$xi[,1],Fit1$xi[,2]+2*Err1[,2],code=3,col="red",angle=90,length=0.07)
  points(Fit2$xi[,1],Fit2$xi[,2],pch=19,col="blue")
  arrows(Fit2$xi[,1]-2*Err2[,1],Fit2$xi[,2],Fit2$xi[,1]+2*Err2[,1],Fit2$xi[,2],code=3,col="blue",angle=90,length=0.07)
  arrows(Fit2$xi[,1],Fit2$xi[,2]-2*Err2[,2],Fit2$xi[,1],Fit2$xi[,2]+2*Err2[,2],code=3,col="blue",angle=90,length=0.07)
  legend("topright",c("SAEM","EM"),pch=c(1,19),col=c("red","blue"))
  if (!is.na(Gen)) {
    points(Gen$XI[,1],Gen$XI[,2],pch=13)    
    AX1<-Fit1$xi[,1]-Gen$XI[,1]
    BX1<-Fit1$xi[,2]-Gen$XI[,2]
    AX2<-Fit2$xi[,1]-Gen$XI[,1]
    BX2<-Fit2$xi[,2]-Gen$XI[,2]
    plot(c(AX1,BX1),pch=1,cex=1.2,col="red",main="Errors", ylim=1.5*range(c(AX1,BX1)))
    arrows(1:length(AX1),AX1-2*Err1[,1],1:length(AX1),AX1+2*Err1[,1],code=3,col="red",angle=90,length=0.07)
    arrows(length(AX1)+1:length(BX1),BX1-2*Err1[,2],length(AX1)+1:length(BX1),BX1+2*Err1[,2],code=3,col="red",angle=90,length=0.07)
    points(1:length(c(AX1,BX1)),c(AX2,BX2),pch=19,col="blue")
    arrows(1:length(AX1),AX2-2*Err2[,1],1:length(AX1),AX2+2*Err2[,1],code=3,col="blue",angle=90,length=0.07)
    arrows(length(AX1)+1:length(BX1),BX2-2*Err2[,2],length(AX1)+1:length(BX1),BX2+2*Err2[,2],code=3,col="blue",angle=90,length=0.07)
    legend("topright",c("SAEM","EM"),pch=c(1,19),col=c("red","blue"))
  }
  
}
