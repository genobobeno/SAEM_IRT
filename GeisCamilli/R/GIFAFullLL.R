GIFAFullLL <-
function(aa,bb,zz,tt,prior=prior) {
  if (length(prior$tmu)<2) {   # zz is NxJ, tt is NxQ
    N<-nrow(zz)
    bb<-as.matrix(bb)
    at<-as.matrix(aa)%*%matrix(tt,1,N)
    a2<-sum(aa^2)
    LL<-(-0.5)*sum(cbind(rowSums(zz^2),as.vector(2*t(bb)%*%t(zz)),sum(bb^2),-2*colSums(t(zz+as.matrix(rep(1,N))%*%t(bb))*at),(a2+1)*tt^2))
  } else {
    N<-nrow(zz)
    bb<-as.matrix(bb)
    at<-as.matrix(aa)%*%t(tt)
    LL<-(-0.5)*sum(cbind(rowSums(zz^2),as.vector(2*t(bb)%*%t(zz)),rep(sum(bb^2),N),-2*colSums(t(zz+as.matrix(rep(1,N))%*%t(bb))*at),colSums(at^2),rowSums(tt^2)))    
  }
#   zz<-as.matrix(t(zz))  # J x N
#   aa<-as.matrix(aa)  # J x D
#   bb<-as.matrix(bb)  # J x 1
#   ifelse(length(prior$tmu)==1,IATA<-as.numeric(1/(prior$tsigma+t(aa)%*%aa)),IATA<-ginv(prior$tsigma+t(aa)%*%aa)) # D x D
#   ZB<-zz+bb%*%t(as.matrix(rep(1,ncol(zz)))) # J X N  
#   ifelse(ncol(aa)>1,alpha<-IATA%*%t(aa)%*%ZB,alpha<-IATA*(t(aa)%*%ZB)) # D x N
#   beta<-IATA+alpha%*%t(alpha) # D x D                     
#   if (ncol(aa)>1) {
#     LL<-0
#     for (i in 1:ncol(zz)) {
#       LL<-LL+t(zz[,i])%*%zz[,i]+2*t(bb)%*%zz[,i]+t(bb)%*%bb-2*t(ZB[,i])%*%aa%*%alpha[,i]+tr(t(aa)%*%aa%*%beta)+tr(beta)-0.5*sum(aa)/settings$Akappa #Final term is a prior on A
#     }
#     LL<-(-0.5)*LL
#   } else {
#     LL<-0
#     for (i in 1:ncol(zz)) {
#       LL<-LL+t(zz[,i])%*%zz[,i]+2*t(bb)%*%zz[,i]+t(bb)%*%bb-2*t(ZB[,i])%*%aa*alpha[,i]+(t(aa)%*%aa+1)*beta-0.5*sum(aa)/settings$Akappa #Final term is a prior on A
#     }    
#     LL<-(-0.5)*LL
#   }
  return(LL)
}
