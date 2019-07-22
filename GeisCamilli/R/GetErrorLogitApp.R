GetErrorLogitApp <-
function(A,B,C=NA,TH,RP,setting=settings) {
  d<-setting$Adim
  N<-nrow(RP)
  J<-ncol(RP)
  Jacob  <- NULL    
  D = 1.7    
  AT<-t(A%*%t(as.matrix(TH))) # t( J x df  %*%  df x N ) = N x J
  Bz<-t(B%*%t(as.matrix(rep(1,nrow(as.matrix(TH))))))  # N x J
  if (!is.na(C)[1]) {
    P = C + (1-C)/(1+exp(-D*(AT-Bz)))
    Hess   <- matrix(rep(0,((d+2)*J)^2),(d+2)*J,(d+2)*J) 
  } else {
    P = 1/(1+exp(-D*(AT-Bz)))
    Hess   <- matrix(rep(0,((d+1)*J)^2),(d+1)*J,(d+1)*J) 
  }
  Q = 1-P
  if (d==1) {
    THM<-as.matrix(TH)%*%matrix(rep(1,J),nrow=1,ncol=J)-Bz
  } else {
    THA<-array(0,dim=c(N,J,d))
    for (i in 1:d) {
      THA[,,i]<-as.matrix(TH[,i])%*%matrix(1,nrow=1,ncol=J)
    }
  }
  if (d==1) {
#    ifelse(setting$guess,P<-ProbIRT(cbind(A,B/A,C),TH),P<-ProbIRT(cbind(A,B/A),TH))  # N x J
#     print("RP-P, THM, P, Q")
#     print((RP-P)[1:5,])
#     print(THM[1:5,])
#     print(Q[1:5,])
    if (is.na(C)[1]) {
      dA = D*colSums((RP-P)*THM)
      dB = D*A*colSums(P-RP)
      ddA = -D^2*colSums(THM*THM*Q*P)
      ddB = -D^2*colSums(P*Q)
      dAdB = -D^2*colSums(P*((RP/P-1)-D*AT*Q*THM))
      Jacob = as.vector(t(cbind(dB,dA)))
      for (j in 1:J) {
        HJ  = rbind(c(ddB[j],dAdB[j]),c(dAdB[j],ddA[j]))
        for (i in 1:nrow(HJ)) 
          Hess[(j-1)*2+i,(j-1)*2+1:2] <-HJ[i,]     
      }
    } else {
      dA = D/(1-C)*colSums((RP-P)*THM*(P-t(matrix(C,nrow=J,ncol=N)))/P)
      dB = D/(1-C)*A*colSums((P-RP)*(P-t(matrix(C,nrow=J,ncol=N)))/P)
      dC = 1/(1-C)*colSums((RP-P)/P)
      ddA = (D/(1-C))^2*colSums(THM*THM*Q*(P-t(matrix(C,nrow=J,ncol=N)))/P*(RP/P*t(matrix(C,nrow=J,ncol=N))-P))
      ddB = (A*D/(1-C))^2*colSums(Q*(P-t(matrix(C,nrow=J,ncol=N)))/P*(RP/P*t(matrix(C,nrow=J,ncol=N))-P))
      ddC = 1/(1-C)^2*colSums((RP-P)/P*(RP*Q/P^2))
      dAdB = -(D/(1-C))*colSums((P-t(matrix(C,nrow=J,ncol=N))) *
                                 ((RP/P-1)+D*AT*THM*Q/P/(1-t(matrix(C,nrow=J,ncol=N))) *
                                    (RP/P*t(matrix(C,nrow=J,ncol=N))-P)))
      dAdC = -(D/(1-C))*colSums((P-t(matrix(C,nrow=J,ncol=N)))*THM*Q*RP/P^2)
      dBdC = (D*A/(1-C))*colSums((P-t(matrix(C,nrow=J,ncol=N)))*Q*RP/P^2)
      Jacob = as.vector(t(cbind(dC,dB,dA)))
      for (j in 1:J) {
        HJ  = rbind(c(ddC[j],dBdC[j],dAdC[j]),c(dBdC[j],ddB[j],dAdB[j]),c(dAdC[j],dAdB[j],ddA[j]))
        for (i in 1:nrow(HJ)) 
          Hess[(j-1)*3+i,(j-1)*3+1:3] <-HJ[i,]     
      }
    }
  } else {
#     ifelse(setting$guess,P<-ProbIRT(cbind(A,B/A,C),TH),P<-ProbIRT(cbind(A,B/A),TH))  # N x J
#     Q = 1-P
    dB = D*colSums(P-RP)
    ddB = -D^2*colSums(P*Q)
    dA = mat.or.vec(J,d)
    ddA = mat.or.vec(J,d)
    dAdB = mat.or.vec(J,d)
    dAdA = mat.or.vec(J,d*(d-1)/2)
    for (i in 1:d) {
      dA[,i] = D*colSums((RP-P)*THA[,,i])
      ddA[,i] = -D^2*colSums(THA[,,i]*THA[,,i]*P*Q)
      dAdB[,i] = D^2*colSums(P*Q*THA[,,i])
    }
    DA2<-combn(1:d,m=2)
    if (length(DA2)==2) {
      dAdA = -D^2*colSums(THA[,,DA2[1]]*THA[,,DA2[2]]*P*Q)
    } else {
      for (i in 1:ncol(DA2)) dAdA[,i] = -D^2*colSums(THA[,,DA2[1,i]]*THA[,,DA2[2,i]]*P*Q)
    }
    Jacob = as.vector(t(cbind(dB,dA)))
    for (j in 1:J) {
      AJ  = diag(ddA[j,])
      if (length(DA2)==2) {
        AJ[1,2]=AJ[2,1]=dAdA[j]
      } else {
        for (i in 1:ncol(DA2)) {
          AJ[DA2[1,i],DA2[2,i]]<-dAdA[j,i]
        }
      }
      HJ  = rbind(c(ddB[j],dAdB[j,]),cbind(dAdB[j,],AJ))
      for (i in 1:nrow(HJ)) 
        Hess[(j-1)*(d+1)+i,(j-1)*(d+1)+1:(d+1)] <- HJ[i,]     
    }    
  }
  return(list(Jacob=Jacob,Hess=Hess))
}
