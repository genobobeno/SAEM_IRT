DrawALowerDiag <- function(covT,covTZ,Q,m,n)
{
  Atemp <- matrix(0,Q,Q)
  for (k in 1:(Q-1)) { 
    E 	= k; IQ	= diag(k)
    TT 	= covT[1:k,1:k]
    TZ 	= covTZ[1:k,k] 
    Atemp[k,1:E] <- solve((2*m)*IQ + n*TT)%*%TZ*n
  }
  E 	= Q
  IQ 	= diag(Q)  #  Q X Q
  TT	= covT     #  Q X Q
  Atemp2<-t(sapply(Q:m,function(k) t(solve((2*m)*IQ + n*TT)%*%(covTZ[,k])*n) ))
  # [Q X Q] * [Q X 1] = Q X 1
  list(Atemp=rbind(Atemp,Atemp2),Avec=NA)
}

