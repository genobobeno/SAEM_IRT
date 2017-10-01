GetHessian <-
function(A,B,Z) { 
  ZZ<-t(Z)                                      # K x N
  K = nrow(ZZ);   N = ncol(ZZ);   Q = ncol(A)
  B<-matrix(B,K,1)                              # K x 1
  sg = matrix(0,K,K)
  H<-matrix(rep(0,((Q+1)*K)^2),(Q+1)*K,(Q+1)*K) 
  IAA<-ginv(diag(Q)+t(A)%*%A)
  Lamb<-IAA%*%t(A)                 
  ZBi<- ZZ + B%*%rep(1,N)           # K x N  ( z + b )_ik
  ZB <- matrix(rowSums(ZBi),K,1)    # K x 1  (sum(z + b)_i)_k
  for (i in 1:N) sg = sg + ZBi[,i]%*%t(ZBi[,i])
  beta <- N*IAA + Lamb%*%sg%*%t(Lamb)
  alphaN <- Lamb %*% as.matrix(ZB)       
  for (j in 1:K)
  {           
    HJ<-rbind(c(-N,alphaN),cbind(alphaN,-beta)) 
    for (i in 1:nrow(HJ)) H[(j-1)*(Q+1)+i,(j-1)*(Q+1)+1:(Q+1)]<-HJ[i,]   
  }
  return(H)
}
