GetJacobian <-
function(A,B,Z) {
  ZZ<-t(Z)                                      # K x N
  K = nrow(ZZ);   N = ncol(ZZ);   Q = ncol(A)
  B<-matrix(B,K,1)                              # K x 1
  db = matrix(0,K,1); da = matrix(0,K,Q); sg = matrix(0,K,K); dLdXi<-vector()
  IAA<-ginv(diag(Q)+t(A)%*%A)
  Lamb<-IAA%*%t(A)
  ZBi<- ZZ + B%*%rep(1,N)           # K x N  ( z + b )_ik
  ZB <- matrix(rowSums(ZBi),K,1)    # K x 1  (sum(z + b)_i)_k
  db <- -ZB + A%*% Lamb%*%ZB
  for (i in 1:N) sg = sg + ZBi[,i]%*%t(ZBi[,i])
  da <- ((diag(K)-A%*%Lamb)%*%sg - N*diag(K))%*%t(Lamb)
  for (j in 1:K)
  {         
    dlda <- da[j,] ; dldb <- db[j]
    dLdXi<-c(dLdXi,dldb,dlda)     # append item k's derivatives to jacobian 
  }
  return(dLdXi) 
}
