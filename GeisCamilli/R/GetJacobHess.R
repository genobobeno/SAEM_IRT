GetJacobHess <-
function(A,B,Z,setting=settings) {
  ZZ<-t(Z); K = nrow(ZZ);   N = ncol(ZZ);   Q = setting$Adim
  B<-matrix(B,K,1)
  db = matrix(0,K,1); da = matrix(0,K,Q); Jac=matrix(0,K*(1+Q),1)
  T1 = 0; T2 = matrix(0,Q,Q); ZT <- matrix(0,K,Q)      
  dLdXi2<-matrix(rep(0,((Q+1)*K)^2),(Q+1)*K,(Q+1)*K) 
  dLdXi<-NULL
  ZBi<- ZZ + B%*%rep(1,N)            
  ZB <- matrix(rowSums(ZBi),K,1)   
  IAA   <- solve(diag(Q) + t(A)%*%A)
  Lamb  <- IAA%*%t(A)
  for (ii in 1:N)  {
    alpha <- Lamb%*%ZBi[,ii]
    T1 <- T1 + alpha
    T2 <- T2 + IAA + alpha%*%t(alpha)
    ZT <- ZT + matrix(ZBi[,ii],K,1)%*%t(alpha) 
  }
  # compute Jacobian             
  db <- -ZB + A %*% T1
  da <-  ZT - A %*% T2
  for (jj in 1:K) {         
    dlda <- da[jj,] ;       dldb <- db[jj]
    dLdXi<-c(dLdXi, dldb, dlda)   # append item k's derivatives to Jac 
  }
  Jac = dLdXi
  #compute Hessian
  for (jj in 1:K)  {           
    HJ<-rbind(c(-N,T1),cbind(T1,-T2)) 
    for (ii in 1:nrow(HJ)) dLdXi2[(jj-1)*(Q+1)+ii,(jj-1)*(Q+1)+1:(Q+1)]<-HJ[ii,]   
  }
  H = dLdXi2 
  return(list(Jacob=Jac,Hess=H))                
}
