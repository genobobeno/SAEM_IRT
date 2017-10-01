GetJacobHessTheta <-
function(A,B,Z,TH,setting=settings) {
  if (setting$Adim>1) {
    ZZ<-t(Z)  # J x N
    TH<-t(TH) # 2 x N
    K = nrow(ZZ);   N = ncol(ZZ);   Q = ncol(A); 
    p<-(Q+1)*K
    B<-matrix(B,K,1)   # J x 1
    da = matrix(0,K,Q); Jac=matrix(0,K*(1+Q),1)
    T1 = 0; T2 = matrix(0,Q,Q); ZT <- matrix(0,K,Q)
    dLdXi2<-matrix(rep(0,(p)^2),p,p) 
    dLdXi<-NULL
    ZBi<- ZZ + B%*%rep(1,N)  # J x N
    ATi<- A %*% TH           # J x N
    db <- matrix(rowSums(ATi - ZBi),K,1)
    da <- (ZBi - ATi) %*% t(TH)
    dldab <- matrix(rowSums(TH),Q,1)
    dda <- TH %*% t(TH)
    ddb <- N*(-1)
    for (jj in 1:K) { 
      dlda <- da[jj,] ; dldb <- db[jj]
      dLdXi<-c(dLdXi, dldb, dlda)   # append item k's derivatives to Jac 
    }
    #compute Hessian
    for (jj in 1:K)  {           
      HJ<-rbind(c(-N,dldab),cbind(dldab,-dda)) 
      for (ii in 1:nrow(HJ)) dLdXi2[(jj-1)*(Q+1)+ii,(jj-1)*(Q+1)+1:(Q+1)]<-HJ[ii,]   
    }
  } else {
    ZZ<-t(Z)  # J x N
    K = nrow(ZZ);   N = ncol(ZZ);  
    TH<-matrix(TH,1,N) # 1 x N
    B<-matrix(B,K,1)   # J x 1
    A = matrix(A,K,1)  # J x 1 
    Jac=matrix(0,K*2,1)
    dLdXi2<-matrix(rep(0,(2*K)^2),2*K,2*K) 
    dLdXi<-NULL
    ZBi<- ZZ + B%*%rep(1,N)  # J x N
    ATi<- A %*% TH           # J x N
    db <- matrix(rowSums(ATi - ZBi),K,1)
    T1<- matrix(1,K,1) %*% TH  # J x 1
    da <- matrix(rowSums((ZBi - ATi)*T1,K,1))
    dldab <- sum(TH)
    dda <- TH %*% t(TH)
    ddb <- N*(-1)
    for (jj in 1:K) { 
      dlda <- da[jj] ; dldb <- db[jj]
      dLdXi<-c(dLdXi, dldb, dlda)   # append item k's derivatives to Jac 
    }
    #compute Hessian
    for (jj in 1:K)  {           
      HJ<-rbind(c(-N,dldab),c(dldab,-dda)) 
      for (ii in 1:nrow(HJ)) dLdXi2[(jj-1)*2+ii,(jj-1)*2+1:2]<-HJ[ii,]   
    }    
  }
  Jac = dLdXi
  H = dLdXi2 
  return(list(Jacob=Jac,Hess=H))                
}
