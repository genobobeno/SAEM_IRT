GIFAEstimate <-
function(aa,bb,zz,tt,settings,gain=NA,w=NA,rp=NA,ez=NA,ezz=NA,EmpT=FALSE) {
  if (is.na(gain)[1]) gain<-1
  ttt<-Sys.time()
  C<-NA
  mT<-colMeans(as.matrix(tt))   # D x 1
  EZ<-colMeans(zz)              # J x 1
  EZZ<-t(zz)%*%zz/(nrow(zz)-1)  # J x J
  zz<-as.matrix(t(zz))          # J x N
  aa<-as.matrix(aa)             # J x D
  bb<-as.matrix(bb)             # J x 1
  N<-ncol(zz)
  J<-nrow(zz)
  #   print("Covariance of Z")
#   print(EZZ)
#   print("mu * transpose(mu)")
#   print(as.matrix(EZ)%*%t(as.matrix(EZ)))
  ifelse(length(settings$tmu)==1,L<-t(aa)/as.numeric(settings$tsigma+t(aa)%*%aa),L<-ginv(settings$tsigma+t(aa)%*%aa)%*%t(aa)) # D x J
#   print("Lambda")
#   print(L)
  #ifelse(ncol(aa)>1,alpha<-L%*%ZB,alpha<-as.matrix(L%*%ZB)) # D x N
  if (tolower(settings$fm)=="camilli" | EmpT) {
    if (is.na(ez)[1]) {
      S<-EZZ-as.matrix(EZ)%*%t(as.matrix(EZ))
#       print("First term")
#       print(S%*%t(L))
#       print("Second term")
#       print(settings$tsigma-L%*%aa)
#       print("Third term")
#       print(L%*%S%*%t(L))
      ifelse(ncol(aa)>1,A<-(S%*%t(L))%*%ginv(settings$tsigma-L%*%aa+L%*%S%*%t(L)),A<-(S%*%t(L))/as.numeric(settings$tsigma-L%*%aa+L%*%S%*%t(L)))
#       print("A loadings")
#       print(A)
    } else {
      #print("Burned, now in FA")
      #print(paste("Iteration",it))
      EZ<-ez+(EZ-ez)*gain
      EZZ<-ezz+(EZZ-ezz)*gain
      S<-EZZ-as.matrix(EZ)%*%t(as.matrix(EZ))
      ifelse(ncol(aa)>1,A<-(S%*%t(L))%*%ginv(settings$tsigma-L%*%aa+L%*%S%*%t(L)),A<-(S%*%t(L))/as.numeric(settings$tsigma-L%*%aa+L%*%S%*%t(L)))
#       print(c((which.max(EZZ-ezz)-1)%%((ncol(aa)+1)*nrow(aa)) + 1,ceiling(which.max(EZZ-ezz)/((ncol(aa)+1)*nrow(aa)))))
    }    
    B<-EZ*(-1)
  } else if (tolower(settings$fm)=="eigen") {
    if (is.na(ez)[1]) {
      S<-EZZ-as.matrix(EZ)%*%t(as.matrix(EZ))
      #ifelse(settings$Adim>1,A<-(S%*%t(L))%*%ginv(settings$tsigma-L%*%aa+L%*%S%*%t(L)),A<-(S%*%t(L))/as.numeric(settings$tsigma-L%*%aa+L%*%S%*%t(L)))
    } else {
      EZ<-ez+(EZ-ez)*gain
      EZZ<-ezz+(EZZ-ezz)*gain
      S<-EZZ-as.matrix(EZ)%*%t(as.matrix(EZ))
    }      
    S     = S - diag(J)    
    out   = eigen( S, symmetric=TRUE)
    Avec  = out$values
    Aload = out$vectors
    if (settings$Adim>1) {
      Aload<-as.matrix(apply(rbind(Aload[,1:settings$Adim],aa),2,
                             function(x) (if (mean(x[1:(length(x)/2)]*x[length(x)/2+1:(length(x)/2)])>0 & sum(x[1:(length(x)/2)])>0) {return(x[1:(length(x)/2)])} 
                                          else {return(-1*x[1:(length(x)/2)])})),nrow(EZZ),settings$Adim)
      A = Aload[,1:settings$Adim]%*%sqrt(diag(Avec[1:settings$Adim]))
    } else {
      ifelse(mean(Aload[,1]*aa)>0,Aload<-Aload[,1],Aload<-(-1)*Aload[,1])
      A = as.matrix(Aload*sqrt(Avec[1]))
    }
    B<-EZ*(-1)
    #ifelse(settings$Adim>1,A<-(S%*%t(L))%*%ginv(settings$tsigma-L%*%aa+L%*%S%*%t(L)),A<-(S%*%t(L))/as.numeric(settings$tsigma-L%*%aa+L%*%S%*%t(L)))
  } else if (tolower(settings$fm)=="new") {
    if (is.na(ez)[1]) {
      S<-EZZ-as.matrix(EZ)%*%t(as.matrix(EZ))
      #ifelse(settings$Adim>1,A<-(S%*%t(L))%*%ginv(settings$tsigma-L%*%aa+L%*%S%*%t(L)),A<-(S%*%t(L))/as.numeric(settings$tsigma-L%*%aa+L%*%S%*%t(L)))
      out = princomp( covmat=S, scores=FALSE, cor=FALSE)
      Avec = out$sd^2
      if (settings$Adim>1) {
        Aload<-as.matrix(apply(rbind(out$loadings[,1:settings$Adim],aa),2,
                             function(x) (if (mean(x[1:(length(x)/2)]*x[length(x)/2+1:(length(x)/2)])>0) {return(x[1:(length(x)/2)])} 
                                          else {return(-1*x[1:(length(x)/2)])})),nrow(EZZ),settings$Adim)
        A = Aload[,1:settings$Adim]%*%sqrt(diag(Avec[1:settings$Adim]))
      } else {
        ifelse(mean(out$loadings[,1]*aa)>0,Aload<-out$loadings[,1],Aload<-(-1)*out$loadings[,1])
        A = as.matrix(Aload*sqrt(Avec[1]))
      }
    } else {
#       print("Burned, now in FA")
#       print(paste("Iteration",it))
      EZ<-ez+(EZ-ez)*gain
      EZZ<-ezz+(EZZ-ezz)*gain
      S<-EZZ-as.matrix(EZ)%*%t(as.matrix(EZ))
#       print("S Eigenvalues")
#       print(princomp(covmat=S,scores=FALSE, cor=FALSE)$sd^2-1)
      #ifelse(ncol(aa)>1,A<-(S%*%t(L))%*%ginv(settings$tsigma-L%*%aa+L%*%S%*%t(L)),A<-(S%*%t(L))/as.numeric(settings$tsigma-L%*%aa+L%*%S%*%t(L)))
      #out = princomp( covmat=A%*%t(A), scores=FALSE, cor=FALSE)
      out = princomp( covmat=S, scores=FALSE, cor=FALSE)
      Avec = out$sd^2
      if (settings$Adim>1) {
        Aload<-as.matrix(apply(rbind(out$loadings[,1:settings$Adim],aa),2,
                             function(x) (if (mean(x[1:(length(x)/2)]*x[length(x)/2+1:(length(x)/2)])>0) {return(x[1:(length(x)/2)])} 
                                          else {return(-1*x[1:(length(x)/2)])})),nrow(EZZ),settings$Adim)
        A = Aload[,1:settings$Adim]%*%sqrt(diag(Avec[1:settings$Adim]))
      } else {
        ifelse(mean(out$loadings[,1]*aa)>0,Aload<-out$loadings[,1],Aload<-(-1)*out$loadings[,1])
        A = as.matrix(Aload*sqrt(Avec[1]))
      }
    }    
    B<-EZ*(-1)
  } else if (tolower(settings$fm)=="old") {
    if (is.na(ez)[1]) {
      S<-EZZ-as.matrix(EZ)%*%t(as.matrix(EZ))
      A<-as.matrix(fa(S,nfactors=ncol(as.matrix(aa)),fm="ml",covar=TRUE,rotate="none",scores="none")$loadings)
      #       print("A loadings")
      #       print(A)
    } else {
      EZ<-ez+(EZ-ez)*gain
      EZZ<-ezz+(EZZ-ezz)*gain
      S<-EZZ-as.matrix(EZ)%*%t(as.matrix(EZ))
      A<-as.matrix(fa(S,nfactors=ncol(as.matrix(aa)),fm="ml",covar=TRUE,rotate="none",scores="none")$loadings)
    }        
    B<-as.vector(A%*%as.matrix(mT)-EZ)  
  } else if (tolower(settings$fm)=="pca") {
    if (is.na(ez)[1]) {
      #  zz is J x N
      zb<-zz+bb%*%rep(1,ncol(zz))
      SZB<-zb%*%t(zb)/(nrow(zz)-1)  # J x J
      S<-EZZ-as.matrix(EZ)%*%t(as.matrix(EZ))#-diag(nrow(aa))
      #A<-as.matrix(princomp(t(zz),covmat=S,scores=FALSE)$loadings)
      out = princomp( covmat=S, scores=FALSE, cor=FALSE)
      Avec = out$sd^2
      Aload<-as.matrix(apply(out$loadings,2,function(x) (if (mean((x>0)+0)>0.5) {return(x)} else {return(-1*x)})),nrow(EZZ),settings$Adim)
      #       Aload = out$loadings
      A = Aload[,1:settings$Adim]%*%sqrt(diag(Avec[1:settings$Adim]-1))
#       print("PCA Loading Scale factors...")
#       print(sqrt(diag(Avec[1:settings$Adim])))
#       print(sqrt(diag(Avec[1:settings$Adim]-1)))
#       out = princomp(covmat=S, scores=FALSE)
#       Avec = out$sd^2
#       Aload = out$loadings
#       A = Aload[,1:settings$Adim]*sqrt(Avec[1:settings$Adim]-1)
#       print("A Loadings")
#       print(A)
    } else {
      EZ<-ez+(EZ-ez)*gain
      EZZ<-ezz+(EZZ-ezz)*gain
      S<-EZZ-as.matrix(EZ)%*%t(as.matrix(EZ))#-diag(nrow(aa))
      #A<-as.matrix(princomp(covmat=S,scores=FALSE)$loadings)

      out = princomp( covmat=S, scores=FALSE, cor=FALSE)
      Avec = out$sd^2
      Aload<-as.matrix(apply(out$loadings,2,function(x) (if (mean((x>0)+0)>0.5) {return(x)} else {return(-1*x)})),nrow(EZZ),settings$Adim)
#       Aload = out$loadings
      A = Aload[,1:settings$Adim]%*%sqrt(diag(Avec[1:settings$Adim]-1))

#       out = princomp( covmat=S, scores=FALSE)
#       Avec = out$sd^2
#       Aload = out$loadings
#       A = Aload[,1:settings$Adim]*sqrt(Avec[1:settings$Adim]-1)
    }     
    B<-EZ*(-1)
    #B<-as.vector(A%*%as.matrix(mT)-EZ)  
  } else {
    if (is.na(ez)[1]) {
      S<-EZZ-as.matrix(EZ)%*%t(as.matrix(EZ))
      A<-as.matrix(fa(S,nfactors=ncol(aa),fm=tolower(settings$fm),covar=TRUE,rotate="none",scores="none")$loadings)      
    } else {
      EZ<-ez+(EZ-ez)*gain
      EZZ<-ezz+(EZZ-ezz)*gain
      S<-EZZ-as.matrix(EZ)%*%t(as.matrix(EZ))
      A<-as.matrix(fa(S,nfactors=ncol(aa),fm=tolower(settings$fm),covar=TRUE,rotate="none",scores="none")$loadings)      
    }             
    B<-EZ*(-1)
  }
  if (!is.na(w)[1]&!is.na(rp)[1]) {
    C<-rep(0,ncol(w))
    s<-rep(0,ncol(w))
    t<-rep(0,ncol(w))
    for (j in 1:ncol(w)) {
      t[j]<-length(which(w[,j]==0))
      s[j]<-sum(rp[which(w[,j]==0),j])
      C[j]<-rbeta(1,s[j]+1,t[j]-s[j]+1)
      if (C[j]>0.45)  { 
        C[j]<-runif(1,min=0.01,max=0.4)
      }
      C[j]<-s[j]/t[j]+(C[j]-s[j]/t[j])*gain
    }
  }
  return(list(A=A,B=B,C=C,EZ=EZ,EZZ=EZZ))
}
