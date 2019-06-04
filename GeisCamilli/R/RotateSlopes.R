RotateSlopes<-function(fit.data=NA,slopes=NA) {
  stopifnot(is.list(fit.data)|!is.na(slopes)[1])
  if (is.list(fit.data)) {
    A<-fit.data$A
  } else {
    A<-slopes
  }
  print("Doing a bifactor rotation")
  ARB<-bifactorT(A, Tmat=diag(ncol(A)), normalize=FALSE, eps=1e-5, maxit=1000)
  ARB$FlippedAfterRotation<-colSums(ARB$loadings)<0
  ARB$loadings<-matrix(rep(ifelse(colSums(ARB$loadings)<0,-1,1),
                           nrow(A)),ncol=ncol(A),nrow=nrow(A),byrow=TRUE)*ARB$loadings
  print(ARB)
  print("Doing a Varimax rotation")
  ARV<-Varimax(A, Tmat=diag(ncol(A)), normalize=FALSE, eps=1e-5, maxit=1000)
  print(ARV)
  print("Doing a Infomax rotation")
  ARV$FlippedAfterRotation<-colSums(ARV$loadings)<0
  ARV$loadings<-matrix(rep(ifelse(colSums(ARV$loadings)<0,-1,1),
                           nrow(A)),ncol=ncol(A),nrow=nrow(A),byrow=TRUE)*ARV$loadings
  ARI<-infomaxT(A, Tmat=diag(ncol(A)), normalize=FALSE, eps=1e-5, maxit=1000)
  print(ARI)
  ARI$FlippedAfterRotation<-colSums(ARI$loadings)<0
  ARI$loadings<-matrix(rep(ifelse(colSums(ARI$loadings)<0,-1,1),
                           nrow(A)),ncol=ncol(A),nrow=nrow(A),byrow=TRUE)*ARI$loadings
  list(Bifactor=ARB,Varimax=ARV,Infomax=ARI)
}
