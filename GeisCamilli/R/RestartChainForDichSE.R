RestartChainForDichSE <- function(FitDATA,IT=1000) {
  cat("\n")
  print(paste("Starting",IT,"iterations of the Chain"))
  rp<-FitDATA$RP
  settings=FitDATA$settings
  library(rlecuyer)		# for rand num generation
  if (get_os()=="windows") library(snow)			# for parallel processing
  library(GPArotation)	# for rotations
  library(mvnfast)		# for function mvrnorm
  library(psych)			# for ML factor analysis
  library(foreach)
  
  if (FitDATA$settings$parallel) { # & settings$cores>1) {
    if (get_os()=="windows") {
      cl <<- parallel::makePSOCKcluster(FitDATA$settings$cores)
      if (is.na(FitDATA$settings$ncat)[1] | FitDATA$settings$ncat==2) {
        clusterCall(cl, cppInit)
      }
      parallel::clusterSetRNGStream(cl, iseed = round(runif(6)*1001))
      #parallel::setDefaultCluster(cl)
      prl<-TRUE
    } else if (get_os()=="linux") {
      cl <<- parallel::makeCluster(FitDATA$settings$cores,type="FORK",outfile="ClusterOutput.txt")
      doParallel::registerDoParallel(cl)
      parallel::clusterSetRNGStream(cl,iseed = round(runif(6)*1001))
      #parallel::setDefaultCluster(cl)
      prl<-TRUE
    } else {
      cat("\nI have only been testing for 'linux' and 'windows' operating systems. 
        Please feel free to edit the code in StartParallel() to add yours.\n")
      prl<-FALSE
      FitDATA$settings$cores<-1
    } 
  } else {
    prl<-FALSE
    FitDATA$settings$cores<-1
  } 
  if (prl) {
    clusterEvalQ(cl,library(mvnfast))
  }
  
  prior<-list(tmu=settings$tmu,tsigma=settings$tsigma)
  J<-ncol(rp)
  N<-nrow(rp)
  rp<-data.matrix(rp)
  CV<-mat.or.vec(J,J)+1
  pxi<-J*(1+settings$Adim+settings$guess)
  indL<-indU<-list()
  library(mcmcse)
  if (settings$Adim>1) {
    if (!is.na(FitDATA$AR)[1]) {
      A=as.matrix(FitDATA$AR$loadings)
      THat=as.matrix(FitDATA$Trot[,1:settings$Adim])
    } else {
      A=as.matrix(FitDATA$A)
      THat=as.matrix(FitDATA$That[,1:settings$Adim])
    }
    B=FitDATA$B
    C=FitDATA$C
    if (settings$fm %in% c("new","pca")) settings$fm <- "camilli"
  } else {
    A=FitDATA$A
    B=FitDATA$B
    C=FitDATA$C
    THat=FitDATA$That[,1]
  }
  
  if (settings$guess) {
    W<-DrawW(aa=A,bb=B,cc=C,tt=THat,rp=rp)
    for (j in 1:J) {
      indL[[j]] <- cbind(1:N,W[,j]+1)
      indU[[j]] <- cbind(1:N,W[,j]+2)
    }
  } else {    
    for (j in 1:J) {
      indL[[j]] <- cbind(1:N,rp[,j]+1)
      indU[[j]] <- cbind(1:N,rp[,j]+2)
    }
  }
  if (settings$parallel) {
    clusterExport(cl,c("J","N","rp","indL","indU"),envir=environment())
    clusterEvalQ(cl,c("pMVNarma"))
  }
    
  Alines<-A
  Blines<-B
  prior<-list(tmu=settings$tmu,tsigma=settings$tsigma)
  J<-ncol(rp)
  N<-nrow(rp)
  ifelse(settings$guess,W<-DrawW(aa=A,bb=B,cc=C,tt=THat,rp=rp),W<-NA)  
  Z<-SampZFast(aa=A,bb=B,that=THat,indL=indL,indU=indU,srp=rp,w=W)    
  ASEiter<-array(A, dim=c(J,settings$Adim,1))    
  BSEiter<-matrix(B, nrow=J, ncol=1)
  ifelse(settings$guess,CSEiter<-matrix(C, nrow=J, ncol=1),CSEiter<-NA)
  LLSEiter<-vector()

  for (it in 1:IT) {
    if (it%%10==1) cat(".")
    if (it%%100==1) cat(":")
    if (it%%500==1) cat(" \n",it,"\t")
    if (settings$guess) {
      W<-DrawW(aa=A,bb=B,cc=C,tt=THat,rp=rp)
      for (j in 1:J) {
        indL[[j]] <- cbind(1:N,W[,j]+1)
        indU[[j]] <- cbind(1:N,W[,j]+2)
      }
      if (settings$parallel) {
        clusterExport(cl,c("indL","indU"),envir=environment())
      }
    } else {
      W<-NA  
    }
    #Z<-SampZ(aa=A,bb=B,that=THat,rp=rp,w=W)    
    Z<-SampZFast(aa=A,bb=B,that=THat,indL=indL,indU=indU,srp=rp,w=W,prl=settings$parallel)    
    LL<-GIFAFullLL(A,B,Z,THat,prior=prior)
    if (settings$guess) {
      PSI<-GIFAEstimate(aa=A,bb=B,zz=Z,tt=THat,settings=settings,w=W,rp=rp,EmpT=TRUE)    
    } else {
      PSI<-GIFAEstimate(aa=A,bb=B,zz=Z,tt=THat,settings=settings,EmpT=TRUE)
    }
    A<-PSI$A
    B<-PSI$B
    C<-PSI$C    
    ASEiter<-abind(ASEiter,as.matrix(A),along=3)    
    BSEiter<-cbind(BSEiter,B)
    if (settings$guess) CSEiter<-cbind(CSEiter,C)
    LLSEiter<-c(LLSEiter,LL)
    THat<-SampT(aa=A,bb=B,zz=Z,rp=rp,prior=prior,prl=settings$parallel,cores=settings$cores) 
  }
  list(Aiter=ASEiter,Biter=BSEiter,Citer=CSEiter,Liter=LLSEiter)
}
