MarcenkoPastur<-function(J,N,EVhat=NA,generate=FALSE,Samples=10,
                         plot=FALSE,dtype=list(dist="norm",pars=c(0,1)),
                         mpTest=FALSE,Tests=NA,nSpectra=5) {
  # N=5000;J=30;plot=TRUE;nSpectra=5
  # dtype=list(dist="mvnorm",pars=c(0,0,0,0),pars2=diag(4),coefs=c(1.5,1.1,1,0.8))
  # generate=TRUE;mpTest=FALSE;Tests=100;Samples=10;EVhat=NA;
  # dtype=list(dist="norm",pars=c(0,1))
  gamma = J/N   # items to examinees, columns(Z)/rows(Z)
  gminus = (1-sqrt(gamma))^2  # analytic minimum expected
  gplus = (1+sqrt(gamma))^2   # analytic maximum expected
  EV = seq(gminus,gplus,length.out = 1000)
  dist<-(1/(2*pi*gamma*EV))*sqrt((EV-gminus)*(gplus-EV))
  HighEV<-ifelse(gplus<2,4,6)
  if (plot) {
    par(mfrow=c(1,1))
    xLimits = c(-0.2,1.3*max(HighEV,ifelse(!is.na(EVhat[1]),max(EVhat),0)))
    plot(EV,dist,xlim=xLimits,ylim=c(0,2*max(dist)),type="n",
         main=paste0("gamma = ",gamma,"; ",Samples," samples"))
    lines(EV,dist,lwd=2)
    if (!is.na(EVhat[1])) lines(density(EVhat),lwd=2)
  }
  gPulls<-seq(gamma/2,2*gamma,length.out=nSpectra)
  library(RColorBrewer)
  cols<-rainbow(nSpectra,start = 0,end = 0.75)
  
  if (generate) {
    j<-J
    if (mpTest) {
      if (is.na(Tests)[1]) Tests<-100
      mpDims<-c()
      nSpectra<-Tests
      gPulls<-gamma
      Samples<-ceiling(J/j)
    }
    for (i in 1:nSpectra) {
      # i=1
      print(paste0("Iteration ",i))
      n<-ceiling(j/ifelse(mpTest,gamma,gPulls[i]))
      # j<-J
      # n<-N
      if (dtype$dist=="pois") {
        EVals=do.call("rbind",lapply(1:Samples,function(y) {    
          rMat=do.call("rbind",lapply(1:n,function(x) rpois(j,lambda = dtype$pars[1])))  # The random matrix
          eigen(cov(rMat))$values
        }))
      } else if (dtype$dist=="unif") {
        EVals=do.call("rbind",lapply(1:Samples,function(y) {    
          rMat=do.call("rbind",lapply(1:n,function(x) runif(j,min = dtype$pars[1],max=dtype$pars[2])))  # The random matrix
          eigen(cov(rMat))$values
        }))
      } else if (dtype$dist=="exp") {
        EVals=do.call("rbind",lapply(1:Samples,function(y) {    
          rMat=do.call("rbind",lapply(1:n,function(x) rexp(j,rate=dtype$pars[1])))  # The random matrix
          eigen(cov(rMat))$values
        }))
      } else if (dtype$dist=="norm") {
        t0<-Sys.time()
        EVals=do.call("rbind",lapply(1:Samples,function(y) {    
          rMat=do.call("rbind",lapply(1:n,function(x) rnorm(j,mean = dtype$pars[1],sd=dtype$pars[2])))  # The random matrix
          cat(":",y)
          eigen(cov(rMat))$values
        }))
        Sys.time()-t0
      } else if (dtype$dist=="mvnorm") {
        if (sum(c("pars2","coefs") %in% names(dtype))==2 & length(dtype$pars)==length(dtype$coefs)) {
          if (sum(dtype$coefs^2)-1>0.01) dtype$coefs<-dtype$coefs/sqrt(sum(dtype$coefs^2))
          print("multivariate norm will result in each element being a SUM of orthogonal factors")
          EVals=do.call("rbind",lapply(1:Samples,function(y) {    
            rMat=do.call("rbind",lapply(1:n,function(x) rowSums(matrix(dtype$coefs,nrow = j,ncol=length(dtype$coefs),byrow = TRUE)*rmvn(j,mu = dtype$pars,sigma=dtype$pars2))))  # The random matrix
            cat(":",y)
            eigen(cov(rMat))$values
          }))
        } else {print("Multivariate norm needs sigma matrix as list element 'pars2'")}
      } else {
        t0<-Sys.time()
        EVals=do.call("rbind",lapply(1:Samples,function(y) {    
          rMat=do.call("rbind",lapply(1:n,function(x) rnorm(j)))  # The random matrix
          cat(":",y)
          eigen(cov(rMat))$values
        }))
        Sys.time()-t0
      }
      if (plot) {
        lines(density(EVals),col=cols[i])
      }
      if (mpTest & !is.na(EVhat)[1]) {
        mpDims<-c(mpDims,sum(EVhat>max(EVals)))
      } else if (mpTest) {
        mpDims<-c(mpDims,sum(EVals>max(EV)))
      }
      rm(EVals)
      gc()
    }
  }
  if (mpTest & !is.na(EVhat)[1]) {
    print(paste0("Returning results of ",Tests," MP tests on submitted eigenvalues."))
    return(mpDims)
  } else if (mpTest) {
    print(paste0("Returning results of ",Tests," MP tests on SAMPLED Random Matrix eigenvalues."))
    return(mpDims)
  }
  if (generate & plot) {
    print(paste0("Returned asymptotic distribution of eigenvalues for gamma = ",gamma))
    return(NULL)
  } else {
    MPDensity<-matrix(c(EV,dist),nrow = 1000,ncol=2)
    colnames(MPDensity)<-c("EV","Density")
    return(MPDensity)
  }
}
  
  
