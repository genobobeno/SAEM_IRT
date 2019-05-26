FixTHatSamples<-function(fitFile) {
  FitList<-readRDS(fitFile)
  if (FitList$settings$parallel & FitList$settings$cores>1) {
    if (get_os()=="windows") {
      cl <<- parallel::makePSOCKcluster(FitList$settings$cores)
      parallel::clusterSetRNGStream(cl, iseed = round(runif(6)*1001))
      #parallel::setDefaultCluster(cl)
      prl<-TRUE
    } else if (get_os()=="linux") {
      cl <<- parallel::makeCluster(FitList$settings$cores,type="FORK",outfile="ClusterOutput.txt")
      doParallel::registerDoParallel(cl)
      parallel::clusterSetRNGStream(cl,iseed = round(runif(6)*1001))
      #parallel::setDefaultCluster(cl)
      prl<-TRUE
    } else {
      cat("\nI have only been testing for 'linux' and 'windows' operating systems. 
        Please feel free to edit the code in StartParallel() to add yours.\n")
      prl<-FALSE
    } 
  }
  library(rlecuyer)		# for rand num generation
  if (get_os()=="windows") library(snow)			# for parallel processing
  library(GPArotation)	# for rotations
  library(mvnfast)		# for function mvrnorm
  library(psych)			# for ML factor analysis
  library(foreach)
  clusterEvalQ(cl,library(mvnfast))
  
  J=ncol(FitList$RP);N=nrow(FitList$RP);Q=FitList$settings$Adim
  # Starting values 
  ncat <- FitList$settings$ncat
  prior<-list(tmu=FitList$settings$tmu,tsigma=FitList$settings$tsigma)
  IQ <- diag(Q) 
  # Set up n x m x (ncat-1) data array of binary values for
  y.a = data.matrix(FitList$RP)
  # cumulative option indicators
  y.b = array(NA,c(J,ncat-1,N))
  for (i in 1:N) { 
    for (j in 1:J) {
      # for missing responses, item propensities = 9 
      if (y.a[i,j]==9)	{
        y.b[j,,i] <- rep(9,(ncat-1))
      } else {
        for (k in 0:(ncat-2))
        {		
          y.b[j,(k+1),i] <- ifelse(y.a[i,j]<=k,0,1) 
        }
      }
    }
  }
  Y <- y.b
  n1cat = (ncat-1)
  MY = lapply(1:(dim(Y)[1]),function(x) (which(Y[x,,]==9,arr.ind = TRUE)))
  MN = lapply(1:(dim(Y)[1]),function(x) (which(Y[x,,]!=9,arr.ind = TRUE)))
  R  = lapply(1:(dim(Y)[1]),function(x) Y[x,,])
  missList = lapply(MY,function(x) (list(miss = nrow(x),nmis = N*n1cat-nrow(x),mcol = (N*n1cat-nrow(x))/n1cat)))
  indL<-list();indU<-list()
  for (j in 1:J) {
    yA<-y.a[,j]
    ### This missing handling needs to be revisited... 
    # Best solution is a sample from a different kind of augmented Z... thresholds off?
    if (sum(yA==9)>0) {
      yA[yA==9]<-floor(mean(yA[yA!=9]))
    }
    yL <- matrix(yA+1,N,1)
    indL[[j]] <- cbind(1:N,yL); 		indU[[j]] <- cbind(1:N,yL + 1)
  }
  clusterExport(cl,c("Y","Q","n1cat","N","J","MY","MN","R","missList","indL","indU"),envir=environment())
  clusterEvalQ(cl,c("WrapX","WrapT","WrapZ","WrapTmv"))
  
  if ("tau" %in% names(FitList)) {
    FitList$settings$thetamap<-FALSE
    That<-GetThetaHat(aa=FitList$A,bb=FitList$B,cc=FitList$C,rp=FitList$RP,tHat=FitList$That,zHat=FitList$EZ,w=NA,
                      prior=prior,setting=FitList$settings,R=NA,
                      TAU=FitList$TAU)
  } else {
    FitList$settings$thetamap<-FALSE
    That<-GetThetaHat(aa=FitList$A,bb=FitList$B,cc=FitList$C,rp=FitList$RP,tHat=FitList$That,zHat=FitList$EZ,w=NA,
                      prior=prior,setting=FitList$settings,R=NA,
                      TAU=NA)
  }
  parallel::stopCluster(cl)
  FitList$That<-That$THETA
  saveRDS(FitList,fitFile)
}
