AnalyzeTestData <-
function(RP,settings=settings,verbose=FALSE,TargetA = NA,simple=FALSE,timed=FALSE) {
  timed<-list(TF=timed)
  if (simple) {
    if (!is.na(TargetA)) {
      settings$Adim<-ncol(TargetA)
      settings$tmu<-rep(0,ncol(TargetA))
      settings$tsigma<-diag(ncol(TargetA))
    }
    settings$esttheta<-FALSE
    settings$nesttheta<-NA
    settings$empiricalse<-FALSE
    settings$thetamap<-FALSE
    settings$record<-FALSE
    settings$ncat<-length(unique(RP))
  }
  stopwatch<-Timing()
  stopifnot(nrow(RP)>ncol(RP))
  if (settings$parallel) { # & settings$cores>1) {
    if (get_os()=="windows") {
      cl <<- parallel::makePSOCKcluster(settings$cores)
      parallel::clusterSetRNGStream(cl, iseed = round(runif(6)*1001))
      #parallel::setDefaultCluster(cl)
      prl<-TRUE
    } else if (get_os()=="linux") {
      cl <<- parallel::makeCluster(settings$cores,type="FORK",outfile="ClusterOutput.txt")
      doParallel::registerDoParallel(cl)
      parallel::clusterSetRNGStream(cl,iseed = round(runif(6)*1001))
      #parallel::setDefaultCluster(cl)
      prl<-TRUE
    } else {
      cat("\nI have only been testing for 'linux' and 'windows' operating systems. 
        Please feel free to edit the code in StartParallel() to add yours.\n")
      prl<-FALSE
      settings$cores<-1
    } 
  } else {
    prl<-FALSE
    settings$cores<-1
  } 
  ### From this ifelse statement on, all missing data should have a value of 9
  if (sum(is.na(RP))>0 | sum(RP==settings$missing,na.rm = TRUE)>0) {
    rp<-CodeMissing(rp=RP,settings=settings)
  } else {
    rp<-RP  
  }
  J<-ncol(rp)
  N<-nrow(rp)
  if (!is.na(settings$missing)[1]) {
    K=max(apply(rp,2,function(x) (length(unique(x[!is.na(x) & x!=settings$missing])))))
  } else {
    K=max(apply(rp,2,function(x) (length(unique(x[!is.na(x)])))))
  }

  if (K>2) settings$ncat <- K
  if (verbose) {
    print(paste("************   Analyzing Test Data, J =",
                J,";  N =",N,"; Output =",settings$estfile,"  *************"))
  }
  Init<-InitializePrior(rp,settings=settings) # returns XI and THat
  if (timed$TF) {
    cat(paste0("Initialized:", Timing(stopwatch)))    
    settings$timed.initialize<-Timing(stopwatch)
  } 
  if (tolower(settings$model)=="gifa" & (is.na(settings$ncat)|settings$ncat==2)) {
    Estimates<-GoGIFA(rp,init=Init,settings=settings,TargetA=TargetA,timed=timed) # returns xi, that, xiErr, thatErr, Arot
  } else if (tolower(settings$model)=="gifa" | (!is.na(settings$ncat)&settings$ncat>2)) {
    Estimates<-GoPolyGIFA(rp,init=Init,settings=settings,TargetA=TargetA,timed=timed) # returns xi, that, xiErr, thatErr, Arot
  } else if (tolower(settings$model)=="irt") {
    Estimates<-GoIRT(rp,init=Init,settings=settings) # returns xi, that, xiErr, thatErr
  } else {
    print(paste("Model",settings$model,"not implemented yet"))
  }
  if (timed$TF) {
    cat(paste0("Initialized:", Timing(stopwatch)))    
    settings$timed.Total<-Timing(stopwatch)
  } 
  if (prl) {
    # cl <- parallel::getDefaultCluster()
    parallel::stopCluster(cl)
  }
  cat("
      Model has been fit
      Time: ",Timing(stopwatch))
  #   if (!is.na(settings$simfile)) {
  #     ifelse(grepl("\\.[Rr][Dd][Aa]",settings$simfile),load(settings$simfile),load(paste(settings$simfile,".rda",sep="")))
  #     plot(as.vector(Estimates$xi-gen.xi),main="Parameter Estimate differences",ylab="XI_hat - XI_gen",xlab="A(1...J), B(1...J)")
  #     plot(as.vector(Estimates$theta-gen.theta),main="Theta Estimate differences",ylab="Theta_hat - Theta_gen",xlab="Theta (1...N)")  
  #   }
  Estimates
}
