AnalyzeTestData <-
function(RP,settings=settings) {
  stopwatch<-Sys.time()
  stopifnot(nrow(RP)>ncol(RP))
  rp<-RP  
  J<-ncol(rp)
  N<-nrow(rp)
  print(paste("*****************   Analyzing Test Data, J =",J,";  N =",N,"; Output =",settings$estfile,"  ******************"))
  Init<-InitializeParams(rp,settings=settings) # returns XI and THat
  if (tolower(settings$model)=="gifa") {
    Estimates<-GoGIFA(rp,init=Init,settings=settings) # returns xi, that, xiErr, thatErr, Arot
  } else if (tolower(settings$model)=="irt") {
    Estimates<-GoIRT(rp,init=Init,settings=settings) # returns xi, that, xiErr, thatErr
  } else {
    print(paste("Model",settings$model,"not implemented yet"))
  }
#   if (!is.na(settings$simfile)) {
#     ifelse(grepl("\\.[Rr][Dd][Aa]",settings$simfile),load(settings$simfile),load(paste(settings$simfile,".rda",sep="")))
#     plot(as.vector(Estimates$xi-gen.xi),main="Parameter Estimate differences",ylab="XI_hat - XI_gen",xlab="A(1...J), B(1...J)")
#     plot(as.vector(Estimates$theta-gen.theta),main="Theta Estimate differences",ylab="Theta_hat - Theta_gen",xlab="Theta (1...N)")  
#   }
}
