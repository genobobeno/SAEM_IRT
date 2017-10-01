GetErrorsIRTApp <-
function(A,B,TH,RP,setting=settings) {
  q<-ncol(as.matrix(TH))
  k<-ncol(RP)
  Jacob  <- NULL    
  Hess   <- matrix(rep(0,((q+1)*k)^2),(q+1)*k,(q+1)*k) 
  D = 1.7    
  for (i in 1:ncol(RP)) {      
    sum1A = 0; sum1b = 0; sum2A = 0; sum2b = 0; sumAb = 0
    P  <- 1/(1+exp( -( D*(A[i]*TH - B[i]))))
    Q  <- 1 - P
    u  <- RP[,i]
    sum1A <- sum(  D*(u - P)*TH  )
    sum1b <- sum(- D*(u - P)        )
    sum2A <- sum(- D*D*P*Q*TH*TH)
    sum2b <- sum(- D*D*P*Q )
    sumAb <- sum(  D*D*P*Q*TH)
    Jacob <- c(Jacob, sum1b, sum1A)       
    HJ    <- rbind(c(sum2b,sumAb),cbind(sumAb,sum2A))
    for (ii in 1:nrow(HJ)) 
      Hess[(i-1)*(q+1)+ii,(i-1)*(q+1)+1:(q+1)] <-HJ[ii,] 
  }    
  return(list(Jacob=Jacob,Hess=Hess))
}
