GenerateC <-
function(j,cdist,cparams) {
  library(truncnorm)
  if (cdist=="norm") {
    stopifnot(length(cparams)==2,cparams[1]>0,cparams[2]<0.3)
    c<-rtruncnorm(j,a=0,b=0.4,cparams[1],cparams[2])
  } else if (cdist=="unif") {
    stopifnot(length(cparams)==2,cparams[1]>=0,cparams[2]<0.7)
    c<-runif(j,cparams[1],cparams[2])
  }
  return(c)  
}
