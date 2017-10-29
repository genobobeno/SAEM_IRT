GenerateTheta <-
function(n,mu,sigma) {
  library("MASS")
  stopifnot(length(mu)==nrow(as.matrix(sigma))) 
  if (length(mu)==1) {
    t<-rnorm(n,mu,sigma)
  } else {
    t<-GenMVNorm(n,mu,sigma)
  }
  return(t)
}
