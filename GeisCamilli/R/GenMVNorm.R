GenMVNorm <-
function(n,mu,sigma) {
  library("MASS")
  stopifnot(min(dim(as.matrix(sigma)))>1,length(mu)>1)
  t<-mvrnorm(n,mu,sigma)
  return(t)
}
