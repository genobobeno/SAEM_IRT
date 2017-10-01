GenerateB <-
function(j,bdist,bparams) {
  if (bdist=="norm") {
    b<-rnorm(j,bparams[1],bparams[2])
  } else if (bdist=="unif") {
    b<-runif(j,bparams[1],bparams[2])
  }
  return(b)
}
