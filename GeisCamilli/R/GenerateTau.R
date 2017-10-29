GenerateTau <-
function(j,ncat,taudist="unif",tauparams=c(-1,1)) {
  # Generate the Tau parameters. J X ncat
  # j=5;ncat=3;taudist="unif";tauparams=c(-1,1)  
  TAU = lapply(rep(ncat-1,j),function(x) {
    if (taudist=="norm") {
      y = rnorm(n = x,tauparams[1],tauparams[2])
    } else if (taudist=="beta") {
      y = rbeta(n = x-1,tauparams[1],tauparams[2])
    } else {
      y = runif(n = x,tauparams[1],tauparams[2])
    } 
    sort(y)
  })
  tau = as.matrix(do.call(rbind,TAU))
  return(tau)
}
