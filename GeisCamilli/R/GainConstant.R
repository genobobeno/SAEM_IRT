GainConstant <-
function(settings=NA) {
  if (settings$annealing == "pseudo") {
    RMwindow<-ceiling(settings$burnin*(0.2))
    c(rep(1,settings$burnin-RMwindow),
      runif(rep(1,ceiling(RMwindow/2)),
            min = 1-(1:(RMwindow/2)/RMwindow)*cos(1:(RMwindow/2))*cos(1:(RMwindow/2)), 
            max = 1.0),
      runif(rep(1,ceiling(RMwindow/2)),
            min = 0.5-(1:(RMwindow/2)/RMwindow)*cos(1:(RMwindow/2))*cos(1:(RMwindow/2)), 
            max = 1-(1:(RMwindow/2)/RMwindow)*cos(1:(RMwindow/2))*cos(1:(RMwindow/2))),
      1/(2:(ceiling(2/settings$eps)))^settings$estgain) 
  } else {
    c(rep(1,settings$burnin),
    1/(1:ceiling(2/settings$eps))^settings$estgain)
  }
}
