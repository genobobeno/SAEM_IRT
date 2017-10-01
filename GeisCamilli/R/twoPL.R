twoPL <-
function(aa,bb,tt) {       # gives back a single probability or a vector if THETA is a vector
  p<-1/(1+exp(-1.7*aa*(tt-bb)))
}
