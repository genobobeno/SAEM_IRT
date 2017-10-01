threePL <-
function(aa,bb,cc,tt) {  # gives back a single probability or a vector if THETA is a vector
  p<-cc+(1-cc)/(1+exp(-1.7*aa*(tt-bb)))
}
