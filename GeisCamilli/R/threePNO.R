threePNO <-
function(aa,bb,cc,tt) { # gives back a single probability or a vector if THETA is a vector
  p<-cc+(1-cc)*pnorm(aa*tt-bb)
}
