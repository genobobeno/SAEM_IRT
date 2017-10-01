GetQuad <-
function(settings) {  # Returns two columns (grid, weights)
  if (tolower(settings$quad=="manual")) {
    Grid<-seq(settings$gridbounds[1],settings$gridbounds[2],length.out=settings$nodes)
    return(cbind(Grid,dnorm(Grid,mean=settings$tmu, sd=settings$tsigma)))
  } else if  (tolower(settings$quad=="gauss-hermite")) {
    Grid<-gaussHermiteData(settings$nodes)
    return(cbind(Grid$x,Grid$w))
  } else {
    print("Quadrature method not chosen, or not among choices")
    print("Please choose manual, or gauss-hermite")
  }
}
