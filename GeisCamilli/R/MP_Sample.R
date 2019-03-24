MP_Sample<-function(gamma,p=NA,n=NA,samples) {
  if (is.na(p)) {p<-as.integer(gamma*n)} else {n<-as.integer(p/gamma)}
  gminus = (1-sqrt(gamma))^2
  gplus = (1+sqrt(gamma))^2
  EV = seq(gminus,gplus,length.out = 1000)
  dist<-(1/(2*pi*gamma*EV))*sqrt((EV-gminus)*(gplus-EV))
  EVals=do.call("rbind",lapply(1:samples,function(y) {    
    rMat=do.call("rbind",lapply(1:n,function(x) rnorm(p)))  # The random matrix
    eigen(cov(rMat))$values
  }))
  max.ev<-apply(EVals,1,max)
  list(ev=EVals,x.mp=EV,y.mp=dist,max.ev=max.ev)
}