###Fitting to Marchenko-Pastur



MPFit<-function(S) {
  MPFunc<-function(x,gamma,sigma) {
    formula(y~(1/(2*pi*sigma*gamma*x))*sqrt((x-sigma*((1-sqrt(gamma))^2))*(sigma*((1+sqrt(gamma))^2)-x)))
  }
  MPdens<-KernSmooth::bkde(eigen(S)$values)
  summary(nls(MPFit,start = list(gamma = dim(RP)[2]/dim(RP)[1],sigma=1),data=data.frame(x=MPdens$x,y=MPdens$y),
    algorithm = "port"))
}