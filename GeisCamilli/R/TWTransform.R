TWTransform<-function(ev,p,n,ptest=FALSE) {
  #ev = 1.1081814; p=100; n=5000
  mu<-1/n*(sqrt(n-1/2)+sqrt(p-1/2))^2
  sig<-sqrt(mu/n)*(1/sqrt(n-1/2)+1/sqrt(p-1/2))^(1/3)
  if (!ptest) {
    (ev-mu)/sig
  } else {
    sapply((ev-mu)/sig,function(x) ptw(q = x,beta = 1))
  }
}