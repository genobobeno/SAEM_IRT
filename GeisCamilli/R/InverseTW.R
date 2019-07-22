InverseTW<-function(prob,p,n) {
  q<-qtw(prob,beta = 1)
  mu<-1/n*(sqrt(n-1/2)+sqrt(p-1/2))^2
  sig<-sqrt(mu/n)*(1/sqrt(n-1/2)+1/sqrt(p-1/2))^(1/3)
  q*sig+mu
}