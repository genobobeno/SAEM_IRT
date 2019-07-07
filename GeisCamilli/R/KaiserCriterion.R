KaiserCriterion<-function(J,N,S){
  sEV<-eigen(S, symmetric=TRUE)
  rEV<-eigen(cov2cor(S),symmetric = TRUE)
  g<-J/N
  s.sumL<-cumsum(c(0,sEV$values))
  r.sumL<-cumsum(c(0,rEV$values))
  s.L.ref = ((J-s.sumL[-length(s.sumL)])/(J-1:J+1))*(1+sqrt(g))^2
  r.L.ref = ((J-r.sumL[-length(r.sumL)])/(J-1:J+1))*(1+sqrt(g))^2
  s.L.ref[s.L.ref<1]<-1
  r.L.ref[r.L.ref<1]<-1
  print(paste("EKC test on Covariance:  n.Factors =",sum(sEV$values>s.L.ref)))
  print(paste("EKC test on Correlation: n.Factors =",sum(rEV$values>r.L.ref)))
}
