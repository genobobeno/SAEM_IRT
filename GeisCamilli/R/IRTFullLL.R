IRTFullLL <-
function(xi,rp,settings) {
  theta<-GetQuad(settings)
  JQ<-mat.or.vec(nrow(theta),1)+1
  L<-rep(0,nrow(rp))
  #Calculate Full Likelihood
  if (tolower(settings$icc)=="ogive") {
    p<-ProbOgive(xi=xi,theta=theta[,1],settings=settings)-as.matrix(0.000001*sin(pi/(2*max(theta[,1]))*theta[,1]))%*%rep(1,ncol(rp))  # Q x J
  } else if (tolower(settings$icc)=="logistic") {
    p<-ProbIRT(xi=xi,theta=theta[,1],settings=settings)  # Q x J
  }
  for (n in 1:nrow(rp)) {
    X<-t(as.matrix(rp)[n,]%*%t(JQ)) # t(J X Q) = Q X J
    umL<-X*log(p)+(1-X)*log(1-p) # Q X J
    umL<-exp(rowSums(umL))*theta[,2] # Q X 1
    L[n]<-sum(umL)
  }
  print("Full Likelihood done")
  return(L)
}
