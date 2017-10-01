Posterior <-
function(L,xi,rp,settings) {
  # N X Q
  theta<-GetQuad(settings)
  NQ<-mat.or.vec(nrow(rp),nrow(theta))
  JQ<-mat.or.vec(nrow(theta),1)+1
  if (tolower(settings$icc)=="ogive") {
    p<-ProbOgive(xi=xi,theta=theta[,1],settings=settings)-as.matrix(0.000001*sin(pi/(2*max(theta[,1]))*theta[,1]))%*%rep(1,ncol(rp))   # Q x J
  } else if (tolower(settings$icc)=="logistic") {
    p<-ProbIRT(xi=xi,theta=theta[,1],settings=settings)  # Q x J
  }
  for (n in 1:nrow(rp)) {
    X<-t(as.matrix(rp)[n,]%*%t(JQ))
    umL<-X*log(p)+(1-X)*log(1-p)
    NQ[n,]<-exp(rowSums(umL))*theta[,2]/L[n] # Q X 1
  } 
  return (NQ)  # N X Q
}
