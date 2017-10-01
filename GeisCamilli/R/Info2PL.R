Info2PL <-
function(aa,bb,Post) {  # Read in estimate of abilities, Theta[n], EAP, MAP, whatev
  NJ<-mat.or.vec(length(Post),1)+1   # N X 1
  #JQ<-mat.or.vec(length(tt),1)+1  # Q X 1
  Z<-outer(Post,bb,'-')     # N X J
  print("dim Z")
  print(dim(Z))
  P<-twoPL(aa,bb,Post)  # N X J
  print("dim P")
  print(dim(P))
  A<-t(aa%*%t(NJ)) 
  print("dim A")
  print(dim(A))
  IA<-(1.7)^2*Z^2*(P)*(1-P)
  IB<-(1.7*A)^2*(P)*(1-P)
  IAB<-(-1)*IA*A/Z
  IA<-colSums(IA)
  IB<-colSums(IB)
  IAB<-colSums(IAB)
  IM<-matrix(c(IA,IB,IAB),nrow=length(aa),ncol=3)
  print("Information Matrix Done")
  print(IM)
  return (IM) # J X 6
}
