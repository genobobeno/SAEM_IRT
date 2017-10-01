Info3PL <-
function(aa,bb,cc,Post) {  # Read in estimate of abilities, Theta[n], EAP, MAP, whatev
  NJ<-mat.or.vec(length(Post),1)+1   # N X 1
  #JQ<-mat.or.vec(length(tt),1)+1  # Q X 1
  Z<-outer(Post,bb,'-')     # N X J
  print("dim Z")
  print(dim(Z))
  P<-threePL(aa,bb,cc,Post)  # N X J
  print("dim P")
  print(dim(P))
  A<-t(aa%*%t(NJ)) 
  print("dim A")
  print(dim(A))
  C<-t(cc%*%t(NJ))
  print("dim C")
  print(dim(C))
  IA<-(1.7/(1-C))^2*Z^2*(P-C)^2*(1-P)/P
  IB<-((1.7*A)/(1-C))^2*(P-C)^2*(1-P)/P
  IC<-1/(1-C)^2*(1-P)/P
  IAB<-(-1)*IA*A/Z
  IAC<-IC*1.7*Z*(P-C)
  IBC<-IC*(-1.7)*A*(P-C)
  IA<-colSums(IA)
  IB<-colSums(IB)
  IC<-colSums(IC)
  IAB<-colSums(IAB)
  IAC<-colSums(IAC)
  IBC<-colSums(IBC)
  IM<-matrix(c(IA,IB,IC,IAB,IAC,IBC),nrow=length(aa),ncol=6)
  print("Information Matrix Done")
  print(IM)
  return (IM) # J X 6
}
