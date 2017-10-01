InfoO2PL <-
function(aa,bb,Post) {  # Read in estimate of abilities, Theta[n], EAP, MAP, whatev
  NQ<-mat.or.vec(length(Post),1)+1    # N X 1
  JQ<-mat.or.vec(length(aa),1)+1  # J X 1
  #JQ<-mat.or.vec(length(tt),1)+1  # Q X 1
  AT<-t(as.matrix(aa)%*%t(as.matrix(Post))) 
  Bz<-NQ%*%t(bb)
  Z<-AT-Bz    
  print("dim Z")
  print(dim(Z))
  P<-OgiveTwoPL(aa,bb,Post)  # N X J
  print("dim P")
  print(dim(P))
  TT<-NQ%*%t(Post)    
  R<-as.matrix(xx)[,jj]%*%t(JQ)  # N X Q
  
  gA<-pL*(R-P)/(P*(1-P))*dnorm(Z)*(TT)
  gB<-pL*(P-R)/(P*(1-P))*dnorm(Z)
  gA<-sum(rowSums(gA))
  gB<-sum(rowSums(gB))
  
  
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
