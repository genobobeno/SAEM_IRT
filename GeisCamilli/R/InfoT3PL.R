InfoT3PL <-
function(aa,bb,cc,Post) {  # Read in estimate of abilities, Theta[n], EAP, MAP, whatev
  NJ<-mat.or.vec(length(Post),1)+1 # N X 1
  P<-threePL(aa,bb,cc,Post)  # N X J
  A<-t(aa%*%t(NJ)) # N X J
  C<-t(cc%*%t(NJ)) # N X J
  IN<-(1.7*A/(1-C))^2*(P-C)*(C-P)*(1-P)/P
  IN<-rowSums(IN) # Sum over J
  return (IN) # N X 1
}
