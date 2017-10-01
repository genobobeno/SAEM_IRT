InfoTO2PL <-
function(aa,bb,Post) {  # Read in estimate of abilities, Theta[n], EAP, MAP, whatev
  NJ<-mat.or.vec(length(Post),1)+1 # N X 1
  P<-OgiveTwoPL(aa,bb,Post)  # N X J
  A<-t(aa%*%t(NJ)) # N X J
  IN<-(1.7*A)^2*(P)*(1-P)
  IN<-rowSums(IN) # Sum over J
  return (IN) # N X 1
}
