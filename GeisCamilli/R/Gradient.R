Gradient <-
function(j,xi,rp,pL,settings) {
  #  PL : N X Q, this returns a sum
  theta<-GetQuad(settings)
  NQ<-mat.or.vec(nrow(rp),1)+1    # N X 1
  JQ<-mat.or.vec(nrow(theta),1)+1  # Q X 1
  if (settings$guess) {
    aa<-xi[,1]
    bb<-xi[,2]
    cc<-xi[,3]
  } else {
    aa<-xi[,1]
    bb<-xi[,2]    
  }
  if (tolower(settings$icc)=="ogive") {
    Z<-NQ%*%t(aa[j]*theta[,1]-bb[j])
    TT<-NQ%*%t(theta[,1]) #N X Q
    p<-pnorm(aa[j]*theta[,1]-bb[j])-0.000001*sin(pi/(2*max(theta[,1]))*theta[,1]) # Q X 1
    P<-NQ%*%t(p)  # N X Q
    R<-as.matrix(rp)[,j]%*%t(JQ)  # N X Q
    gA<-pL*(R-P)/(P*(1-P))*dnorm(Z)*(TT)
    gB<-pL*(P-R)/(P*(1-P))*dnorm(Z)
    gA<-sum(rowSums(gA))
    gB<-sum(rowSums(gB))
    if (is.nan(gA)|is.nan(gB)) {
      print("Something was recorded as an Infinity...")
      print("p")
      print(p)
      print("A")
      print(aa)
      print("B")
      print(bb)
      print("Posterior")
      print(pL[1:10,])
    }    
    dL<-c(gA,gB)
  }
  if (tolower(settings$icc)=="logistic") {
    Z<-NQ%*%t(theta[,1]-bb[j])
    if (settings$guess) {
      p<-cc[j]+(1-cc[j])*exp(1.7*aa[j]*(theta[,1]-bb[j]))/(1+exp(1.7*aa[j]*(theta[,1]-bb[j]))) # Q X 1
    } else {
      p<-exp(1.7*aa[j]*(theta[,1]-bb[j]))/(1+exp(1.7*aa[j]*(theta[,1]-bb[j]))) # Q X 1
    }
    P<-NQ%*%t(p)  # N X Q
    R<-as.matrix(rp)[,j]%*%t(JQ)  # N X Q
    gA<-pL*(R-P)*(1.7)*Z
    gB<-pL*(R-P)*(-1.7)*aa[j]
    gA<-sum(rowSums(gA))
    gB<-sum(rowSums(gB))
    dL<-c(gA,gB)
    if (settings$guess) {
      gC<-pL*(1-P)/(1-cc[j])
      gC<-sum(rowSums(gC)) 
      dL<-c(dL,gC)
    }
  }
#   print("Gradient Calculated")
#   print(dL)
  return (dL)
}
