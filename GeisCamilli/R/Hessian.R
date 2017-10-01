Hessian <-
function(j,xi,rp,pL,settings) {
  #  PL : N X Q, this returns a sum
  theta<-GetQuad(settings)
  NQ<-mat.or.vec(nrow(rp),1)+1   # N X 1
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
    TT<-NQ%*%t(theta[,1])
    p<-pnorm(aa[j]*theta[,1]-bb[j])-0.000001*sin(pi/(2*max(theta[,1]))*theta[,1]) # Q X 1
#     print("Hess-P")
#     print(c(aa[j],bb[j]))
#     print(p)
    P<-NQ%*%t(p)  # N X Q
    R<-as.matrix(rp)[,j]%*%t(JQ)  # N X Q
    dA<-dnorm(Z)*(TT)
    dB<-(-1)*dnorm(Z)
    dAdB<-Z*TT*dnorm(Z)
    ddA<-(-1)*Z^2*TT^2*dnorm(Z)
    ddB<-Z*dnorm(Z)
    ddLA<-pL/(P*(1-P))*((R-P)*ddA-((1-2*P)*(R-P)/(P*(1-P)) + 1)*dA^2)
    ddLB<-pL/(P*(1-P))*((R-P)*ddB-((1-2*P)*(R-P)/(P*(1-P)) + 1)*dB^2)
    ddLAB<-pL/(P*(1-P))*((R-P)*dAdB-((1-2*P)*(R-P)/(P*(1-P)) + 1)*dB*dA)
#     print("Derivatives")
#     print(ddLA[1:5,])
#     print(ddLB[1:5,])
#     print(ddLAB[1:5,])
    ddLA<-sum(rowSums(ddLA))     
    ddLB<-sum(rowSums(ddLB))     
    ddLAB<-sum(rowSums(ddLAB))     
    ddL<-c(ddLA,ddLB,ddLAB)
  } else if (tolower(settings$icc)=="logistic") {
    Z<-NQ%*%t(theta[,1]-bb[j])    # N X Q
    p<-exp(1.7*aa[j]*(theta[,1]-bb[j]))/(1+exp(1.7*aa[j]*(theta[,1]-bb[j]))) # Q X 1
    P<-NQ%*%t(p)  # N X Q
    R<-as.matrix(rp)[,j]%*%t(JQ)  # N X Q 
    if (settings$guess) {
      #   gA<-(1.7/(1-cc[j]))^2*pL*(R*cc[j]/P^2-1)*Z^2*(P-cc[j])*(1-P)                # My derivatives
      #   gB<-(1.7*aa[j]/(1-cc[j]))^2*pL*(R*cc[j]/P^2-1)*(P-cc[j])*(1-P)             # My derivatives
      #   gC<-pL/(1-cc[j])^2*(R/P-1-R*(1-P)/P^2)                                       # My derivatives
      #   gAB<-(-1.7*1.7*aa[j])/(P*(1-cc[j]))^2*pL*                                   # My derivatives
      #      (P-cc[j])*(1-P)*(Z*(R*cc[j]-P^2)-P*(R-P)*(P-cc[j]))                     # My derivatives
      #   gAC<-pL*(-1.7)*Z/(P*(1-cc[j]))^2*(R*cc[j]*(1-2*P)+P^2*(R+cc[j]-1))         # My derivatives
      #   gBC<-pL*1.7*aa[j]/(P*(1-cc[j]))^2*(R*cc[j]*(1-2*P)+P^2*(R-cc[j]-1))       # My derivatives
      gA<-(1.7/(1-cc[j]))^2*pL*(R*cc[j]/P^2-1)*Z^2*(P-cc[j])*(1-P)
      gB<-(1.7*aa[j]/(1-cc[j]))^2*pL*(R*cc[j]/P^2-1)*(P-cc[j])*(1-P)
      gC<-pL/(1-cc[j])^2*(R/P-1)*(R*(1-P)/P^2)
      gAB<-(-1.7)/(1-cc[j])*pL*(P-cc[j])*((R/P-1)+1.7*aa[j]/(1-cc[j])*Z*(1-P)/P*(R*cc[j]/P-P))
      gAC<-(-1.7)/(1-cc[j])^2*pL*Z*(P-cc[j])*(1-P)/P^2*R
      gBC<-1.7*aa[j]/(1-cc[j])^2*pL*(P-cc[j])*(1-P)/P^2*R
      ddA<-sum(rowSums(gA))
      ddB<-sum(rowSums(gB))
      ddC<-sum(rowSums(gC))
      dAdB<-sum(rowSums(gAB))
      dAdC<-sum(rowSums(gAC))
      dBdC<-sum(rowSums(gBC))
      ddL<-c(ddA,ddB,ddC,dAdB,dAdC,dBdC)
    } else {
      #   gA<-pL*(R-P)*Z                                                    # Jimmys 
      #   ddA<-(-1.7*1.7)*sum(rowSums(gA)^2)                                # Jimmys 
      #   gB<-pL*(R-P)                                                      # Jimmys 
      #   ddB<-(-1.7*1.7)*aa[jj]^2*sum(rowSums(gB)^2)                       # Jimmys 
      #   gAB1<-pL*(R-P)*Z                                                  # Jimmys 
      #   gAB2<-pL*(R-P)                                                    # Jimmys 
      #   dAdB<-(-1.7*1.7)*aa[jj]*sum(rowSums(gAB1)*rowSums(gAB2))          # Jimmys 
      gA<-pL*Z^2*P*(1-P)                              # EMAN
      ddA<-(-1.7*1.7)*sum(rowSums(gA))                # EMAN
      gB<-pL*P*(1-P)*aa[j]^2                         # EMAN
      ddB<-(-1.7*1.7)*sum(rowSums(gB))                # EMAN
      gAB<-pL*((R-P)-1.7*aa[j]*Z*P*(1-P))            # EMAN
      dAdB<-(-1.7)*sum(rowSums(gAB))                  # EMAN
      ddL<-c(ddA,ddB,dAdB)
    }
  }
#   print("Hessian computed")
#   print(ddL)
  return (ddL)
}
