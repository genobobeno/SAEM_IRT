GenerateRP <-
function(xi,theta,structure=structure,tau=NA) {
  RP<-mat.or.vec(nrow(as.matrix(theta)),nrow(as.matrix(xi)))
  if (tolower(structure$icc)=="ogive") {
    if (!is.na(structure$ncat) && structure$ncat!=2) {
      P<-ProbOgive(xi=xi,theta=theta,guess=structure$guess,tau=tau)
    } else {
      P<-ProbOgive(xi=xi,theta=theta,guess=structure$guess)
    }
  } else if (tolower(structure$icc)=="logistic") {
    if (!is.na(structure$ncat) && structure$ncat!=2) {
      P<-ProbIRT(xi,theta,guess=structure$guess)
    } else {
      P<-ProbIRT(xi,theta,guess=structure$guess,tau)
    }
  }  
  if (is.na(tau)) {
    for (i in 1:nrow(RP)) for (j in 1:ncol(RP)) RP[i,j]<-rbinom(1,size=1,prob=P[i,j])
  } else {
    for (i in 1:nrow(RP)) for (j in 1:ncol(RP)) {
      RP[i,j]<-which(rmultinom(1,size=1,prob=P[i,j,])==1)
    }
  }
  return(RP)
}
