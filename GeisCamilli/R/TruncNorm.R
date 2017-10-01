TruncNorm <-
function(Means,HighLow,prl) {
  #   Ez[which(w[,j]==1),j]<-rtruncnorm( length(Ez[which(w[,j]==1),j]), a=0,  b=Inf, mean=Eta[which(w[,j]==1),j], sd=1)
  #   Ez[which(w[,j]==0),j]<-rtruncnorm( length(Ez[which(w[,j]==0),j]), a=-Inf, b=0, mean=Eta[which(w[,j]==0),j], sd=1)
  N<-length(Means)
  if (prl) {
    if (HighLow=="high") {
      Q<-unlist(mclapply(Means,function(x) {if(x>4) {LB=0.0000001+pnorm(0-x)} else {LB=pnorm(0-x)}
                                          return(x+qnorm(runif(1,LB,1)-0.00000001))}))
      #      x <- qnorm( runif(1, pnorm(0-Means[i]), 1) )
    } else {
      Q<-unlist(mclapply(Means,function(x) {if(x<(-4)) {UB=pnorm(0-x)-0.0000001} else {UB=pnorm(0-x)}
                                          return(x+qnorm(runif(1,0,UB)+0.00000001))}))
    }
  } else {
    if (HighLow=="high") {
      Q<-unlist(lapply(Means,function(x) {if(x>4) {LB=0.0000001+pnorm(0-x)} else {LB=pnorm(0-x)}
                                            return(x+qnorm(runif(1,LB,1)-0.00000001))}))
      #      x <- qnorm( runif(1, pnorm(0-Means[i]), 1) )
    } else {
      Q<-unlist(lapply(Means,function(x) {if(x<(-4)) {UB=pnorm(0-x)-0.0000001} else {UB=pnorm(0-x)}
                                            return(x+qnorm(runif(1,0,UB)+0.00000001))}))
    }    
  }
  return(Q)
}
