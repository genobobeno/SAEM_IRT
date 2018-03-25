DrawA <- function(covZ,Q,a) {
  # DrawA generates A factor coefficients w/ eigenanalysis
  out 		<- eigen( covZ, symmetric=TRUE)
  if (!is.na(a)) {
    if (Q>1) {
      Aload<-as.matrix(apply(rbind(out$vectors[,1:Q],a),2,
                             function(x) (if (mean(x[1:(length(x)/2)]*x[length(x)/2+1:(length(x)/2)])>0 & sum(x[1:(length(x)/2)])>0) {return(x[1:(length(x)/2)])} 
                                          else {return(-1*x[1:(length(x)/2)])})),nrow(covZ),Q)
      Atemp = out$vectors[,1:Q]%*%sqrt(diag(out$values[1:Q]))
    } else {
      ifelse(mean(out$vectors[,1]*a)>0,Aload<-out$vectors[,1],Aload<-(-1)*out$vectors[,1])
      Atemp = as.matrix(Aload*sqrt(out$values[1]))
    }
  } else {
    if (Q>1) {
      Atemp = out$vectors[,1:Q]%*%sqrt(diag(out$values[1:Q]))
    } else {
      Atemp = out$vectors[,1]*sqrt(out$values[1])
    }
  }
  list(Atemp=Atemp,Avec=out$values)
}
