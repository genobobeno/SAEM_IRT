DrawAEigen <- function(covZ,Q)
{
  # Subtract diag(m)
  out 		<- eigen( covZ, symmetric=TRUE)
  if (sum(out$vectors[,1])<0) out$vectors[,1] <- -1*out$vectors[,1]
  if (Q==1) {
    list(Atemp=out$vectors[,1]*sqrt(out$values[1]),Avec=out$values) 
  } else {
    list(Atemp=out$vectors[,1:Q]%*%sqrt(diag(out$values[1:Q])),Avec=out$values)
  }
}