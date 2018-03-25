DrawAEigen <- function(covZ,Q)
{
  # Subtract diag(m)
  out 		<- eigen( covZ, symmetric=TRUE)
  list(Atemp=out$vectors[,1:Q]%*%sqrt(diag(out$values[1:Q])),Avec=out$values)
}