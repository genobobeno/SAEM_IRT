WrapT <- function(i,A,Z,BTB_INV,b) { 
  rmvn(1,BTB_INV%*%(t(A)%*%(Z[i,] + b)),BTB_INV)
}
