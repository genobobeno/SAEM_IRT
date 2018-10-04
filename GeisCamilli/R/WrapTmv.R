WrapTmv <- function(i,A,Z,BTB_INV,b,dbltrunc) {
  if (dbltrunc) {
    rmvn(1,BTB_INV%*%(t(A)%*%(Z[i,])),BTB_INV) # formerly  Z + b
  } else {
    rmvn(1,BTB_INV%*%(t(A)%*%(Z[i,]+b)),BTB_INV) # formerly  Z + b
  }
}
