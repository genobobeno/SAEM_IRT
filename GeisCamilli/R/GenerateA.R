GenerateA <-
function(j,Adim,Adist,Aparams) {
  library(truncnorm)
  if (Adist=="unif") {
    stopifnot(length(Aparams)==2) # Aparams are the limits to uniform
    A<-runif(j,min=Aparams[1],max=Aparams[2]) #1D, N=5000, J=50
    if (Adim>1) {
      for (i in 2:Adim) {
        A<-cbind(A,runif(j,min=Aparams[1]-0.1,max=1))
        if (i==2) A[1:((Adim-1)*ceiling(j/Adim)),i]<-0
        if (i>=3) {
          A[1:((Adim-(i-1))*ceiling(j/Adim)),i]<-0
          A[(((Adim-(i-2))*ceiling(j/Adim))+1):j,i]<-0    
        }
      }
    }
  } else if (Adist=="norm") {
    stopifnot(length(Aparams)==2)
    A<-rtruncnorm(j,a=0.2,b=2.0,Aparams[1],Aparams[2]) #1D, N=5000, J=50
    if (Adim>1) {
      for (i in 2:Adim) {
        A<-cbind(A,rtruncnorm(j,a=0.2,b=2.0,Aparams[1],Aparams[2]))
        if (i==2) A[1:((Adim-1)*ceiling(j/Adim)),i]<-0
        if (i>=3) {
          A[1:((Adim-(i-1))*ceiling(j/Adim)),i]<-0
          A[(((Adim-(i-2))*ceiling(j/Adim))+1):j,i]<-0    
        }
      }
    }
  } else if (Adist=="beta") {
    stopifnot(length(Aparams)==2,Aparams[1]>0,Aparams[2]<2.5)
    A<-rbeta(j,2.5,3.0) #1D, N=5000, J=50
    A<-(Aparams[2]-Aparams[1])*A+Aparams[1]
    if (Adim>1) {
      for (i in 2:Adim) {
        A<-cbind(A,(1.0-0.1)*rbeta(j,2.5,3.0)+0.1)
        if (i==2) A[1:((Adim-1)*ceiling(j/Adim)),i]<-0
        if (i>=3) {
          A[1:((Adim-(i-1))*ceiling(j/Adim)),i]<-0
          A[(((Adim-(i-2))*ceiling(j/Adim))+1):j,i]<-0    
        }
      }
    }
  }  
  return(A)
}
