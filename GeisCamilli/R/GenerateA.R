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
  } else if (Adist=="bifactor") {
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
  } else if (Adist=="bifactor2") {
    stopifnot(length(Aparams)==2,Aparams[1]>0,Aparams[2]<2.5)
    A<-rbeta(j,2.5,3.0) #1D, N=5000, J=50
    A<-(Aparams[2]-Aparams[1])*A+Aparams[1]
    if (Adim>1) {
      j.bi<-1:floor(j/(Adim-1))
      for (i in 2:Adim) {
        bi.TF<-1:j==(i-2)*floor(j/(Adim-1))+j.bi
        A<-cbind(A,(1.0-0.1)*rbeta(j,2.5,3.0)+0.1)
        A[!bi.TF,i]<-0
      }
    }
  } else if (Adist=="subscale") {
    if (Adim>1) {
      A<-mat.or.vec(j,Adim)
      for (i in 1:Adim) {
        A[((i-1)*ceiling(j/Adim)+1):(i*ceiling(j/Adim)),i]<-rbeta(ceiling(j/Adim),2.5,3.0) #1D, N=5000, J=50
        A[which(A[,i]!=0),i]<-(Aparams[2]-Aparams[1])*A[which(A[,i]!=0),i]+Aparams[1]
      }
    } else {
      A<-rbeta(j,2.5,3.0) #1D, N=5000, J=50
      A<-(Aparams[2]-Aparams[1])*A+Aparams[1]
    }
  }
  return(A)
}
