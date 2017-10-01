NumHessian <-
function(j,xi,rp,pL,settings) {
  XI<-xi
  dp<-0.001
  if (!settings$guess) {
    XI[,1]<-xi[,1]-dp
    lowA<-Gradient(j,XI,rp,pL,settings)
    XI[,1]<-xi[,1]+dp
    highA<-Gradient(j,XI,rp,pL,settings)
    XI<-xi
    XI[,2]<-xi[,2]-dp
    lowB<-Gradient(j,XI,rp,pL,settings)
    XI[,2]<-xi[,2]+dp
    highB<-Gradient(j,XI,rp,pL,settings)
    ddA<-0.5*(highA-lowA)/dp
    ddB<-0.5*(highB-lowB)/dp
    ddL<-rbind(ddA,ddB)
  } else {
    XI[,1]<-xi[,1]-dp
    lowA<-Gradient(j,XI,rp,pL,settings)
    XI[,1]<-xi[,1]+dp
    highA<-Gradient(j,XI,rp,pL,settings)
    XI<-xi
    XI[,2]<-xi[,2]-dp
    lowB<-Gradient(j,XI,rp,pL,settings)
    XI[,2]<-xi[,2]+dp
    highB<-Gradient(j,XI,rp,pL,settings)
    XI<-xi
    XI[,3]<-xi[,3]-dp
    lowC<-Gradient(j,XI,rp,pL,settings)
    XI[,3]<-xi[,3]-dp
    highC<-Gradient(j,XI,rp,pL,settings)
    ddA<-0.5*(highA-lowA)/dp
    ddB<-0.5*(highB-lowB)/dp
    ddC<-0.5*(highC-lowC)/dp
    ddL<-rbind(ddA,rbind(ddB,ddC))    
  }
  return (ddL)
}
