CheckBurn <-
function(PArray=PArray,Lag=Lag) {
  ind<-nrow(PArray)-(Lag-1):0
  PArray<-PArray[ind,]
  MN<-apply(PArray,2,sd)
  return(MN)
}
