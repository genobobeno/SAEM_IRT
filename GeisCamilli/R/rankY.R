rankY<-function(v) { #returns rank array; opposite of order()
  O<-cbind(1:length(v),order(v))
  O[order(O[,2]),1]
}