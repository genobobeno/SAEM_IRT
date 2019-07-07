vsplit<-function(x,f) {
  setNames(lapply(f,function (y) x[ceiling((1:length(x))/(length(x)/max(f)))==y]),f)
}
