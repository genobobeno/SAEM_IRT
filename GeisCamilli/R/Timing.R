Timing<-function(start=NA) {
  if (is.na(start)) {Sys.time()} else {Sys.time()-start}
}