AutocovMC<-function(x,end=800,start=600,items=1:5,title.append=NA) {
  #par(mfrow=c(5,2),mar=c(3,3,3,1))
  if (is.na(items)[1]) items<-dim(x$Aiter)[1]
  for (j in items) {
    acf(x$Aiter[j,1,start:end],main=paste0(ifelse(is.na(title.append),"",title.append),"Slope, Item ",j),xlab=NA)
    acf(x$Biter[j,start:end],main=paste0(ifelse(is.na(title.append),"",title.append),"Intercept, Item ",j),xlab=NA)
  }
}
