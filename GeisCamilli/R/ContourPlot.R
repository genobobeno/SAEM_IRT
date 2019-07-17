ContourPlot <- function(df = NA, var1 = NA, var2 = NA, bins=100, ...) {
  library(Hmisc)
  cat("\n\nThis function takes either,
      1. A two-column data.frame: ContourPlot(df=DF), or
      2. Two vectors of the same length: ContourPlot(var1 = V1, var2 = V2)\n")
  stopifnot(is.na(df) | (ncol(df)==2 & is.na(var1) & is.na(var2)))
  if (!is.na(df)) {
    var1 = df[,1] #; XN = names(df)[1]
    var2 = df[,2] #; YN = names(df)[2]
  } 
  # else {
  #   XN = "Var1"; YN = "Var2"
  # }
  if (is.integer(var1)) {
    if (length(unique(var1))>bins-1 & length(unique(var1))<2*bins) {
      v1c=c(min(var1,na.rm=TRUE)+0:(length(unique(var1))/2)*diff(range(var1,na.rm=TRUE))/(length(unique(var1))/2))
    } else if (length(unique(var1))>2*bins) {
      v1c=c(min(var1,na.rm=TRUE)+0:bins*diff(range(var1,na.rm=TRUE))/bins)
    } else if (length(unique(var1))<21) {
      v1c=c(min(var1,na.rm=TRUE):(max(var1,na.rm=TRUE)+1))
    } else {
      v1c=c(min(var1,na.rm=TRUE):max(var1,na.rm=TRUE))
    }
  } else if (is.character(var1)) {
    var1 = factor(var1)
    v1c = levels(var1)
  } else {
    v1c=c(min(var1,na.rm=TRUE)-0.0001+0:50*diff(range(var1,na.rm=TRUE))/50+c(rep(0,50),0.0001))
  }
  if (is.integer(var2)) {
    if (length(unique(var2))>bins-1 & length(unique(var2))<2*bins) {
      v2c=c(min(var2,na.rm=TRUE)+0:(length(unique(var2))/2)*diff(range(var2,na.rm=TRUE))/(length(unique(var2))/2))
    } else if (length(unique(var2))>2*bins) {
      v2c=c(min(var2,na.rm=TRUE)+0:bins*diff(range(var2,na.rm=TRUE))/bins)
    } else if (length(unique(var2))<21) {
      v2c=c(min(var2,na.rm=TRUE):(max(var2,na.rm=TRUE)+1))
    } else {
      v2c=c(min(var2,na.rm=TRUE):max(var2,na.rm=TRUE))
    }
  } else if (is.character(var2)) {
    var2 = factor(var2)
    v2c = levels(var2)
  } else {
    v2c=c(min(var2,na.rm=TRUE)-0.0001+0:50*diff(range(var2,na.rm=TRUE))/50+c(rep(0,50),0.0001))
  }
  v1n = unique(as.character(cut2(seq(mean(v1c[1:2]),mean(v1c[length(v1c)-1:0]),length.out = 1000),cuts=signif(v1c,6))))
  v2n = unique(as.character(cut2(seq(mean(v2c[1:2]),mean(v2c[length(v2c)-1:0]),length.out = 1000),cuts=signif(v2c,6))))
  v1f = cut2(var1,cuts=signif(v1c,6))
  v2f = cut2(var2,cuts=signif(v2c,6))
  vtab<-table(v1f,v2f,useNA = "ifany")/length(var1)
  if (length(colnames(vtab))>(length(v2c)-1)) {
    vtab<-vtab[,1:(length(v2c)-1)]
  }
  image(v1c,v2c,log(vtab),col=rainbow(1000,start = 0.2,end = 0.75),...)
}
