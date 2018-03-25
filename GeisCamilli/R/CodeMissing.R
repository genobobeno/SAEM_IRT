CodeMissing<-function(rp,settings) {
  K=max(apply(rp,2,function(x) (length(unique(x[!is.na(x) | x!=settings$missing])))))
  if (!is.numeric(rp)) {
    NJ<-dim(rp)
    if (is.na(settings$missing)) {
      rp[is.na(rp)]<-'9'
    } else {
      if (sum(is.na(rp))>0) rp[is.na(rp)]<-'9'
      if (sum(rp[!rp==settings$missing] %in% c(letters,LETTERS))>0) {
        cat(paste("\nFunctions are expecting codings of 0 thru",K-1,"as responses in the data.
        You should recode the data or write a function to transform it.
        Things will probably break now...\n\n"))
      }
      rp[rp==settings$missing]<-'9'
    }
    rp<-matrix(as.numeric(rp),nrow=NJ[1],ncol=NJ[2])
  } else {
    if (is.na(settings$missing)) {
      rp[is.na(rp)]<-9
    } else {
      if (sum(is.na(rp))>0) rp[is.na(rp)]<-9
      rp[rp==settings$missing]<-9
    }
  }
  rp
}