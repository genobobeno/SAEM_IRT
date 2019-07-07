ResponseCurves<-function(responses,scores,prows=5,pcols=3,correct=NA,j.legend=c(24,6),score.bins=NA,...) {
  # responses<-b.RP; scores<-apply(Response,1,mean); correct<-Correct
  # stopifnot(sum(is.numeric(data.matrix(responses)))==length(data.matrix(responses)))
  library(Hmisc)
  J=ncol(responses)
  N=nrow(responses)
  S<-sort(unique(scores))
  if (length(S)>(J+1)) {
    bins<-J+1
    S<-cut2(scores,g=bins)
  } else {
    bins<-length(S)
  }
  ncat<-max(apply(responses,2, function(x) length(unique(x[!is.na(x)]))))
  PCH<-c(0,1,2,6,3,7,9)[1:ncat]
#  par(mfrow=c(prows,pcols),...)
  par(mfrow=c(prows,pcols),cex=0.7)
  l.fits<-list()
  for (j in 1:J) {
    if ((j%%pcols)==1 & ceiling((((j-1)%%(prows*pcols))+1)/pcols)!=prows) {
      par(mar=c(1,4,1,1))
    } else if ((j%%pcols)!=1 & ceiling((((j-1)%%(prows*pcols))+1)/pcols)!=prows) {
      par(mar=c(1,1,1,1))
    } else if ((j%%pcols)==1 & ceiling((((j-1)%%(prows*pcols))+1)/pcols)==prows) {
      par(mar=c(4,4,1,1))
    } else {
      par(mar=c(4,1,1,1))
    }
    plot(range(scores),c(0,1),type="n",main=NA,
         ylab=ifelse((j%%pcols)==1,"percent",NA),
         yaxt="n",xaxt="n",
         xlab=ifelse(ceiling((((j-1)%%(prows*pcols))+1)/pcols)==prows,"score",NA))
    if ((j%%pcols)==1) axis(2)
    if (ceiling((((j-1)%%(prows*pcols))+1)/pcols)==prows) axis(1)
    rf<-factor(responses[,j],levels=1:ncat)
    if (is.numeric(S)) {
      p<-sapply(S,function(x) {
           if (is.na(rf[scores==x]) == length(rf[scores==x])) {
             rep(0,length(levels(rf)))
           } else {
             sapply(levels(rf),function(y) mean(rf[scores==x]==y,na.rm=TRUE))
           }
         })
      for (b in 1:nrow(p)) {
        points(x = S,y = p[b,],pch=PCH[b],col=c(1,2)[ifelse(!is.na(correct)[1],ifelse(b==correct[j],2,1),1)])
      }
      if (!is.na(correct)[1]) {
        fit2 <- glm(p[correct[j],] ~ S, family=binomial)
        lines(S, fit2$fitted, type="l", col="red")
        l.fits[[paste0("J",j)]]<-fit2
      }
      text(0.5,0.95,paste("Item",j))
    } else {
      p<-sapply(levels(S),function(x) {
        if (is.na(rf[S==x]) == length(rf[S==x])) {
          rep(0,length(levels(rf)))
        } else {
          sapply(levels(rf),function(y) mean(rf[S==x]==y,na.rm=TRUE))
        }
      })
      x.data = sapply(levels(S),function(x) mean(scores[S==x]))
      for (b in 1:nrow(p)) {
        points(x = x.data,
               y = p[b,],pch=PCH[b],col=c(1,2)[ifelse(!is.na(correct)[1],ifelse(b==correct[j],2,1),1)])
      }
      if (!is.na(correct)[1]) {
        fit2 <- glm(p[correct[j],] ~ x.data, family=binomial)
        lines(x.data, fit2$fitted, type="l", col="red")
        l.fits[[paste0("J",j)]]<-fit2
      }
      text(mean(scores),0.95,paste("Item",j))
    }
    if (j %in% j.legend) {
      legend("right",legend = LETTERS[1:ncat],pch=PCH)
    }
  }
  if (!is.na(correct)[1]) return(l.fits)
}
