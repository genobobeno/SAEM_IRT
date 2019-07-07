ResponseCurves<-function(responses,scores,prows=5,pcols=3,correct=NA,j.legend=c(24,6),score.bins=NA,missing=9,...) {
  # responses<-b.RP; scores<-apply(Response,1,mean); correct<-Correct
  # stopifnot(sum(is.numeric(data.matrix(responses)))==length(data.matrix(responses)))
  #responses = R.QOL-2;scores = rowMeans(R.QOL-2,na.rm = TRUE)
  #prows = 6;pcols = 4;correct = rep(2,ncol(R.QOL));j.legend = c(6,18);score.bins = 25
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
  if (sum(apply(responses,2, function(x) sum(x==missing,na.rm=TRUE)))>0) {
    responses[responses==missing]<-NA
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
    if (ncat==2) {
      rf<-factor(responses[,j],levels=1:ncat)
    } else {
      Lvls<-unique(responses[!is.na(responses[,j]),j])
      rf<-factor(responses[,j],levels=Lvls)
    }
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
      text(mean(range(scores,na.rm=TRUE)),0.95,paste("Item",j))
    } else {
      p<-sapply(levels(S),function(x) {
        if (is.na(rf[S==x]) == length(rf[S==x])) {
          rep(0,length(levels(rf)))
        } else {
          sapply(levels(rf),function(y) mean(rf[S==x]==y,na.rm=TRUE))
        }
      })
      x.data = sapply(levels(S),function(x) mean(scores[S==x]))
      if (ncat>2 && sum(grepl("(-)?[0-9]+",row.names(p)))==nrow(p)) {
        p<-p[as.character(sort(as.numeric(row.names(p)))),]
        for (b in 1:nrow(p)) {
          points(x = x.data,
                 y = p[b,],pch=PCH[b],col=c(1,2)[ifelse(!is.na(correct)[1],ifelse(row.names(p)[b]==correct[j],2,1),1)])
        }
      } else {
        for (b in 1:nrow(p)) {
          points(x = x.data,
                 y = p[b,],pch=PCH[b],col=c(1,2)[ifelse(!is.na(correct)[1],ifelse(b==correct[j],2,1),1)])
        }
      }
      if (!is.na(correct)[1]) {
        if (ncat>2 && sum(grepl("(-)?[0-9]+",row.names(p)))==nrow(p)) {
          fit2 <- glm(p[as.character(correct[j]),] ~ x.data, family=binomial)
        } else {
          fit2 <- glm(p[correct[j],] ~ x.data, family=binomial)
        }
        lines(x.data, fit2$fitted, type="l", col="red")
        l.fits[[paste0("J",j)]]<-fit2
      }
      text(mean(range(scores)),0.95,paste("Item",j))
    }
    if (j %in% j.legend) {
      if (ncat==2) {
        legend("right",legend = LETTERS[1:ncat],pch=PCH)
      } else if (sum(grepl("(-)?[0-9]+",row.names(p)))==nrow(p)) {
        legend("right",legend = sort(as.numeric(row.names(p))),pch=PCH)
      } else {
        legend("right",legend = Lvls,pch=PCH)
      }
    }
  }
  if (!is.na(correct)[1]) return(l.fits)
}
