CLASS<-function(testcut="P",PrePost="Pre",extracut=NA) {
  print("testcut = \"P\" for physics, 'C' for chemistry")
  print("PrePost = \"Pre\" for pretest, \"Post\" for posttest")

  if (testcut=="C") {
    N=50
    Acats<-outer(c("Over","All","PI","RW","Gen","Con","Soph","SE","CC","CL","Per"),c("F","U"),function(x,y) paste(x,y,sep=""))
  } else {
    N=42
    Acats<-outer(c("Over","All","PI","RW","Gen","Con","Soph","SE","CU","AC"),c("F","U"),function(x,y) paste(x,y,sep=""))
  }
  if (PrePost=="Post") {
    Qcols<-paste("q",tolower(testcut),1:N,sep="")
    Acols<-paste(tolower(testcut),"class",t(Acats),sep="")
  } else {
    Qcols<-paste("Q",testcut,1:N,sep="")
    Acols<-paste(testcut,"CLASS",t(Acats),sep="")
  }
  Data<-sapply(IDF[,Acols], function(x) x/100)
  Opin<-data.frame(o1=Data[,1]+Data[,2])
  Polr<-data.frame(p1=Data[,1]-Data[,2])
  Map<-data.frame(m1=log(Data[,1]/Data[,2]))
  names(Opin)[1]<-substr(Acols[1], 1, nchar(Acols[1])-1)
  names(Polr)[1]<-substr(Acols[1], 1, nchar(Acols[1])-1)
  names(Map)[1]<-substr(Acols[1], 1, nchar(Acols[1])-1)
  for (i in 2:(length(Acols)/2)) {
    Opin<-cbind(Opin,Data[,2*i-1]+Data[,2*i])    
    Polr<-cbind(Polr,Data[,2*i-1]-Data[,2*i])
    Map<-cbind(Map,log(Data[,2*i-1]/Data[,2*i]))
    names(Opin)[i]<-substr(Acols[2*i], 1, nchar(Acols[2*i])-1)
    names(Polr)[i]<-substr(Acols[2*i], 1, nchar(Acols[2*i])-1)
    names(Map)[i]<-substr(Acols[2*i], 1, nchar(Acols[2*i])-1) 
  }
  panel.hist <- function(x, ...)
  {
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(usr[1:2], 0, 1.5) )
    h <- hist(x, plot = FALSE)
    breaks <- h$breaks; nB <- length(breaks)
    y <- h$counts; y <- y/max(y)
    rect(breaks[-nB], 0, breaks[-1], y, col = "cyan", ...)
  }
  pairs(Opin, panel = panel.smooth,
        cex = 1, pch = 1, bg = "light blue",
        diag.panel = panel.hist, cex.labels = 1, font.labels = 2)
  pairs(Polr, panel = panel.smooth,
        cex = 1, pch = 1, bg = "light blue",
        diag.panel = panel.hist, cex.labels = 1, font.labels = 2)
  pairs(Map, panel = panel.smooth,
        cex = 1, pch = 1, bg = "light blue",
        diag.panel = panel.hist, cex.labels = 1, font.labels = 2)
  print("Opinionated?")
  for (i in 1:length(Opin)) {
    COR<-cor(cbind(IDF$Grade123,Opin[,i]),use="pairwise.complete.obs")      
    print(paste(colnames(Opin)[i],COR[1,2]))
  }
  print("Polarized?")
  for (i in 1:length(Polr)) {
    COR<-cor(cbind(IDF$Grade123,Polr[,i]),use="pairwise.complete.obs")      
    print(paste(colnames(Polr)[i],COR[1,2]))
  }
  print("Mapped onto 1D?")
  for (i in 1:length(Map)) {
    COR<-cor(cbind(IDF$Grade123,Map[,i]),use="pairwise.complete.obs")      
    print(paste(colnames(Map)[i],COR[1,2]))
  }
  
  plot(IDF[,Acols[2*i-1]]~IDF[,Acols[2*i]],main=paste(Acols[2*i-1],"vs.",Acols[2*i]),ylim=c(0,1),xlim=c(0,1))


}


CLASSCols<-function() {
print(" The following is a list of the Physics CLASS categories and ")
print("   their column NAMES and NUMBERS. The list is organized as")
print("   follows:")
print("")
print("             COLUMN CATEGORY       : PretestName  NUMBER: PosttestName NUMBER")
print("=============================================================================")
print("                  Overall - Fav    :   PCLASSOverF   403:   pclassOverF  804")
print("                  Overall - Unfav  :   PCLASSOverU   404:   pclassOverU  805")
print("                      All - Fav    :   PCLASSAllF    405:   pclassAllF   806")
print("                      All - Unfav  :   PCLASSAllU    406:   pclassAllU   807")
print("        Personal Interest - Fav    :   PCLASSPIF     407:   pclassPIF    808")
print("        Personal Interest - Unfav  :   PCLASSPIU     408:   pclassPIU    809")
print("    Real World Connection - Fav    :   PCLASSRWF     409:   pclassRWF    810")
print("    Real World Connection - Unfav  :   PCLASSRWU     410:   pclassRWU    811")
print("                  General - Fav    :   PCLASSGenF    411:   pclassGenF   812")
print("                  General - Unfav  :   PCLASSGenU    412:   pclassGenU   813")
print("               Confidence - Fav    :   PCLASSConF    413:   pclassConF   814")
print("               Confidence - Unfav  :   PCLASSConU    414:   pclassConU   815")
print("           Sophistication - Fav    :   PCLASSSophF   415:   pclassSophF  816")
print("           Sophistication - Unfav  :   PCLASSSophU   416:   pclassSophU  817")
print("      Sense Making/Effort - Fav    :   PCLASSSEF     417:   pclassSEF    818")
print("      Sense Making/Effort - Unfav  :   PCLASSSEU     418:   pclassSEU    819")
print(" Conceptual Understanding - Fav    :   PCLASSCUF     419:   pclassCUF    820")
print(" Conceptual Understanding - Unfav  :   PCLASSCUU     420:   pclassCUU    821")
print("    Applied Understanding - Fav    :   PCLASSACF     421:   pclassACF    822")
print("    Applied Understanding - Unfav  :   PCLASSACU     422:   pclassACU    823")
print(" ..............................................................................")
print("                 FOR CHEMISTRY, CHANGE THE 'P' to a 'C'   ")
print(" ..............................................................................")
print("")
print("      CHEMISTRY CLASS CATEGORY     :     PRETEST     :     POSTTEST    ")
print("==============================================================================")
print("                  Overall - Fav    :   CCLASSOverF   338:   cclassOverF  739")
print("                  Overall - Unfav  :   CCLASSOverU   339:   cclassOverU  740")
print("                      All - Fav    :   CCLASSAllF    340:   cclassAllF   741")
print("                      All - Unfav  :   CCLASSAllU    341:   cclassAllU   742")
print("        Personal Interest - Fav    :   CCLASSPIF     342:   cclassPIF    743")
print("        Personal Interest - Unfav  :   CCLASSPIU     343:   cclassPIU    744")
print("    Real World Connection - Fav    :   CCLASSRWF     344:   cclassRWF    745")
print("    Real World Connection - Unfav  :   CCLASSRWU     345:   cclassRWU    746")
print("                  General - Fav    :   CCLASSGenF    346:   cclassGenF   747")
print("                  General - Unfav  :   CCLASSGenU    347:   cclassGenU   748")
print("               Confidence - Fav    :   CCLASSConF    348:   cclassConF   749")
print("               Confidence - Unfav  :   CCLASSConU    349:   cclassConU   750")
print("           Sophistication - Fav    :   CCLASSSophF   350:   cclassSophF  751")
print("           Sophistication - Unfav  :   CCLASSSophU   351:   cclassSophU  752")
print("      Sense Making/Effort - Fav    :   CCLASSSEF     352:   cclassSEF    753")
print("      Sense Making/Effort - Unfav  :   CCLASSSEU     353:   cclassSEU    754")
print("   Conceptual Connections - Fav    :   CCLASSCCF     354:   cclassCCF    755")
print("   Conceptual Connections - Unfav  :   CCLASSCCU     355:   cclassCCU    756")
print("      Conceptual Learning - Fav    :   CCLASSCLF     356:   cclassCLF    757")
print("      Conceptual Learning - Unfav  :   CCLASSCLU     357:   cclassCLU    758")
print("Atom/Molecule Perspective - Fav    :   CCLASSPerF    358:   cclassPerF   759")
print("Atom/Molecule Perspective - Unfav  :   CCLASSPerU    359:   cclassPerU   760")
}


AnalyzeCLASS<-function(testcut="P",PrePost="Pre",extracut=NA) {
  
  rowProds <- function(X){ apply(X,1,FUN="prod") }
  
  print("testcut = 'P' for physics, 'C' for chemistry")
  print("PrePost = 'Pre' for pretest, 'Post' for posttest, 'Change' for Post-Pre")
  if (!is.na(extracut)) {
    IDF<-IDF[eval(parse(text=extracut)),] 
  }
  PrePost="Change"
  if (testcut=="C") {
    Cats<-outer(c("Over","All","PI","RW","Gen","Con","Soph","SE","CC","CL","Per"),c("F","U"),function(x,y) paste(x,y,sep=""))
  } else {
    Cats<-outer(c("Over","All","PI","RW","Gen","Con","Soph","SE","CU","AC"),c("F","U"),function(x,y) paste(x,y,sep=""))
  }
  PFname<-paste(testcut,"CLASS",sep="")
  if (PrePost=="Pre") { 
    CCols<-which(colnames(IDF) %in% paste(PFname,t(Cats),sep=""))
  } else if (PrePost=="Post") { 
    CCols<-which(colnames(IDF) %in% paste(tolower(PFname),t(Cats),sep=""))
  } else {
    CCols<-c(which(colnames(IDF) %in% paste(PFname,t(Cats),sep="")),which(colnames(IDF) %in% paste(tolower(PFname),t(Cats),sep="")))
  }

  IDF[,CCols]<-IDF[,CCols]/100
  
  if (PrePost=="Pre"|PrePost=="Post") {
    par(mfrow=c(3,4))
    for (i in 1:(length(CCols)/2)) {
      plot(IDF[,CCols[2*i]],IDF[,CCols[2*i-1]],ylim=c(0,1),xlim=c(0,1),main=paste(PrePost,substr(Cats[i,1],1,nchar(Cats[i,1])-1),"Fav vs. Unfav",sep=" : "))
      abline(a=1,b=-1)
    }
    par(mfrow=c(3,4))
    for (i in 1:(length(CCols)/2)) {
      CONTOUR(IDF[,CCols[2*i]],IDF[,CCols[2*i-1]],25,paste(PrePost,substr(Cats[i,1],1,nchar(Cats[i,1])-1),"Fav vs. Unfav",sep=" : "))
      abline(a=1,b=-1)
    }
  } else {
    MCols<-matrix(CCols,nrow=length(CCols)/2,ncol=2)
    par(mfrow=c(3,4))
    for (i in 1:(nrow(MCols)/2)) {
      hist(IDF[,MCols[2*i-1,2]]-IDF[,MCols[2*i-1,1]],breaks=8,main=paste(PrePost,substr(Cats[i,1],1,nchar(Cats[i,1])-1),"Fav Difference",sep=" : "))
      abline(v=0)
      print(paste(substr(Cats[i,1],1,nchar(Cats[i,1])-1),"Mean Fav Difference",format(mean(IDF[,MCols[2*i-1,2]]-IDF[,MCols[2*i-1,1]],na.rm=T),digits=4),"SD",format(sd(IDF[,MCols[2*i-1,2]]-IDF[,MCols[2*i-1,1]],na.rm=T),digits=4),sep=" : "))
    }
    par(mfrow=c(3,4))
    for (i in 1:(nrow(MCols)/2)) {
      hist(IDF[,MCols[2*i,2]]-IDF[,MCols[2*i,1]],breaks=8,main=paste(PrePost,substr(Cats[i,1],1,nchar(Cats[i,1])-1),"Unfav Difference",sep=" : "))
      abline(v=0)
      print(paste(substr(Cats[i,1],1,nchar(Cats[i,1])-1),"Mean Unfav Difference",format(mean(IDF[,MCols[2*i,2]]-IDF[,MCols[2*i,1]],na.rm=T),digits=4),"SD",format(sd(IDF[,MCols[2*i,2]]-IDF[,MCols[2*i,1]],na.rm=T),digits=4),sep=" : "))
    }
    par(mfrow=c(3,4))
    for (i in 1:(nrow(MCols)/2)) {
      plot(IDF[,MCols[2*i,2]]-IDF[,MCols[2*i,1]],IDF[,MCols[2*i-1,2]]-IDF[,MCols[2*i-1,1]],main=paste(PrePost,substr(Cats[i,1],1,nchar(Cats[i,1])-1),"Fav vs. Unfav Difference",sep=" : "))
      abline(h=0,v=0)
    }
#     par(mfrow=c(3,4))
#     CONTOUR<-function(x=x,y=y,nbins=nbins,title=NA) {
#       xy <- cbind(x,y)
#       x.bin <- seq(floor(min(xy[,1],na.rm=T)), ceiling(max(xy[,1],na.rm=T)), length=nbins)
#       y.bin <- seq(floor(min(xy[,2],na.rm=T)), ceiling(max(xy[,2],na.rm=T)), length=nbins)
#       
#       freq <-  as.data.frame(table(findInterval(xy[,1], x.bin),findInterval(xy[,2], y.bin)))
#       freq[,1] <- as.numeric(freq[,1])
#       freq[,2] <- as.numeric(freq[,2])
#       
#       freq2D <- diag(nbins)*0
#       freq2D[cbind(freq[,1], freq[,2])] <- freq[,3]
#       
#       image(x.bin, y.bin, freq2D, col=topo.colors(max(freq2D)),main=title)
#       contour(x.bin, y.bin, freq2D, add=TRUE, col=rgb(1,1,1,.7))
#     }
#     for (i in 1:(nrow(MCols)/2)) {
#       CONTOUR(IDF[,MCols[2*i,2]]-IDF[,MCols[2*i,1]],IDF[,MCols[2*i-1,2]]-IDF[,MCols[2*i-1,1]],25,paste(PrePost,substr(Cats[i,1],1,nchar(Cats[i,1])-1),"Fav vs. Unfav Difference",sep=" : "))
#       abline(h=0,v=0)
#     }
    par(mfrow=c(3,4))
    for (i in 1:(nrow(MCols)/2)) {
      plot(IDF[,MCols[2*i,1]],IDF[,MCols[2*i-1,1]],main=paste(PrePost,substr(Cats[i,1],1,nchar(Cats[i,1])-1),"Fav vs. Unfav Changes",sep=" : "),type="n",ylim=c(0,1),xlim=c(0,1))
      ar<-cbind(IDF[,MCols[2*i,1]],IDF[,MCols[2*i-1,1]],IDF[,MCols[2*i,2]],IDF[,MCols[2*i-1,2]])      
      ar<-ar[which(rowProds(ar)>0),]
      abline(a=1,b=-1,col=2)
      arrows(ar[,1],ar[,2],ar[,3],ar[,4],length=0.05)
    }
    par(mfrow=c(3,4))
    for (i in 1:(nrow(MCols)/2)) {
      ar<-cbind(IDF[,MCols[2*i,1]],IDF[,MCols[2*i-1,1]],IDF[,MCols[2*i,2]],IDF[,MCols[2*i-1,2]])      
      ar<-ar[which(rowProds(ar)>0),]
      delta<-(ar[,3]-ar[,1])*1+(ar[,4]-ar[,2])*-1
      hist(delta,main=paste(substr(Cats[i,1],1,nchar(Cats[i,1])-1),"Projected Change toward Favorable"))
      print(paste(substr(Cats[i,1],1,nchar(Cats[i,1])-1),"Mean",format(mean(delta),digits=4),"SD",format(sd(delta),digits=4),sep=" : "))
    }
  }
}

  