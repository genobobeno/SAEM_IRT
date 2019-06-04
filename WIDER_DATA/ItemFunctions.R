#################################
# This is the code for the functions
#    1. ItemStats()
#    2. PrePostItem()
#    3. XBubble()

ItemStats <- function(it,yearcut=NA,schoolcut=NA,test="Pre") {
  #   it=78
  #   test="Pre"
  #   classcut=NA
  #   yearcut=NA
  #   schoolcut=NA
  #   formcut=NA
  #   # 4-Year  5-School  6-Class   8-Form  
  stopifnot(test=="Pre"|test=="Post"|is.numeric(it))
  if (test=="Pre") {
    Istr = c("C","B")
    Test<-TestX
  }
  if (test=="Post") {
    Istr = c("c","b")
    Test<-TestY
  }
  ## Put in a stopifnot(classcut %in% Class, yearcut %in% Year, etc... )
  # it = 40
  print("Item indexed in Master List\n")
  descript = paste("Item",it)
  print(paste(descript,Test[it,"Item"],sep=" : "))
  # Name of columns we want to analyze
  DF<-IDF
  CutIndex <- rep(0,4)
  CutString = c("Year = ","School = ","Class = ")
  ci = paste(Istr[1],it,sep="")
  bi = paste(Istr[2],it,sep="")
  
  if (!is.na(yearcut)) {
    DF<-DF[which(DF[,"Year"]==yearcut),]    
    CutString[1]<-paste(CutString[1],yearcut,sep="")
    CutIndex[1]<-1
  }
  if (!is.na(schoolcut)) {
    DF<-DF[which(DF[,"School"]==schoolcut),]    
    CutString[2]<-paste(CutString[2],schoolcut,sep="")
    CutIndex[2]<-2
  }
  
  #   if (!is.na(testcut)) {
  #     if (!is.character(test)) print("Try your cut again but put parentheses around your test name")
  #     DF<-DF[which(DF[,"Test"]==testcut),]    
  #     CutString[4]<-paste(CutString[4],testcut,sep="")
  #     CutIndex[4]<-4
  #   }
  #   if (!is.na(formcut)) {
  #     if (is.na(testcut)) print("You cut on the form but not on the Test... Numbers will probably not be interpretable")
  #     DF<-DF[which(DF[,"Form"]==formcut),]    
  #     CutString[4]<-paste(CutString[4],formcut,sep="")
  #     CutIndex[4]<-4
  #   }
  
  if (length(DF[,1])<5) print("Cuts are too strict or poorly defined, sample size is less than five")
  if (test=="Pre") {
    print(ci)
    COLB <- DF[which(!is.na(DF[,ci])),c(ci,bi,"Year","School","Score.x")]
    COLB[which(is.na(COLB[,2])),2]<-0
    for(i in 3:4)  COLB[,i]<-factor(COLB[,i])
    COLB[,2]<-factor(COLB[,2],levels=c("0","1","2","3","4","5"))
    COLN <- DF[which(!is.na(DF[,ci])&!is.na(DF[,bi])),c(ci,bi,"Year","School","Score.x")]
  }
  else {
    COLB <- DF[which(!is.na(DF[,ci])),c(ci,bi,"Year","School","Score.y")]
    COLB[which(is.na(COLB[,2])),2]<-0
    for(i in 3:4)  COLB[,i]<-factor(COLB[,i])
    COLB[,2]<-factor(COLB[,2],levels=c("0","1","2","3","4","5"))
    COLN <- DF[which(!is.na(DF[,ci])&!is.na(DF[,bi])),c(ci,bi,"Year","School","Score.y")]    
  }
  colnames(COLB)[5]<-"Score"
  colnames(COLN)[5]<-"Score"
  remove(DF)
  pCB = sum(COLB[,1])/length(COLB[,1])
  pCN = sum(COLN[,1])/length(COLN[,1])
  sdB = sqrt(pCB*(1-pCB)/length(COLB[,1]))
  sdN = sqrt(pCN*(1-pCN)/length(COLN[,1]))
  nB = length(COLB[,1])-length(COLN[,1])
  Bubs = c("A","B","C","D","E")
  CB = as.numeric(unique(COLB[which(COLB[,1]==1),2]))-1
  pBA = length(which(COLB[,2]==1))/length(COLB[,2])  
  pBB = length(which(COLB[,2]==2))/length(COLB[,2])  
  pBC = length(which(COLB[,2]==3))/length(COLB[,2])  
  pBD = length(which(COLB[,2]==4))/length(COLB[,2])  
  pBE = length(which(COLB[,2]==5))/length(COLB[,2])  
  
  print(paste("        Sample Size (Blanks included) : ",length(COLB[,1])))
  print(paste("        Students left this item blank : ",nB))
  print(paste("    Percent Correct (Blanks included) : ",format(100*pCB,digits=2,nsmall=2)))
  print(paste("       Binomial -  Standard Deviation : ",format(100*sdB,digits=2,nsmall=2)))
  print(paste("Percent Correct (Blanks Not included) : ",format(100*pCN,digits=2,nsmall=2)))
  print(paste("       Binomial -  Standard Deviation : ",format(100*sdN,digits=2,nsmall=2)))
  print(paste("      Percent of students who chose A : ",format(100*pBA,digits=2,nsmall=2)))
  print(paste("      Percent of students who chose B : ",format(100*pBB,digits=2,nsmall=2)))
  print(paste("      Percent of students who chose C : ",format(100*pBC,digits=2,nsmall=2)))
  print(paste("      Percent of students who chose D : ",format(100*pBD,digits=2,nsmall=2)))
  print(paste("      Percent of students who chose E : ",format(100*pBE,digits=2,nsmall=2)))
  
  par(mfrow=c(2,2))
  tmp<-length(COLB[,2])
  title = paste(Test[it,"ItemIndex"],paste(Test[it,"Item"],": Bubbled Responses\n"),sep="-")
  sampsize<-paste("N =",tmp)
  correct = paste("Correct Response is",paste(Bubs[CB],sampsize,sep="; "))
  title = paste(title,correct)
  plot(COLB[,2],xlab = "Student Answer", ylab = "Responses (%)", main = title, xaxt="n",yaxt="n",ylim=c(0,tmp))
  # draw an axis on the left
  axis(1, at=0.7+1.2*0:5, labels=c("Blank","A","B","C","D","E"))
  axis(2, at=c(0,tmp/2,tmp),labels=c(0,0.5,1.0))
  PropCor<-paste(format(100*pCB,digits=2,nsmall=2),"% Correct",sep="")
  title = paste(paste(paste(Test[it,"ItemIndex"],Test[it,"Item"],sep="-"),"Proportion correct\n",sep=" : "),PropCor)
  barplot(c(pCN,pCB), ylab = "Percent Correct", main = title, xaxt="n",ylim=c(0,1.0))
  axis(1, at=0.7+1.2*0:1, labels=c("No Blanks","Blank = Wrong"))
  COLB$ZS<-qqnorm(COLB$Score,main="Q-Q test for Normality of Scores")$x
  qqline(COLB$Score)
  
  breaks1=unique(quantile(COLB$ZS,c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0),na.rm=TRUE))
  breaks2=unique(quantile(COLB$Score,c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0),na.rm=TRUE))
  COLB$ZBreaks<-cut(COLB$ZS,breaks=breaks1)
  COLB$Breaks<-cut(COLB$Score,breaks=breaks2)
  
  #   sums<-by(COLB[,1], COLB$Breaks, sum)
  #   lengths<-by(COLB[,1], COLB$Breaks, length)
  #   plot(sums/lengths~c(-0.05+0.1*1:10),xlab="Quantiles (Z)")
  #   fit <- glm(cbind(sums, lengths-sums) ~ c(-0.05+0.1*1:10), family=binomial)
  #   lines(c(-0.05+0.1*1:10), fit$fitted, type="l", col="red")
  
  sums2<-by(COLB[,1], COLB$ZBreaks, sum)
  lengths2<-by(COLB[,1], COLB$ZBreaks, length)
  plot(sums2/lengths2~qnorm(-(0.5/length(sums2))+(1.0/length(sums2))*1:length(sums2)), ylim=c(0,1.0), main = "Logistic by Score", xlab="Standardized Quantiles by Score", ylab="Percent Correct")
  fit2 <- glm(cbind(sums2, lengths2-sums2) ~ c(qnorm(-(0.5/length(sums2))+(1.0/length(sums2))*1:length(sums2))), family=binomial)
  lines(qnorm(c(-0.05+0.1*1:10)), fit2$fitted, type="l", col="red")
  
}

XBubble <-function(it,test="Pre",yearcut=NA,schoolcut=NA) {
  stopifnot(test=="Pre"|test=="Post")
  #it=76
  if (test=="Pre") {
    Istr = c("C","B")
    Test<-TestX
  }
  if (test=="Post") {
    Istr = c("c","b")
    Test<-TestY
  }
  DF<-IDF
  CutIndex <- rep(0,4)
  CutString = c("Year = ","School = ","Class = ")
  ci = paste(Istr[1],it,sep="")
  bi = paste(Istr[2],it,sep="")
  
  Bubs = c("A","B","C","D","E")
  CB = as.numeric(unique(DF[which(DF[,ci]==1),bi]))
  print(CB)
  print(paste("Correct answer for Item",Test[it,"Item"],"is",Bubs[CB]))
  
  if (!is.na(yearcut)) {
    DF<-DF[which(DF[,"Year"]==yearcut),]    
    CutString[1]<-paste(CutString[1],yearcut,sep="")
    CutIndex[1]<-1
  }
  if (!is.na(schoolcut)) {
    DF<-DF[which(DF[,"School"]==schoolcut),]    
    CutString[2]<-paste(CutString[2],schoolcut,sep="")
    CutIndex[2]<-2
  }
  
  #   if (!is.na(testcut)) {
  #     if (!is.character(test)) print("Try your cut again but put parentheses around your test name")
  #     DF<-DF[which(DF[,"Test"]==testcut),]    
  #     CutString[4]<-paste(CutString[4],testcut,sep="")
  #     CutIndex[4]<-4
  #   }
  #   if (!is.na(formcut)) {
  #     if (is.na(testcut)) print("You cut on the form but not on the Test... Numbers will probably not be interpretable")
  #     DF<-DF[which(DF[,"Form"]==formcut),]    
  #     CutString[4]<-paste(CutString[4],formcut,sep="")
  #     CutIndex[4]<-4
  #   }
  if (length(DF[,1])<5) print("Cuts are too strict or poorly defined, sample size is less than five")
  dfXItem<-data.frame(Item=rep(it,5),Bubble=Bubs)
  if (it>1 & it<nrow(Test)) XItem<-c(1:(it-1),(it+1):nrow(Test))
  if (it==1) XItem<-c(2:nrow(Test))
  if (it==nrow(Test)) XItem<-c(1:nrow(Test))
  
  hitD<-rep(0,5*5*(ncol(Test)-1))
  hitN<-rep(0,ncol(Test)-1)
  for (j in 1:(nrow(Test)-1)) {
    hits<-table(DF[,bi],DF[,paste(Istr[2],XItem[j],sep="")])
    if (j==80) print(hits)
    hitN[j]<-sum(hits,na.rm=TRUE)
    hits<-as.matrix(hits)/sum(hits,na.rm=TRUE)
    if (length(hits)==25) hitD[25*j+1:25]<-as.vector(hits)
    for (i in 0:4) {
      dfXItem<-cbind(dfXItem,hits[1:5+i*5])
      colnames(dfXItem)[ncol(dfXItem)]<-paste("Q",j,"_",Bubs[i+1],sep="")
    }
  }
  par(mfrow=c(2,2))
  hist(hitD[which(hitD>0.0)],main="paired hits on distractors")
  plot(hitN~XItem, main=paste("Item",it))
  #   Items = 3
  #   Qs<-paste("Q",XItem[order(-hitN)[1:10]],sep="")
  #   As<-paste("_",Bubs,sep="")
  #   PolyA<-outer(Qs,As,FUN=paste ,sep=""))
  BList <- vector(mode="numeric", length=0)
  RM<-c(CB+0:4*5)
  for (i in 1:ceiling(nrow(Test)/4)) {
    SBub<-dfXItem[,paste("Q",XItem[order(-hitN)[i]],"_",Bubs,sep="")]
    UL<-unlist(SBub)[-RM]    
    BList<-c(BList,UL)
  }
  plot(BList~c(1:length(BList)))
  Points=10
  BList<-BList[order(-BList)[1:Points]]
  plot(BList~c(1:Points))
  text(c(1:Points),BList, names(BList), cex=0.6, pos=4, col=2) 
  return(list(XBub=dfXItem,Hits=hitN,Bubs=BList))
}

PrePostItem <- function(it,classcut=NA,yearcut=NA,schoolcut=NA,diff=1) {

  print("The function calls for PrePostItem may include:")
  print("   0. first argument is the 1-129 integer index of the 129 test questions.")
  print("   1. classcut=integer, e.g. 115, 123, 159...; Course number")
  print("   2. yearcut=integer, e.g. 15, 16, 17, ...;   Year of graduation")
  print("   3. schoolcut=integer, e.g. 1, 11, 14, ...;  School code, i.e. Engineering = 14")
  print("   4. diff=1, choices: 1 for Pre Post, 2 for Post Postpost, 3 for Pre Postpost")
  print("         This allows comparison between any two pairs of cognitive diagnostic tests.")
  if (diff==1) {
    Ci = paste("C",it,sep="")
    Bi = paste("B",it,sep="")
    ci = paste("c",it,sep="")
    bi = paste("b",it,sep="")    
    tscore<-c("Score.x","Score.y")
    tstring<-c("PRE-TEST","POST-TEST")
  } else if (diff==2) {
    Ci = paste("c",it,sep="")
    Bi = paste("b",it,sep="")
    ci = paste("cc",it,sep="")
    bi = paste("bb",it,sep="")    
    tscore<-c("Score.y","Score.z")
    tstring<-c("POST-TEST","POST-POST-TEST")
  } else if (diff==3) {
    Ci = paste("C",it,sep="")
    Bi = paste("B",it,sep="")
    ci = paste("cc",it,sep="")
    bi = paste("bb",it,sep="")    
    tscore<-c("Score.x","Score.z")
    tstring<-c("PRE-TEST","POST-POST-TEST")
  }  
  print("Item indexed in Master List")
  descript = paste("Item",it)
  print(paste(descript,TestX[it,"Item"],sep=" : "))
  # Name of columns we want to analyze
  DF<-IDF
  CutIndex <- rep(0,4)
  CutString = c("Year = ","School = ","Class = ")
      
  if (!is.na(yearcut)) {
    DF<-DF[which(as.numeric(DF[,"Year"])==as.numeric(yearcut)),]    
    print(paste(CutString[1],yearcut,sep=""))    
  }
  if (!is.na(schoolcut)) {
    DF<-DF[which(as.numeric(DF[,"School"])==as.numeric(schoolcut)),]    
    print(paste(CutString[2],schoolcut,sep=""))
  }
  if (!is.na(classcut)) {
    CCut<-paste("Grade",as.character(classcut),sep="")
    DF<-DF[!is.na(DF[,CCut]),]    
    print(paste(CutString[3],classcut,sep=""))
  }
  if (nrow(DF)<5) print("Cuts are too strict or poorly defined, sample size is less than five")
  
  COLB <- DF[which(!is.na(DF[,Ci])),c(Ci,Bi,"Year","School",tscore[1],"MHI","GENDER_CD","ETHNIC_CD","SAT_V","SAT_M")]
  COLB[which(is.na(COLB[,2])),2]<-0
  for(i in 3:4)  COLB[,i]<-factor(COLB[,i])
  COLB[,2]<-factor(COLB[,2],levels=c("0","1","2","3","4","5"))
  COLN <- DF[which(!is.na(DF[,Ci])&!is.na(DF[,Bi])),c(Ci,Bi,"Year","School",tscore[1],"MHI","GENDER_CD","ETHNIC_CD","SAT_V","SAT_M")]
  
  COLb <- DF[which(!is.na(DF[,ci])),c(ci,bi,"Year","School",tscore[2],"MHI","GENDER_CD","ETHNIC_CD","SAT_V","SAT_M")]
  COLb[which(is.na(COLb[,2])),2]<-0
  for(i in 3:4)  COLb[,i]<-factor(COLb[,i])
  COLb[,2]<-factor(COLb[,2],levels=c("0","1","2","3","4","5"))
  COLn <- DF[which(!is.na(DF[,ci])&!is.na(DF[,bi])),c(ci,bi,"Year","School",tscore[2],"MHI","GENDER_CD","ETHNIC_CD","SAT_V","SAT_M")]
  
  remove(DF)
  pCB = sum(COLB[,1])/length(COLB[,1])
  pCb = sum(COLb[,1])/length(COLb[,1])
  pCN = sum(COLN[,1])/length(COLN[,1])
  pCn = sum(COLn[,1])/length(COLn[,1])
  sdB = sqrt(pCB*(1-pCB)/length(COLB[,1]))
  sdb = sqrt(pCb*(1-pCb)/length(COLb[,1]))
  sdN = sqrt(pCN*(1-pCN)/length(COLN[,1]))
  sdn = sqrt(pCn*(1-pCn)/length(COLn[,1]))
  nB = length(COLB[,1])-length(COLN[,1])
  nb = length(COLb[,1])-length(COLn[,1])
  CB = as.numeric(unique(COLB[which(COLB[,1]==1),2]))-1
  Cb = as.numeric(unique(COLb[which(COLb[,1]==1),2]))-1
  pBA = length(which(COLB[,2]==1))/length(COLB[,2])  
  pBB = length(which(COLB[,2]==2))/length(COLB[,2])  
  pBC = length(which(COLB[,2]==3))/length(COLB[,2])  
  pBD = length(which(COLB[,2]==4))/length(COLB[,2])  
  pBE = length(which(COLB[,2]==5))/length(COLB[,2])  
  pBa = length(which(COLb[,2]==1))/length(COLb[,2])  
  pBb = length(which(COLb[,2]==2))/length(COLb[,2])  
  pBc = length(which(COLb[,2]==3))/length(COLb[,2])  
  pBd = length(which(COLb[,2]==4))/length(COLb[,2])  
  pBe = length(which(COLb[,2]==5))/length(COLb[,2])  
  
  print("*****************************************************************")
  print(paste(tstring[1],"Results for Item:",it))
  print(paste("        Sample Size (Blanks included) : ",length(COLB[,1])))
  print(paste("        Students left this item blank : ",nB))
  print(paste("    Percent Correct (Blanks included) : ",format(100*pCB,digits=2,nsmall=2)))
  print(paste("       Binomial -  Standard Deviation :  ",format(100*sdB,digits=2,nsmall=2)))
  print(paste("Percent Correct (Blanks Not included) : ",format(100*pCN,digits=2,nsmall=2)))
  print(paste("       Binomial -  Standard Deviation :  ",format(100*sdN,digits=2,nsmall=2)))
  print("       ***********************************")
  print(paste("      Percent of students who chose A : ",format(100*pBA,digits=2,nsmall=2)))
  print(paste("      Percent of students who chose B : ",format(100*pBB,digits=2,nsmall=2)))
  print(paste("      Percent of students who chose C : ",format(100*pBC,digits=2,nsmall=2)))
  print(paste("      Percent of students who chose D : ",format(100*pBD,digits=2,nsmall=2)))
  print(paste("      Percent of students who chose E : ",format(100*pBE,digits=2,nsmall=2)))
  
  print("*****************************************************************")
  print(paste(tstring[2],"Results for Item:",it))
  print(paste("        Sample Size (Blanks included) : ",length(COLb[,1])))
  print(paste("        Students left this item blank : ",nb))
  print(paste("    Percent Correct (Blanks included) : ",format(100*pCb,digits=2,nsmall=2)))
  print(paste("       Binomial -  Standard Deviation :  ",format(100*sdb,digits=2,nsmall=2)))
  print(paste("Percent Correct (Blanks Not included) : ",format(100*pCn,digits=2,nsmall=2)))
  print(paste("       Binomial -  Standard Deviation :  ",format(100*sdn,digits=2,nsmall=2)))
  print("       ***********************************")
  print(paste("      Percent of students who chose A : ",format(100*pBa,digits=2,nsmall=2)))
  print(paste("      Percent of students who chose B : ",format(100*pBb,digits=2,nsmall=2)))
  print(paste("      Percent of students who chose C : ",format(100*pBc,digits=2,nsmall=2)))
  print(paste("      Percent of students who chose D : ",format(100*pBd,digits=2,nsmall=2)))
  print(paste("      Percent of students who chose E : ",format(100*pBe,digits=2,nsmall=2)))
  print("*****************************************************************")
  par(mfrow=c(2,2))
  ##############  GRAPH 1  ######################  
  TMP<-length(COLB[,2])
  title = paste(tstring[1],"\n",TestX[it,"ItemIndex"],"-",TestX[it,"Item"],": Bubbled Responses\n",sep=" ")
  Bubs = c("A","B","C","D","E")
  sampsize<-paste("N =",TMP)
  correct = paste("Correct Response is",paste(Bubs[CB],sampsize,sep="; "))
  title = paste(title,correct)
  plot(COLB[,2],xlab = "Student Answer", ylab = "Responses (%)", main = title, xaxt="n",yaxt="n",ylim=c(0,TMP))
  # draw an axis on the left
  axis(1, at=0.7+1.2*0:5, labels=c("Blank","A","B","C","D","E"))
  axis(2, at=c(0,TMP/2,TMP),labels=c(0,0.5,1.0))
  ##############  GRAPH 2  ######################  
  tmp<-length(COLb[,2])
  title = paste(tstring[2],"\n",TestY[it,"ItemIndex"],"-",TestY[it,"Item"],": Bubbled Responses\n",sep=" ")
  Bubs = c("A","B","C","D","E")
  sampsize<-paste("N =",tmp)
  correct = paste("Correct Response is",paste(Bubs[Cb],sampsize,sep="; "))
  title = paste(title,correct)
  plot(COLb[,2],xlab = "Student Answer", ylab = "Responses (%)", main = title, xaxt="n",yaxt="n",ylim=c(0,tmp))
  # draw an axis on the left
  axis(1, at=0.7+1.2*0:5, labels=c("Blank","A","B","C","D","E"))
  axis(2, at=c(0,tmp/2,tmp),labels=c(0,0.5,1.0))
  
  #   PropCor<-paste(format(100*pCB,digits=2,nsmall=2),"% Correct",sep="")
  #   title = paste(paste(paste(Test[it,"ItemIndex"],Test[it,"Item"],sep="-"),"Proportion correct\n",sep=" : "),PropCor)
  #   barplot(c(pCN,pCB), ylab = "Percent Correct", main = title, xaxt="n",ylim=c(0,1.0))
  #   axis(1, at=0.7+1.2*0:1, labels=c("No Blanks","Blank = Wrong"))
  
  #  COLB$ZS<-qqnorm(COLB$Score,main="Q-Q test for Normality of Scores")$x
  #  qqline(COLB$Score)
  
  breaksMath=unique(quantile(COLB$SAT_M,c(0,0.2,0.4,0.6,0.8,1.0),na.rm=TRUE))
  breaksMHI=unique(quantile(COLB$MHI,c(0,0.2,0.4,0.6,0.8,1.0),na.rm=TRUE))
  COLB$MathBreaks<-cut(COLB$SAT_M,breaks=breaksMath)
  COLB$MHIBreaks<-cut(COLB$MHI,breaks=breaksMHI)
  sumsMath<-by(COLB[,1], COLB$MathBreaks, sum)
  lengthsMath<-by(COLB[,1], COLB$MathBreaks, length)
  sumsMHI<-by(COLB[,1], COLB$MHIBreaks, sum)
  lengthsMHI<-by(COLB[,1], COLB$MHIBreaks, length)
  
  breaksmath=unique(quantile(COLb$SAT_M,c(0,0.2,0.4,0.6,0.8,1.0),na.rm=TRUE))
  breaksmhi=unique(quantile(COLb$MHI,c(0,0.2,0.4,0.6,0.8,1.0),na.rm=TRUE))
  COLb$MathBreaks<-cut(COLb$SAT_M,breaks=breaksmath)
  COLb$MHIBreaks<-cut(COLb$MHI,breaks=breaksmhi)
  sumsmath<-by(COLb[,1], COLb$MathBreaks, sum)
  lengthsmath<-by(COLb[,1], COLb$MathBreaks, length)
  sumsmhi<-by(COLb[,1], COLb$MHIBreaks, sum)
  lengthsmhi<-by(COLb[,1], COLb$MHIBreaks, length)
  
  print(paste(tstring[1],"SAT Quintiles"))
  print(breaksMath)
  print(paste(tstring[2],"SAT Quintiles"))
  print(breaksmath)
  print(paste(tstring[1],"MHI Quintiles"))
  print(breaksMHI)
  print(paste(tstring[2],"MHI Quintiles"))
  print(breaksmhi)
  
  #   sums<-by(COLB[,1], COLB$Breaks, sum)
  #   lengths<-by(COLB[,1], COLB$Breaks, length)
  #   plot(sums/lengths~c(-0.05+0.1*1:10),xlab="Quantiles (Z)")
  #   fit <- glm(cbind(sums, lengths-sums) ~ c(-0.05+0.1*1:10), family=binomial)
  #   lines(c(-0.05+0.1*1:10), fit$fitted, type="l", col="red")
  
  plot(sumsMath/lengthsMath~c(0.1+0:4*0.2), ylim=c(0,1.0), main = paste("Item Success by SAT Quintile\n Red = ",tstring[1],", Blue = ",tstring[2],sep=""), xlab="SAT Quintiles", ylab="Percent Correct",pch=16,col=2)
  points(sumsmath/lengthsmath~c(0.1+0:4*0.2),pch=16,col=4)
  plot(sumsMHI/lengthsMHI~c(0.1+0:4*0.2), ylim=c(0,1.0), main = paste("Item Success by MHI Quintile\n Red = ",tstring[1],", Blue = ",tstring[2],sep=""), xlab="SAT Quintiles", ylab="Percent Correct",pch=16,col=2)
  points(sumsmhi/lengthsmhi~c(0.1+0:4*0.2),pch=16,col=4)
}
