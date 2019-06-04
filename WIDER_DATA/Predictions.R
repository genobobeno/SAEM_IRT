#################################################################
#  This file holds the functions for predicting Passing by grade, or some other measure
#     1. PredictFailing()  <- Needs cleaning
#     2. PredictPassing()  <- Cleaned 

PredictFailing <- function(yearcut=17,classcut=100,gradesplit="C",IncludeW=TRUE) {
  #   gradesplit="C"
  #   IncludeW=TRUE
  GradeList<-c("A","B+","B","C+","C","D","F","W")
  #label<-"Target"
  if (IncludeW==FALSE) { 
    GradeList<-GradeList[-length(GradeList)] 
  }
  par(mfrow=c(2,2))
  plot(as.factor(IDF$LetterGrade115),main="115 Grade Distribution")
  plot(as.factor(IDF$LetterGrade123),main="123 Grade Distribution")
  plot(as.factor(IDF$LetterGrade159),main="159 Grade Distribution")
  plot(as.factor(IDF$LetterGrade227),main="227 Grade Distribution")
  par(mfrow=c(1,2))
  plot(IDF$Grade159,IDF$Grade123)
  plot(IDF$Grade159,IDF$Grade115)
  #   print("Rows are 123, Columns are 159") 
  #   table(IDF$LetterGrade123,IDF$LetterGrade159)
  GCut<-match(gradesplit,GradeList)
  Win<-GradeList[1:GCut]
  print(paste("Passing Grades are:",Win))
  Lose<-GradeList[(GCut+1):length(GradeList)]
  label<-""
  for (i in (GCut+1):length(GradeList)) {label<-paste(label,GradeList[i],sep="")}
  Predictors<-c("MHI","SAT_M","SAT_V","GENDER_CD","ETHNIC_CD","Score.x",paste("C",1:nrow(TestX),sep=""))
  Predict<-c("MHI","SAT_M","SAT_V","Score.x","FCI.x","MATH_REASON.x","CCI.x","MBL.x","ACS.x",paste("C",1:nrow(TestX),sep=""))  
  Q<-paste("C",1:nrow(TestX),sep="")
  QDisc<-data.frame(Q=1:nrow(TestX),Fr=rep(0,nrow(TestX)),Meet=rep(0,nrow(TestX)),G115=rep(0,nrow(TestX)),G123=rep(0,nrow(TestX)),G159=rep(0,nrow(TestX)),G227=rep(0,nrow(TestX)))
  
  if (classcut==100) {
    Meet<-IDF[which(!is.na(IDF$Grade123)&!is.na(IDF$Grade159)),c(Predictors,"Grade159","Grade123","FCI.x","MATH_REASON.x","CCI.x","MBL.x","ACS.x","LetterGrade115") %in% colnames(IDF)]
    Meet$Pass<-(Meet$Grade123>=67|Meet$Grade159>=315) + 0
    for (i in 1:nrow(TestX)) {
      test<-Meet[which(!is.na(Meet[,Q[i]])),c("Pass",Q[i])]
      if (nrow(test)/length(Meet$Pass)>0.3) { QDisc[i,"Meet"]<-mean(test[which(test$Pass==1),Q[i]],na.rm=T)-mean(test[which(test$Pass==0),Q[i]],na.rm=T) }
    }
    QS<-QDisc[order(-QDisc$Meet),c("Q","Meet")]    
    FROSH<-Meet
    label<-"Target"
    title<-"Target"
    print(Meet[,paste("C",QS[1:20,1],sep="")])
    
  }
  else if (classcut==200) {
    Meet<-IDF[which(!is.na(IDF$Grade123)&!is.na(IDF$Grade159)),c(Predictors,"Grade159","Grade123","FCI.x","MATH_REASON.x","CCI.x","MBL.x","ACS.x","LetterGrade115") %in% colnames(IDF)]
    Meet$Pass<-(Meet$Grade123<=73|Meet$Grade159<=336) + 0
    for (i in 1:nrow(TestX)) {
      test<-Meet[which(!is.na(Meet[,Q[i]])),c("Pass",Q[i])]
      if (nrow(test)/length(Meet$Pass)>0.3) { QDisc[i,"Meet"]<-mean(test[which(test$Pass==1),Q[i]],na.rm=T)-mean(test[which(test$Pass==0),Q[i]],na.rm=T) }
    }
    QS<-QDisc[order(-QDisc$Meet),c("Q","Meet")]    
    FROSH<-Meet
    label<-"AB"
    title<-"AB"
    print(Meet[,paste("C",QS[1:20,1],sep="")])
    
  }
  else if (classcut==115) {
    G115<-IDF[which(!is.na(IDF$LetterGrade115)),c(Predictors,"FCI.x","MATH_REASON.x","CCI.x","MBL.x","ACS.x","LetterGrade115") %in% colnames(IDF)]
    G115$Pass<-(G115$LetterGrade115 %in% Win) + 0
    for (i in 1:nrow(TestX)) {
      test115<-G115[which(!is.na(G115[,Q[i]])),c("Pass",Q[i])]
      if (nrow(test115)/length(G115$Pass)>0.25) { QDisc[i,"G115"]<-mean(test115[which(test115$Pass==1),Q[i]],na.rm=T)-mean(test115[which(test115$Pass==0),Q[i]],na.rm=T) }
    }
    QS<-QDisc[order(-QDisc$G115),c("Q","G115")]    
    FROSH<-G115
    title<-"115"
    label<-paste(label,title)
    print(G115[,paste("C",QS[1:20,1],sep="")])
  }
  else if (classcut==123) {
    G123<-IDF[which(!is.na(IDF$LetterGrade123)),c(Predictors,"FCI.x","MATH_REASON.x","CCI.x","MBL.x","ACS.x","LetterGrade123") %in% colnames(IDF)]
    G123$Pass<-(G123$LetterGrade123 %in% Win) + 0
    for (i in 1:nrow(TestX)) {
      test123<-G123[which(!is.na(G123[,Q[i]])),c("Pass",Q[i])]
      if (nrow(test123)/length(G123$Pass)>0.25) { QDisc[i,"G123"]<-mean(test123[which(test123$Pass==1),Q[i]],na.rm=T)-mean(test123[which(test123$Pass==0),Q[i]],na.rm=T) }
    }
    QS<-QDisc[order(-QDisc$G123),c("Q","G123")]    
    FROSH<-G123
    title<-"123"
    label<-paste(label,title)
  }
  else if (classcut==159) {
    G159<-IDF[which(!is.na(IDF$LetterGrade159)),c(Predictors,"FCI.x","MATH_REASON.x","CCI.x","MBL.x","ACS.x","LetterGrade159") %in% colnames(IDF)]
    G159$Pass<-(G159$LetterGrade159 %in% Win) + 0
    for (i in 1:nrow(TestX)) {
      test159<-G159[which(!is.na(G159[,Q[i]])),c("Pass",Q[i])]
      if (nrow(test159)/length(G159$Pass)>0.25) { QDisc[i,"G159"]<-mean(test159[which(test159$Pass==1),Q[i]],na.rm=T)-mean(test159[which(test159$Pass==0),Q[i]],na.rm=T) }
    }
    QS<-QDisc[order(-QDisc$G159),c("Q","G159")]
    print(G159[,paste("C",QS[1:20,1],sep="")])
    FROSH<-G159
    title<-"159"
    label<-paste(label,title)
  }
  else if (classcut==227) {
    G227<-IDF[which(!is.na(IDF$LetterGrade227)),c(Predictors,"CSEM.x","MATH_REASON.x","LetterGrade227") %in% colnames(IDF)]
    G227$Pass<-(G227$LetterGrade227 %in% Win) + 0
    for (i in 1:nrow(TestX)) {
      test227<-G227[which(!is.na(G227[,Q[i]])),c("Pass",Q[i])]
      if (nrow(test227)/length(G227$Pass)>0.25) { QDisc[i,"G227"]<-mean(test227[which(test227$Pass==1),Q[i]],na.rm=T)-mean(test227[which(test227$Pass==0),Q[i]],na.rm=T) }
    }
    QS<-QDisc[order(-QDisc$G227),c("Q","G227")]
    FROSH<-G227
    title<-"227"
    label<-paste(label,title)
  }
  else {  
    FRESH<-IDF[which(!is.na(IDF$LetterGrade123)&!is.na(IDF$LetterGrade159)),c(Predictors,"FCI.x","MATH_REASON.x","CCI.x","MBL.x","ACS.x","Grade123","Grade159","LetterGrade123","LetterGrade159") %in% colnames(IDF)]
    FRESH$Pass<-((FRESH$LetterGrade123 %in% Win) + 0)*((FRESH$LetterGrade159 %in% Win) + 0)
    for (i in 1:nrow(TestX)) {
      testF<-FRESH[which(!is.na(FRESH[,Q[i]])),c("Pass",Q[i])]
      if (nrow(testF)>40) { QDisc[i,"Fr"]<-mean(testF[which(testF$Pass==1),Q[i]],na.rm=T)-mean(testF[which(testF$Pass==0),Q[i]],na.rm=T) }
    }
    QS<-QDisc[order(-QDisc$Fr),c("Q","Fr")]
    FROSH<-FRESH
    title<-"123 & 159"  
    label<-paste(label,title)
  }
  
  print(QS[1:20,])
  BestQ<-paste("C",QS$Q[1:12],sep="")
  GoodQ<-paste("C",QS$Q[1:20],sep="")
  FROSH$T1<-rowMeans(FROSH[,BestQ],na.rm=TRUE)
  FROSH$T2<-rowMeans(FROSH[,GoodQ],na.rm=TRUE)
  
  #BestQ<-paste("C",c(2,6,8,68,71,73,74,83,89,99,116,124),sep="")
  #GoodQ<-paste("C",c(2,6,8,68,71,73,74,83,89,99,116,124,4,76,82,91,110,113,117,123),sep="")
  # Most significant - 73,71
  
  #logitF <- glm(Pass ~ C2+C6+C8+C68+C71+C73+C74+C83+C89+C99+C116+C124, data = test, family = "binomial")
  par(mfrow=c(2,4))
  if (classcut!=227) {
    plot(density(FROSH$CCI.x[which(FROSH$Pass==0)],na.rm=TRUE),main=paste("CCI Scores \nBlue - Passed",title,"\n Red -",label),col=2,xlab="CCI Pretest Score",xlim=c(0,1.0))
    lines(density(FROSH$CCI.x[which(FROSH$Pass==1)],na.rm=TRUE),col=4)
  }
  else {
    plot(density(FROSH$CSEM.x[which(FROSH$Pass==0)],na.rm=TRUE),main=paste("CSEM Scores \nBlue - Passed",title,"\n Red -",label),col=2,xlab="CSEM Pretest Score",xlim=c(0,1.0))
    lines(density(FROSH$CSEM.x[which(FROSH$Pass==1)],na.rm=TRUE),col=4)
  }    
  plot(density(FROSH$FCI.x[which(FROSH$Pass==0)],na.rm=TRUE),main=paste("FCI Scores \nBlue - Passed",title,"\n Red -",label),col=2,xlab="FCI Pretest Score",xlim=c(0,1.0))
  lines(density(FROSH$FCI.x[which(FROSH$Pass==1)],na.rm=TRUE),col=4)
  plot(density(FROSH$T1[which(FROSH$Pass==0)],na.rm=TRUE),main=paste("12 Question Test \nBlue - Passed",title,"\n Red -",label),col=2,xlab="12Q Score",xlim=c(0,1.0))
  lines(density(FROSH$T1[which(FROSH$Pass==1)],na.rm=TRUE),col=4)
  plot(density(FROSH$T2[which(FROSH$Pass==0)],na.rm=TRUE),main=paste("20 Question Test \nBlue - Passed",title,"\n Red -",label),col=2,xlab="20Q Score",xlim=c(0,1.0))
  lines(density(FROSH$T2[which(FROSH$Pass==1)],na.rm=TRUE),col=4)
  if (classcut!=227) {
    hist(FROSH$CCI.x[which(FROSH$Pass==1)],breaks=c(0:10*0.1),main=paste("CCI Scores \nBlue - Passed",title,"\n Red -",label),col=4,xlab="CCI Score",xlim=c(0,1.0))
    hist(FROSH$CCI.x[which(FROSH$Pass==0)],breaks=c(0:10*0.1),add=T,col=2)
  }
  else {
    hist(FROSH$CSEM.x[which(FROSH$Pass==1)],breaks=c(0:10*0.1),main=paste("CSEM Scores \nBlue - Passed",title,"\n Red -",label),col=4,xlab="CSEM Score",xlim=c(0,1.0))
    hist(FROSH$CSEM.x[which(FROSH$Pass==0)],add=T,col=2,breaks=c(0:10*0.1))
  }
  hist(FROSH$FCI.x[which(FROSH$Pass==1)],breaks=c(0:10*0.1),main=paste("FCI Scores \nBlue - Passed",title,"\n Red -",label),col=4,xlab="FCI Score",xlim=c(0,1.0))
  hist(FROSH$FCI.x[which(FROSH$Pass==0)],add=T,breaks=c(0:10*0.1),col=2)
  hist(FROSH$T1[which(FROSH$Pass==1)],breaks=c(0:10*0.1),main=paste("12 Question Test \nBlue - Passed",title,"\n Red -",label),col=4,xlab="12Q Score",xlim=c(0,1.0))
  hist(FROSH$T1[which(FROSH$Pass==0)],add=T,breaks=c(0:10*0.1),col=2)
  hist(FROSH$T2[which(FROSH$Pass==1)],breaks=c(0:10*0.1),main=paste("20 Question Test \nBlue - Passed",title,"\n Red -",label),col=4,xlab="20Q Score",xlim=c(0,1.0))
  hist(FROSH$T2[which(FROSH$Pass==0)],add=T,breaks=c(0:10*0.1),col=2)
  
  
  par(mfrow=c(2,2))
  FROSH<-FROSH[which(!is.na(FROSH$T1)),]
  N<-nrow(FROSH)
  Passed<-length(which(FROSH$Pass==1))
  WDF<-length(which(FROSH$Pass==0))
  print(paste("Total % of students who did NOT pass:",WDF/N))
  CID<-vector(length=0,mode="numeric")
  Type1<-vector(length=0,mode="numeric")
  Type2<-vector(length=0,mode="numeric")
  Stud<-vector(length=0,mode="numeric")
  Fails<-vector(length=0,mode="numeric")
  for (i in 0:length(BestQ)) {
    Test1<-subset(FROSH,T1>=1/length(BestQ)*i)
    Stud<-c(Stud,nrow(Test1)/N) # % passing cut
    CID<-c(CID,(length(which(FROSH$Pass==0))-length(which(Test1$Pass==0)))/(N-nrow(Test1)))
    Type1<-c(Type1,(length(which(FROSH$Pass==1))-length(which(Test1$Pass==1)))/(N-nrow(Test1)))
    Type2<-c(Type2,length(which(Test1$Pass==0))/nrow(Test1))
    Fails<-c(Fails,1-length(which(Test1$Pass==0))/WDF)
  }
  plot(c(0:(length(Type1)-1))/(length(Type1)-1),Type1,type="n",main=paste("%",label,"as a function of \n12-question cutoff Score"),xlab="Cut Grade",ylim=c(0,1.0),ylab=paste("%",label))
  lines(c(0:(length(Type1)-1))/(length(Type1)-1),Stud,col=1,lty=1)
  lines(c(0:(length(Type1)-1))/(length(Type1)-1),CID,col=3,lty=1)
  lines(c(0:(length(Type1)-1))/(length(Type1)-1),Type1,col=2,lty=2)
  lines(c(0:(length(Type1)-1))/(length(Type1)-1),Type2,col=6,lty=2)
  lines(c(0:(length(Type1)-1))/(length(Type1)-1),Fails,col=4,lty=5)
  legend(0.65,0.75,c("Pass","CorrID","Type1","Type2","Fail"),lty=c(1,1,2,2,5),col=c(1,3,2,6,4))
  
  CID<-vector(length=0,mode="numeric")
  Type1<-vector(length=0,mode="numeric")
  Type2<-vector(length=0,mode="numeric")
  Stud<-vector(length=0,mode="numeric")
  Fails<-vector(length=0,mode="numeric")
  for (i in 0:length(GoodQ)) {
    Test1<-subset(FROSH,T2>=1/length(GoodQ)*i)
    Stud<-c(Stud,nrow(Test1)/N) # % passing cut
    CID<-c(CID,(length(which(FROSH$Pass==0))-length(which(Test1$Pass==0)))/(N-nrow(Test1)))
    Type1<-c(Type1,(length(which(FROSH$Pass==1))-length(which(Test1$Pass==1)))/(N-nrow(Test1)))
    Type2<-c(Type2,length(which(Test1$Pass==0))/nrow(Test1))
    Fails<-c(Fails,1-length(which(Test1$Pass==0))/WDF)
  }
  plot(c(0:(length(Type1)-1))/(length(Type1)-1),Type1,type="n",main=paste("%",label,"as a function of \n20-question cutoff Score"),xlab="Cut Grade",ylim=c(0,1.0),ylab=paste("%",label))
  lines(c(0:(length(Type1)-1))/(length(Type1)-1),Stud,col=1,lty=1)
  lines(c(0:(length(Type1)-1))/(length(Type1)-1),CID,col=3,lty=1)
  lines(c(0:(length(Type1)-1))/(length(Type1)-1),Type1,col=2,lty=2)
  lines(c(0:(length(Type1)-1))/(length(Type1)-1),Type2,col=6,lty=2)
  lines(c(0:(length(Type1)-1))/(length(Type1)-1),Fails,col=4,lty=5)
  legend(0.65,0.75,c("Pass","CorrID","Type1","Type2","Fail"),lty=c(1,1,2,2,5),col=c(1,3,2,6,4))
  
  CID<-vector(length=0,mode="numeric")
  Type1<-vector(length=0,mode="numeric")
  Type2<-vector(length=0,mode="numeric")
  Stud<-vector(length=0,mode="numeric")
  Fails<-vector(length=0,mode="numeric")
  for (i in 0:nrow(FCTest)) {
    Test1<-subset(FROSH,FCI.x>=1/nrow(FCTest)*i)
    Stud<-c(Stud,nrow(Test1)/N) # % passing cut
    CID<-c(CID,(length(which(FROSH$Pass==0))-length(which(Test1$Pass==0)))/(N-nrow(Test1)))
    Type1<-c(Type1,(length(which(FROSH$Pass==1))-length(which(Test1$Pass==1)))/(N-nrow(Test1)))
    Type2<-c(Type2,length(which(Test1$Pass==0))/nrow(Test1))
    Fails<-c(Fails,1-length(which(Test1$Pass==0))/WDF)
  }
  plot(c(0:(length(Type1)-1))/(length(Type1)-1),Type1,type="n",main=paste("%",label,"as a function of \nFCI cutoff Score"),xlab="Cut Grade",ylim=c(0,1.0),ylab=paste("%",label))
  lines(c(0:(length(Type1)-1))/(length(Type1)-1),Stud,col=1,lty=1)
  lines(c(0:(length(Type1)-1))/(length(Type1)-1),CID,col=3,lty=1)
  lines(c(0:(length(Type1)-1))/(length(Type1)-1),Type1,col=2,lty=2)
  lines(c(0:(length(Type1)-1))/(length(Type1)-1),Type2,col=6,lty=2)
  lines(c(0:(length(Type1)-1))/(length(Type1)-1),Fails,col=4,lty=5)
  legend(0.65,0.75,c("Pass","CorrID","Type1","Type2","Fail"),lty=c(1,1,2,2,5),col=c(1,3,2,6,4))
  
  CID<-vector(length=0,mode="numeric")
  Type1<-vector(length=0,mode="numeric")
  Type2<-vector(length=0,mode="numeric")
  Stud<-vector(length=0,mode="numeric")
  Fails<-vector(length=0,mode="numeric")
  if (classcut!=227) {
    for (i in 0:length(which(TestX$Category=="CCI"))) {
      Test1<-subset(FROSH,CCI.x>=1/length(which(TestX$Category=="CCI"))*i)
      Stud<-c(Stud,nrow(Test1)/N) # % passing cut
      CID<-c(CID,(length(which(FROSH$Pass==0))-length(which(Test1$Pass==0)))/(N-nrow(Test1)))
      Type1<-c(Type1,(length(which(FROSH$Pass==1))-length(which(Test1$Pass==1)))/(N-nrow(Test1)))
      Type2<-c(Type2,length(which(Test1$Pass==0))/nrow(Test1))
      Fails<-c(Fails,1-length(which(Test1$Pass==0))/WDF)
    }
    plot(c(0:(length(Type1)-1))/(length(Type1)-1),Type1,type="n",main=paste("%",label,"as a function of \nCCI cutoff Score"),xlab="Cut Grade",ylim=c(0,1.0),ylab=paste("%",label))
  }
  else {
    for (i in 0:length(which(TestX$Category=="CSEM"))) {
      Test1<-subset(FROSH,CSEM.x>=1/length(which(TestX$Category=="CSEM"))*i)
      Stud<-c(Stud,nrow(Test1)/N) # % passing cut
      CID<-c(CID,(length(which(FROSH$Pass==0))-length(which(Test1$Pass==0)))/(N-nrow(Test1)))
      Type1<-c(Type1,(length(which(FROSH$Pass==1))-length(which(Test1$Pass==1)))/(N-nrow(Test1)))
      Type2<-c(Type2,length(which(Test1$Pass==0))/nrow(Test1))
      Fails<-c(Fails,1-length(which(Test1$Pass==0))/WDF)
    }
    plot(c(0:(length(Type1)-1))/(length(Type1)-1),Type1,type="n",main=paste("%",label,"as a function of \nCSEM cutoff Score"),xlab="Cut Grade",ylim=c(0,1.0),ylab=paste("%",label))
  }
  lines(c(0:(length(Type1)-1))/(length(Type1)-1),Stud,col=1,lty=1)
  lines(c(0:(length(Type1)-1))/(length(Type1)-1),CID,col=3,lty=1)
  lines(c(0:(length(Type1)-1))/(length(Type1)-1),Type1,col=2,lty=2)
  lines(c(0:(length(Type1)-1))/(length(Type1)-1),Type2,col=6,lty=2)
  lines(c(0:(length(Type1)-1))/(length(Type1)-1),Fails,col=4,lty=5)
  legend(0.65,0.75,c("Pass","CorrID","Type1","Type2","Fail"),lty=c(1,1,2,2,5),col=c(1,3,2,6,4))
  
  ###########  12 Question Test Analysis  #################### 
  # SAT_MATH
  par(mfrow=c(2,2))
  CID<-vector(length=0,mode="numeric")
  Type1<-vector(length=0,mode="numeric")
  Type2<-vector(length=0,mode="numeric")
  Stud<-vector(length=0,mode="numeric")
  Fails<-vector(length=0,mode="numeric")
  for (i in 0:length(BestQ)) {
    Test1<-subset(FROSH,T1>=1/length(BestQ)*i)
    Stud<-c(Stud,nrow(Test1)/N) # % passing cut
    CID<-c(CID,(length(which(FROSH$Pass==0))-length(which(Test1$Pass==0)))/(N-nrow(Test1)))
    Type1<-c(Type1,(length(which(FROSH$Pass==1))-length(which(Test1$Pass==1)))/(N-nrow(Test1)))
    Type2<-c(Type2,length(which(Test1$Pass==0))/nrow(Test1))
    Fails<-c(Fails,1-length(which(Test1$Pass==0))/WDF)
  }
  plot(c(0:(length(Type1)-1))/(length(Type1)-1),Type1,type="n",main=paste("%",label,"as a function of \n12-Question Score cutoff"),xlab="Cut Grade",ylim=c(0,1.0),ylab=paste(" %",label))
  lines(c(0:(length(Type1)-1))/(length(Type1)-1),Stud,col=1,lty=1)
  lines(c(0:(length(Type1)-1))/(length(Type1)-1),CID,col=3,lty=1)
  lines(c(0:(length(Type1)-1))/(length(Type1)-1),Type1,col=2,lty=2)
  lines(c(0:(length(Type1)-1))/(length(Type1)-1),Type2,col=6,lty=2)
  lines(c(0:(length(Type1)-1))/(length(Type1)-1),Fails,col=4,lty=5)
  legend(0.65,0.75,c("Pass","CorrID","Type1","Type2","Fail"),lty=c(1,1,2,2,5),col=c(1,3,2,6,4))
  
  plot(FROSH$SAT_M[which(FROSH$Pass==1)],FROSH$T1[which(FROSH$Pass==1)],col=4,pch=16,main="12 Question Test vs. SAT MATH",xlab="SAT Math Scores",ylab="20 Question Test Score")
  points(FROSH$SAT_M[which(FROSH$Pass==0)],FROSH$T1[which(FROSH$Pass==0)],col=2,pch=16)
  abline(b=-0.17/10,a=1.2)
  abline(b=-0.3/10,a=2.0)
  
  A<--0.03
  B<-2.0
  CID<-vector(length=0,mode="numeric")
  Type1<-vector(length=0,mode="numeric")
  Type2<-vector(length=0,mode="numeric")
  Stud<-vector(length=0,mode="numeric")
  Fails<-vector(length=0,mode="numeric")
  for (i in 0:60) {
    Test1<-subset(FROSH,T1>=SAT_M*A+B+0.01*i)
    Stud<-c(Stud,nrow(Test1)/N) # % passing cut
    CID<-c(CID,(length(which(FROSH$Pass==0))-length(which(Test1$Pass==0)))/(N-nrow(Test1)))
    Type1<-c(Type1,(length(which(FROSH$Pass==1))-length(which(Test1$Pass==1)))/(N-nrow(Test1)))
    Type2<-c(Type2,length(which(Test1$Pass==0))/nrow(Test1))
    Fails<-c(Fails,1-length(which(Test1$Pass==0))/WDF)
  }
  plot(c(0:(length(Type1)-1))*0.01+B,Type1,type="n",main=paste("%",label,"as a function of \n12-Question Score cutoff"),xlab=paste("Cut Intercept, Slope =",A),ylim=c(0,1.0),ylab=paste("%",label))
  lines(c(0:(length(Type1)-1))*0.01+B,Stud,col=1,lty=1)
  lines(c(0:(length(Type1)-1))*0.01+B,CID,col=3,lty=1)
  lines(c(0:(length(Type1)-1))*0.01+B,Type1,col=2,lty=2)
  lines(c(0:(length(Type1)-1))*0.01+B,Type2,col=6,lty=2)
  lines(c(0:(length(Type1)-1))*0.01+B,Fails,col=4,lty=5)
  legend(0.65,0.75,c("Pass","CorrID","Type1","Type2","Fail"),lty=c(1,1,2,2,5),col=c(1,3,2,6,4))
  
  A<--0.017
  B<-1.2
  CID<-vector(length=0,mode="numeric")
  Type1<-vector(length=0,mode="numeric")
  Type2<-vector(length=0,mode="numeric")
  Stud<-vector(length=0,mode="numeric")
  Fails<-vector(length=0,mode="numeric")
  for (i in 0:60) {
    Test1<-subset(FROSH,T1>=SAT_M*A+B+0.01*i)
    Stud<-c(Stud,nrow(Test1)/N) # % passing cut
    CID<-c(CID,(length(which(FROSH$Pass==0))-length(which(Test1$Pass==0)))/(N-nrow(Test1)))
    Type1<-c(Type1,(length(which(FROSH$Pass==1))-length(which(Test1$Pass==1)))/(N-nrow(Test1)))
    Type2<-c(Type2,length(which(Test1$Pass==0))/nrow(Test1))
    Fails<-c(Fails,1-length(which(Test1$Pass==0))/WDF)
  }
  plot(c(0:(length(Type1)-1))*0.01+B,Type1,type="n",main=paste("%",label,"as a function of \n12-Question Score cutoff"),xlab=paste("Cut Intercept, Slope =",A),ylim=c(0,1.0),ylab=paste("%",label))
  lines(c(0:(length(Type1)-1))*0.01+B,Stud,col=1,lty=1)
  lines(c(0:(length(Type1)-1))*0.01+B,CID,col=3,lty=1)
  lines(c(0:(length(Type1)-1))*0.01+B,Type1,col=2,lty=2)
  lines(c(0:(length(Type1)-1))*0.01+B,Type2,col=6,lty=2)
  lines(c(0:(length(Type1)-1))*0.01+B,Fails,col=4,lty=5)
  legend(0.65,0.75,c("Pass","CorrID","Type1","Type2","Fail"),lty=c(1,1,2,2,5),col=c(1,3,2,6,4))
  
  # SAT_TOTAL
  par(mfrow=c(2,2))
  CID<-vector(length=0,mode="numeric")
  Type1<-vector(length=0,mode="numeric")
  Type2<-vector(length=0,mode="numeric")
  Stud<-vector(length=0,mode="numeric")
  Fails<-vector(length=0,mode="numeric")
  for (i in 0:length(BestQ)) {
    Test1<-subset(FROSH,T1>=1/length(BestQ)*i)
    Stud<-c(Stud,nrow(Test1)/N) # % passing cut
    CID<-c(CID,(length(which(FROSH$Pass==0))-length(which(Test1$Pass==0)))/(N-nrow(Test1)))
    Type1<-c(Type1,(length(which(FROSH$Pass==1))-length(which(Test1$Pass==1)))/(N-nrow(Test1)))
    Type2<-c(Type2,length(which(Test1$Pass==0))/nrow(Test1))
    Fails<-c(Fails,1-length(which(Test1$Pass==0))/WDF)
  }
  plot(c(0:(length(Type1)-1))/(length(Type1)-1),Type1,type="n",main=paste("%",label,"as a function of \n12-Question Score cutoff"),xlab="Cut Grade",ylim=c(0,1.0),ylab=paste("%",label))
  lines(c(0:(length(Type1)-1))/(length(Type1)-1),Stud,col=1,lty=1)
  lines(c(0:(length(Type1)-1))/(length(Type1)-1),CID,col=3,lty=1)
  lines(c(0:(length(Type1)-1))/(length(Type1)-1),Type1,col=2,lty=2)
  lines(c(0:(length(Type1)-1))/(length(Type1)-1),Type2,col=6,lty=2)
  lines(c(0:(length(Type1)-1))/(length(Type1)-1),Fails,col=4,lty=5)
  legend(0.65,0.75,c("Pass","CorrID","Type1","Type2","Fail"),lty=c(1,1,2,2,5),col=c(1,3,2,6,4))
  
  plot((FROSH$SAT_M[which(FROSH$Pass==1)]+FROSH$SAT_V[which(FROSH$Pass==1)]),FROSH$T1[which(FROSH$Pass==1)],col=4,pch=16,main="12 Question Test vs. SAT TOTAL",xlab="SAT Total Scores",ylab="20 Question Test Score")
  points((FROSH$SAT_M[which(FROSH$Pass==0)]+FROSH$SAT_V[which(FROSH$Pass==0)]),FROSH$T1[which(FROSH$Pass==0)],col=2,pch=16)
  abline(b=-0.17/10,a=2.4)
  abline(b=-0.3/10,a=3.9)
  
  A<--0.03
  B<-3.7
  CID<-vector(length=0,mode="numeric")
  Type1<-vector(length=0,mode="numeric")
  Type2<-vector(length=0,mode="numeric")
  Stud<-vector(length=0,mode="numeric")
  Fails<-vector(length=0,mode="numeric")
  for (i in 0:60) {
    Test1<-subset(FROSH,T1>=(SAT_M+SAT_V)*A+B+0.01*i)
    Stud<-c(Stud,nrow(Test1)/N) # % passing cut
    CID<-c(CID,(length(which(FROSH$Pass==0))-length(which(Test1$Pass==0)))/(N-nrow(Test1)))
    Type1<-c(Type1,(length(which(FROSH$Pass==1))-length(which(Test1$Pass==1)))/(N-nrow(Test1)))
    Type2<-c(Type2,length(which(Test1$Pass==0))/nrow(Test1))
    Fails<-c(Fails,1-length(which(Test1$Pass==0))/WDF)
  }
  plot(c(0:(length(Type1)-1))*0.01+B,Type1,type="n",main=paste("%",label,"as a function of \n12-Question Score cutoff"),xlab=paste("Cut Intercept, Slope =",A),ylim=c(0,1.0),ylab=paste("%",label))
  lines(c(0:(length(Type1)-1))*0.01+B,Stud,col=1,lty=1)
  lines(c(0:(length(Type1)-1))*0.01+B,CID,col=3,lty=1)
  lines(c(0:(length(Type1)-1))*0.01+B,Type1,col=2,lty=2)
  lines(c(0:(length(Type1)-1))*0.01+B,Type2,col=6,lty=2)
  lines(c(0:(length(Type1)-1))*0.01+B,Fails,col=4,lty=5)
  legend(0.65,0.75,c("Pass","CorrID","Type1","Type2","Fail"),lty=c(1,1,2,2,5),col=c(1,3,2,6,4))
  
  A<--0.017
  B<-2.1
  CID<-vector(length=0,mode="numeric")
  Type1<-vector(length=0,mode="numeric")
  Type2<-vector(length=0,mode="numeric")
  Stud<-vector(length=0,mode="numeric")
  Fails<-vector(length=0,mode="numeric")
  for (i in 0:60) {
    Test1<-subset(FROSH,T1>=(SAT_M+SAT_V)*A+B+0.01*i)
    Stud<-c(Stud,nrow(Test1)/N) # % passing cut
    CID<-c(CID,(length(which(FROSH$Pass==0))-length(which(Test1$Pass==0)))/(N-nrow(Test1)))
    Type1<-c(Type1,(length(which(FROSH$Pass==1))-length(which(Test1$Pass==1)))/(N-nrow(Test1)))
    Type2<-c(Type2,length(which(Test1$Pass==0))/nrow(Test1))
    Fails<-c(Fails,1-length(which(Test1$Pass==0))/WDF)
  }
  plot(c(0:(length(Type1)-1))*0.01+B,Type1,type="n",main=paste("%",label,"as a function of \n12-Question Score cutoff"),xlab=paste("Cut Intercept, Slope =",A),ylim=c(0,1.0),ylab=paste("%",label))
  lines(c(0:(length(Type1)-1))*0.01+B,Stud,col=1,lty=1)
  lines(c(0:(length(Type1)-1))*0.01+B,CID,col=3,lty=1)
  lines(c(0:(length(Type1)-1))*0.01+B,Type1,col=2,lty=2)
  lines(c(0:(length(Type1)-1))*0.01+B,Type2,col=6,lty=2)
  lines(c(0:(length(Type1)-1))*0.01+B,Fails,col=4,lty=5)
  legend(0.65,0.75,c("Pass","CorrID","Type1","Type2","Fail"),lty=c(1,1,2,2,5),col=c(1,3,2,6,4))
  
  #   
  #   start<-c(0.01,10)/sqrt(0.1^2+10^2)
  #   end<-c(0.3,10)/sqrt(0.3^2+10^2)
  #   astart<-atan(start[2]/start[1])
  #   aend<-atan(end[2]/end[1])
  #   steps<-(astart-aend)/16
  #   par(mfrow=c(4,4))
  #   TMean<-mat.or.vec(16,4)
  #   for (i in 1:16) {
  #     FROSH$V1<-FROSH$T1*sin(astart-steps*i)+(FROSH$SAT_M+FROSH$SAT_V)*cos(astart-steps*i)
  #     TMean[i,1]<-mean(FROSH$V1[which(FROSH$Pass==0)],na.rm=TRUE)
  #     TMean[i,3]<-mean(FROSH$V1[which(FROSH$Pass==1)],na.rm=TRUE)
  #     TMean[i,2]<-sd(FROSH$V1[which(FROSH$Pass==0)],na.rm=TRUE)
  #     TMean[i,4]<-sd(FROSH$V1[which(FROSH$Pass==1)],na.rm=TRUE)
  #     hist(FROSH$V1[which(FROSH$Pass==1)],main=paste("Slope Angle",i),col=4)
  #     hist(FROSH$V1[which(FROSH$Pass==0)],add=T,col=2)
  #     
  #     #     breaks1=unique(quantile(FROSH$V1,c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0),na.rm=TRUE))
  #     #     FROSH$Breaks<-cut(FROSH$V1,breaks=breaks1)  
  #     #     sums2<-by(FROSH[,"Pass"], FROSH$Breaks, sum)
  #     #     lengths2<-by(FROSH[,"Pass"], FROSH$Breaks, length)
  #     #     plot(sums2/lengths2~qnorm(-(0.5/length(sums2))+(1.0/length(sums2))*1:length(sums2)), ylim=c(0,1.0), main = "Proportion WDF", xlab="Quantiles by T(Score)", ylab="Score")
  #     #     fit2 <- glm(cbind(sums2, lengths2-sums2) ~ c(qnorm(-(0.5/length(sums2))+(1.0/length(sums2))*1:length(sums2))), family=binomial)
  #     #     lines(qnorm(c(-0.05+0.1*1:10)), fit2$fitted, type="l", col="red")
  #   }
  #   print(TMean[,3]/TMean[,4]-TMean[,1]/TMean[,2])
  
  
  ###########  20 Question Test Analysis  #################### 
  # SAT_MATH
  par(mfrow=c(2,2))
  CID<-vector(length=0,mode="numeric")
  Type1<-vector(length=0,mode="numeric")
  Type2<-vector(length=0,mode="numeric")
  Stud<-vector(length=0,mode="numeric")
  Fails<-vector(length=0,mode="numeric")
  for (i in 0:length(GoodQ)) {
    Test1<-subset(FROSH,T2>=1/length(GoodQ)*i)
    Stud<-c(Stud,nrow(Test1)/N) # % passing cut
    CID<-c(CID,(length(which(FROSH$Pass==0))-length(which(Test1$Pass==0)))/(N-nrow(Test1)))
    Type1<-c(Type1,(length(which(FROSH$Pass==1))-length(which(Test1$Pass==1)))/(N-nrow(Test1)))
    Type2<-c(Type2,length(which(Test1$Pass==0))/nrow(Test1))
    Fails<-c(Fails,1-length(which(Test1$Pass==0))/WDF)
  }
  plot(c(0:(length(Type1)-1))/(length(Type1)-1),Type1,type="n",main=paste("%",label,"as a function of \n20-Question Score cutoff"),xlab="Cut Grade",ylim=c(0,1.0),ylab=paste("%",label))
  lines(c(0:(length(Type1)-1))/(length(Type1)-1),Stud,col=1,lty=1)
  lines(c(0:(length(Type1)-1))/(length(Type1)-1),CID,col=3,lty=1)
  lines(c(0:(length(Type1)-1))/(length(Type1)-1),Type1,col=2,lty=2)
  lines(c(0:(length(Type1)-1))/(length(Type1)-1),Type2,col=6,lty=2)
  lines(c(0:(length(Type1)-1))/(length(Type1)-1),Fails,col=4,lty=5)
  legend(0.65,0.75,c("Pass","CorrID","Type1","Type2","Fail"),lty=c(1,1,2,2,5),col=c(1,3,2,6,4))
  
  plot(FROSH$SAT_M[which(FROSH$Pass==1)],FROSH$T2[which(FROSH$Pass==1)],col=4,pch=16,main="20 Question Test vs. SAT MATH",xlab="SAT Math Scores",ylab="20 Question Test Score")
  points(FROSH$SAT_M[which(FROSH$Pass==0)],FROSH$T2[which(FROSH$Pass==0)],col=2,pch=16)
  abline(b=-0.17/10,a=1.2)
  abline(b=-0.3/10,a=2.1)
  
  A<--0.03
  B<-2.1
  CID<-vector(length=0,mode="numeric")
  Type1<-vector(length=0,mode="numeric")
  Type2<-vector(length=0,mode="numeric")
  Stud<-vector(length=0,mode="numeric")
  Fails<-vector(length=0,mode="numeric")
  Fails<-vector(length=0,mode="numeric")
  for (i in 0:60) {
    Test1<-subset(FROSH,T2>=SAT_M*A+B+0.01*i)
    Stud<-c(Stud,nrow(Test1)/N) # % passing cut
    CID<-c(CID,(length(which(FROSH$Pass==0))-length(which(Test1$Pass==0)))/(N-nrow(Test1)))
    Type1<-c(Type1,(length(which(FROSH$Pass==1))-length(which(Test1$Pass==1)))/(N-nrow(Test1)))
    Type2<-c(Type2,length(which(Test1$Pass==0))/nrow(Test1))
    Fails<-c(Fails,1-length(which(Test1$Pass==0))/WDF)
  }
  plot(c(0:(length(Type1)-1))*0.01+B,Type1,type="n",main=paste("%",label,"as a function of \n20-Question Score cutoff"),xlab=paste("Cut Intercept, Slope =",A),ylim=c(0,1.0),ylab=paste("%",label))
  lines(c(0:(length(Type1)-1))*0.01+B,Stud,col=1,lty=1)
  lines(c(0:(length(Type1)-1))*0.01+B,CID,col=3,lty=1)
  lines(c(0:(length(Type1)-1))*0.01+B,Type1,col=2,lty=2)
  lines(c(0:(length(Type1)-1))*0.01+B,Type2,col=6,lty=2)
  lines(c(0:(length(Type1)-1))*0.01+B,Fails,col=4,lty=5)
  legend(0.65,0.75,c("Pass","CorrID","Type1","Type2","Fail"),lty=c(1,1,2,2,5),col=c(1,3,2,6,4))
  
  A<--0.017
  B<-1.2
  CID<-vector(length=0,mode="numeric")
  Type1<-vector(length=0,mode="numeric")
  Type2<-vector(length=0,mode="numeric")
  Stud<-vector(length=0,mode="numeric")
  Fails<-vector(length=0,mode="numeric")
  for (i in 0:60) {
    Test1<-subset(FROSH,T2>=SAT_M*A+B+0.01*i)
    Stud<-c(Stud,nrow(Test1)/N) # % passing cut
    CID<-c(CID,(length(which(FROSH$Pass==0))-length(which(Test1$Pass==0)))/(N-nrow(Test1)))
    Type1<-c(Type1,(length(which(FROSH$Pass==1))-length(which(Test1$Pass==1)))/(N-nrow(Test1)))
    Type2<-c(Type2,length(which(Test1$Pass==0))/nrow(Test1))
    Fails<-c(Fails,1-length(which(Test1$Pass==0))/WDF)
  }
  plot(c(0:(length(Type1)-1))*0.01+B,Type1,type="n",main=paste("%",label,"as a function of \n20-Question Score cutoff"),xlab=paste("Cut Intercept, Slope =",A),ylim=c(0,1.0),ylab=paste("%",label))
  lines(c(0:(length(Type1)-1))*0.01+B,Stud,col=1,lty=1)
  lines(c(0:(length(Type1)-1))*0.01+B,CID,col=3,lty=1)
  lines(c(0:(length(Type1)-1))*0.01+B,Type1,col=2,lty=2)
  lines(c(0:(length(Type1)-1))*0.01+B,Type2,col=6,lty=2)
  lines(c(0:(length(Type1)-1))*0.01+B,Fails,col=4,lty=5)
  legend(0.65,0.75,c("Pass","CorrID","Type1","Type2","Fail"),lty=c(1,1,2,2,5),col=c(1,3,2,6,4))
  
  
  # SAT_TOTAL
  par(mfrow=c(2,2))
  CID<-vector(length=0,mode="numeric")
  Type1<-vector(length=0,mode="numeric")
  Type2<-vector(length=0,mode="numeric")
  Stud<-vector(length=0,mode="numeric")
  Fails<-vector(length=0,mode="numeric")
  for (i in 0:length(GoodQ)) {
    Test1<-subset(FROSH,T2>=1/length(GoodQ)*i)
    Stud<-c(Stud,nrow(Test1)/N) # % passing cut
    CID<-c(CID,(length(which(FROSH$Pass==0))-length(which(Test1$Pass==0)))/(N-nrow(Test1)))
    Type1<-c(Type1,(length(which(FROSH$Pass==1))-length(which(Test1$Pass==1)))/(N-nrow(Test1)))
    Type2<-c(Type2,length(which(Test1$Pass==0))/nrow(Test1))
    Fails<-c(Fails,1-length(which(Test1$Pass==0))/WDF)
  }
  plot(c(0:(length(Type1)-1))/(length(Type1)-1),Type1,type="n",main=paste("%",label,"as a function of \n20-Question Score cutoff"),xlab="Cut Grade",ylim=c(0,1.0),ylab=paste("%",label))
  lines(c(0:(length(Type1)-1))/(length(Type1)-1),Stud,col=1,lty=1)
  lines(c(0:(length(Type1)-1))/(length(Type1)-1),CID,col=3,lty=1)
  lines(c(0:(length(Type1)-1))/(length(Type1)-1),Type1,col=2,lty=2)
  lines(c(0:(length(Type1)-1))/(length(Type1)-1),Type2,col=6,lty=2)
  lines(c(0:(length(Type1)-1))/(length(Type1)-1),Fails,col=4,lty=5)
  legend(0.65,0.75,c("Pass","CorrID","Type1","Type2","Fail"),lty=c(1,1,2,2,5),col=c(1,3,2,6,4))
  
  plot((FROSH$SAT_M[which(FROSH$Pass==1)]+FROSH$SAT_V[which(FROSH$Pass==1)]),FROSH$T2[which(FROSH$Pass==1)],col=4,pch=16,main="20 Question Test vs. SAT TOTAL",xlab="SAT Total Scores",ylab="20 Question Test Score")
  points((FROSH$SAT_M[which(FROSH$Pass==0)]+FROSH$SAT_V[which(FROSH$Pass==0)]),FROSH$T2[which(FROSH$Pass==0)],col=2,pch=16)
  abline(b=-0.17/10,a=2.4)
  abline(b=-0.3/10,a=3.9)
  
  A<--0.03
  B<-3.7
  CID<-vector(length=0,mode="numeric")
  Type1<-vector(length=0,mode="numeric")
  Type2<-vector(length=0,mode="numeric")
  Stud<-vector(length=0,mode="numeric")
  Fails<-vector(length=0,mode="numeric")
  for (i in 0:60) {
    Test1<-subset(FROSH,T2>=(SAT_M+SAT_V)*A+B+0.01*i)
    Stud<-c(Stud,nrow(Test1)/N) # % passing cut
    CID<-c(CID,(length(which(FROSH$Pass==0))-length(which(Test1$Pass==0)))/(N-nrow(Test1)))
    Type1<-c(Type1,(length(which(FROSH$Pass==1))-length(which(Test1$Pass==1)))/(N-nrow(Test1)))
    Type2<-c(Type2,length(which(Test1$Pass==0))/nrow(Test1))
    Fails<-c(Fails,1-length(which(Test1$Pass==0))/WDF)
  }
  plot(c(0:(length(Type1)-1))*0.01+B,Type1,type="n",main=paste("%",label,"as a function of \n20-Question Score cutoff"),xlab=paste("Cut Intercept, Slope =",A),ylim=c(0,1.0),ylab=paste("%",label))
  lines(c(0:(length(Type1)-1))*0.01+B,Stud,col=1,lty=1)
  lines(c(0:(length(Type1)-1))*0.01+B,CID,col=3,lty=1)
  lines(c(0:(length(Type1)-1))*0.01+B,Type1,col=2,lty=2)
  lines(c(0:(length(Type1)-1))*0.01+B,Type2,col=6,lty=2)
  lines(c(0:(length(Type1)-1))*0.01+B,Fails,col=4,lty=5)
  legend(0.65,0.75,c("Pass","CorrID","Type1","Type2","Fail"),lty=c(1,1,2,2,5),col=c(1,3,2,6,4))
  
  A<--0.017
  B<-2.1
  CID<-vector(length=0,mode="numeric")
  Type1<-vector(length=0,mode="numeric")
  Type2<-vector(length=0,mode="numeric")
  Stud<-vector(length=0,mode="numeric")
  Fails<-vector(length=0,mode="numeric")
  for (i in 0:60) {
    Test1<-subset(FROSH,T2>=(SAT_M+SAT_V)*A+B+0.01*i)
    Stud<-c(Stud,nrow(Test1)/N) # % passing cut
    CID<-c(CID,(length(which(FROSH$Pass==0))-length(which(Test1$Pass==0)))/(N-nrow(Test1)))
    Type1<-c(Type1,(length(which(FROSH$Pass==1))-length(which(Test1$Pass==1)))/(N-nrow(Test1)))
    Type2<-c(Type2,length(which(Test1$Pass==0))/nrow(Test1))
    Fails<-c(Fails,1-length(which(Test1$Pass==0))/WDF)
  }
  plot(c(0:(length(Type1)-1))*0.01+B,Type1,type="n",main=paste("%",label,"as a function of \n20-Question Score cutoff"),xlab=paste("Cut Intercept, Slope =",A),ylim=c(0,1.0),ylab=paste("%",label))
  lines(c(0:(length(Type1)-1))*0.01+B,Stud,col=1,lty=1)
  lines(c(0:(length(Type1)-1))*0.01+B,CID,col=3,lty=1)
  lines(c(0:(length(Type1)-1))*0.01+B,Type1,col=2,lty=2)
  lines(c(0:(length(Type1)-1))*0.01+B,Type2,col=6,lty=2)
  lines(c(0:(length(Type1)-1))*0.01+B,Fails,col=4,lty=5)
  legend(0.65,0.75,c("Pass","CorrID","Type1","Type2","Fail"),lty=c(1,1,2,2,5),col=c(1,3,2,6,4))
  
  #   start<-c(0.01,10)/sqrt(0.1^2+10^2)
  #   end<-c(0.3,10)/sqrt(0.3^2+10^2)
  #   astart<-atan(start[2]/start[1])
  #   aend<-atan(end[2]/end[1])
  #   steps<-(astart-aend)/16
  #   par(mfrow=c(4,4))
  #   TMean<-mat.or.vec(16,4)
  #   for (i in 1:16) {
  #     FROSH$V1<-FROSH$T2*sin(astart-steps*i)+(FROSH$SAT_M+FROSH$SAT_V)*cos(astart-steps*i)
  #     TMean[i,1]<-mean(FROSH$V1[which(FROSH$Pass==0)],na.rm=TRUE)
  #     TMean[i,3]<-mean(FROSH$V1[which(FROSH$Pass==1)],na.rm=TRUE)
  #     TMean[i,2]<-sd(FROSH$V1[which(FROSH$Pass==0)],na.rm=TRUE)
  #     TMean[i,4]<-sd(FROSH$V1[which(FROSH$Pass==1)],na.rm=TRUE)
  #     hist(FROSH$V1[which(FROSH$Pass==1)],main=paste("Slope Angle",i),col=4)
  #     hist(FROSH$V1[which(FROSH$Pass==0)],add=T,col=2)
  # 
  # #     breaks1=unique(quantile(FROSH$V1,c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0),na.rm=TRUE))
  # #     FROSH$Breaks<-cut(FROSH$V1,breaks=breaks1)  
  # #     sums2<-by(FROSH[,"Pass"], FROSH$Breaks, sum)
  # #     lengths2<-by(FROSH[,"Pass"], FROSH$Breaks, length)
  # #     plot(sums2/lengths2~qnorm(-(0.5/length(sums2))+(1.0/length(sums2))*1:length(sums2)), ylim=c(0,1.0), main = "Proportion WDF", xlab="Quantiles by T(Score)", ylab="Score")
  # #     fit2 <- glm(cbind(sums2, lengths2-sums2) ~ c(qnorm(-(0.5/length(sums2))+(1.0/length(sums2))*1:length(sums2))), family=binomial)
  # #     lines(qnorm(c(-0.05+0.1*1:10)), fit2$fitted, type="l", col="red")
  #   }
  #   print(TMean[,3]/TMean[,4]-TMean[,1]/TMean[,2])
  
  #   par(mfrow=c(2,2))
  #   plot(FROSH$SAT_V[which(FROSH$Pass==1)],FROSH$V1[which(FROSH$Pass==1)],col=4)  
  #   points(FROSH$SAT_V[which(FROSH$Pass==0)],FROSH$V1[which(FROSH$Pass==0)],col=2)
  # 
  #   plot(FROSH$SAT_V[which(FROSH$Pass==1)]+FROSH$SAT_M[which(FROSH$Pass==1)],FROSH$T2[which(FROSH$Pass==1)],col=4)  
  #   points(FROSH$SAT_V[which(FROSH$Pass==0)]+FROSH$SAT_M[which(FROSH$Pass==0)],FROSH$T2[which(FROSH$Pass==0)],col=2)
  # 
  #   
  #   plot(FROSH$T2~FROSH$SAT_V)
  #   plot(FROSH$T2~FROSH$SAT_V)
  
}


PredictPassing <- function(yearcut=17,classcut=200,gradesplit="B",IncludeW=FALSE) {
  #   gradesplit="C"
  #   IncludeW=TRUE
  GradeList<-c("A","B+","B","C+","C","D","F","W")
  #label<-"Target"
  if (IncludeW==FALSE) { 
    GradeList<-GradeList[-length(GradeList)] 
  }
  par(mfrow=c(2,2))
  plot(as.factor(IDF$LetterGrade115),main="115 Grade Distribution")
  plot(as.factor(IDF$LetterGrade123),main="123 Grade Distribution")
  plot(as.factor(IDF$LetterGrade159),main="159 Grade Distribution")
  plot(as.factor(IDF$LetterGrade227),main="227 Grade Distribution")
  par(mfrow=c(1,2))
  plot(IDF$Grade159,IDF$Grade123)
  plot(IDF$Grade159,IDF$Grade115)
  #   print("Rows are 123, Columns are 159") 
  #   table(IDF$LetterGrade123,IDF$LetterGrade159)
  GCut<-match(gradesplit,GradeList)
  Win<-GradeList[1:GCut]
  print(paste("Passing Grades are:",Win))
  Lose<-GradeList[(GCut+1):length(GradeList)]
  label<-""
  for (i in 1:GCut) {label<-paste(label,GradeList[i],sep="")}
  Predictors<-c("MHI","SAT_M","SAT_V","GENDER_CD","ETHNIC_CD","Score.x",paste("C",1:nrow(TestX),sep=""))
  Predict<-c("MHI","SAT_M","SAT_V","Score.x","FCI.x","MATH_REASON.x","CCI.x","MBL.x","ACS.x",paste("C",1:nrow(TestX),sep=""))  
  Q<-paste("C",1:nrow(TestX),sep="")
  QDisc<-data.frame(Q=1:nrow(TestX),Fr=rep(0,nrow(TestX)),Meet=rep(0,nrow(TestX)),G115=rep(0,nrow(TestX)),G123=rep(0,nrow(TestX)),G159=rep(0,nrow(TestX)),G227=rep(0,nrow(TestX)))
  
  if (classcut==100) {  #### Separate analysis for 123 OR 159 with A, B+, B 
    Meet<-IDF[which(!is.na(IDF$Grade123)&!is.na(IDF$Grade159)),c(Predictors,"Grade159","Grade123","FCI.x","MATH_REASON.x","CCI.x","MBL.x","ACS.x","LetterGrade115") %in% colnames(IDF)]
    Meet$Pass<-(Meet$Grade123>=67|Meet$Grade159>=315) + 0
    for (i in 1:nrow(TestX)) {
      test<-Meet[which(!is.na(Meet[,Q[i]])),c("Pass",Q[i])]
      if (nrow(test)/length(Meet$Pass)>0.3) { QDisc[i,"Meet"]<-mean(test[which(test$Pass==1),Q[i]],na.rm=T)-mean(test[which(test$Pass==0),Q[i]],na.rm=T) }
    }
    QS<-QDisc[order(-QDisc$Meet),c("Q","Meet")]    
    FROSH<-Meet
    label<-"Target"
    title<-"Target"
    print(Meet[,paste("C",QS[1:20,1],sep="")])
    
  }
  else if (classcut==200) {  #### Separate analysis for 123 AND 159 passing with A, B+, B
    Meet<-IDF[which(!is.na(IDF$Grade123)&!is.na(IDF$Grade159)),c(Predictors,"Grade159","Grade123","FCI.x","MATH_REASON.x","CCI.x","MBL.x","ACS.x","LetterGrade115") %in% colnames(IDF)]
    Meet$Pass<-(Meet$Grade123>=73&Meet$Grade159>=336) + 0
    for (i in 1:nrow(TestX)) {
      test<-Meet[which(!is.na(Meet[,Q[i]])),c("Pass",Q[i])]
      if (nrow(test)/length(Meet$Pass)>0.3) { QDisc[i,"Meet"]<-mean(test[which(test$Pass==1),Q[i]],na.rm=T)-mean(test[which(test$Pass==0),Q[i]],na.rm=T) }
    }
    QS<-QDisc[order(-QDisc$Meet),c("Q","Meet")]    
    FROSH<-Meet
    label<-"A B+ B"
    title<-"A B+ B"
    label0<-"C+ & Lower"
    print(Meet[,paste("C",QS[1:20,1],sep="")])    
  }
  else if (classcut==115) {
    G115<-IDF[which(!is.na(IDF$LetterGrade115)),c(Predictors,"FCI.x","MATH_REASON.x","CCI.x","MBL.x","ACS.x","LetterGrade115") %in% colnames(IDF)]
    G115$Pass<-(G115$LetterGrade115 %in% Win) + 0
    for (i in 1:nrow(TestX)) {
      test115<-G115[which(!is.na(G115[,Q[i]])),c("Pass",Q[i])]
      if (nrow(test115)/length(G115$Pass)>0.25) { QDisc[i,"G115"]<-mean(test115[which(test115$Pass==1),Q[i]],na.rm=T)-mean(test115[which(test115$Pass==0),Q[i]],na.rm=T) }
    }
    QS<-QDisc[order(-QDisc$G115),c("Q","G115")]    
    FROSH<-G115
    title<-"115"
    label<-paste(label,title)
    print(G115[,paste("C",QS[1:20,1],sep="")])
  }
  else if (classcut==123) {
    G123<-IDF[which(!is.na(IDF$LetterGrade123)),c(Predictors,"FCI.x","MATH_REASON.x","CCI.x","MBL.x","ACS.x","LetterGrade123") %in% colnames(IDF)]
    G123$Pass<-(G123$LetterGrade123 %in% Win) + 0
    for (i in 1:nrow(TestX)) {
      test123<-G123[which(!is.na(G123[,Q[i]])),c("Pass",Q[i])]
      if (nrow(test123)/length(G123$Pass)>0.25) { QDisc[i,"G123"]<-mean(test123[which(test123$Pass==1),Q[i]],na.rm=T)-mean(test123[which(test123$Pass==0),Q[i]],na.rm=T) }
    }
    QS<-QDisc[order(-QDisc$G123),c("Q","G123")]    
    FROSH<-G123
    title<-"123"
    label<-paste(label,title)
  }
  else if (classcut==159) {
    G159<-IDF[which(!is.na(IDF$LetterGrade159)),c(Predictors,"FCI.x","MATH_REASON.x","CCI.x","MBL.x","ACS.x","LetterGrade159") %in% colnames(IDF)]
    G159$Pass<-(G159$LetterGrade159 %in% Win) + 0
    for (i in 1:nrow(TestX)) {
      test159<-G159[which(!is.na(G159[,Q[i]])),c("Pass",Q[i])]
      if (nrow(test159)/length(G159$Pass)>0.25) { QDisc[i,"G159"]<-mean(test159[which(test159$Pass==1),Q[i]],na.rm=T)-mean(test159[which(test159$Pass==0),Q[i]],na.rm=T) }
    }
    QS<-QDisc[order(-QDisc$G159),c("Q","G159")]
    print(G159[,paste("C",QS[1:20,1],sep="")])
    FROSH<-G159
    title<-"159"
    label<-paste(label,title)
  }
  else if (classcut==227) {
    G227<-IDF[which(!is.na(IDF$LetterGrade227)),c(Predictors,"CSEM.x","MATH_REASON.x","LetterGrade227") %in% colnames(IDF)]
    G227$Pass<-(G227$LetterGrade227 %in% Win) + 0
    for (i in 1:nrow(TestX)) {
      test227<-G227[which(!is.na(G227[,Q[i]])),c("Pass",Q[i])]
      if (nrow(test227)/length(G227$Pass)>0.25) { QDisc[i,"G227"]<-mean(test227[which(test227$Pass==1),Q[i]],na.rm=T)-mean(test227[which(test227$Pass==0),Q[i]],na.rm=T) }
    }
    QS<-QDisc[order(-QDisc$G227),c("Q","G227")]
    FROSH<-G227
    title<-"227"
    label<-paste(label,title)
  }
  else {  
    FRESH<-IDF[which(!is.na(IDF$LetterGrade123)&!is.na(IDF$LetterGrade159)),c(Predictors,"FCI.x","MATH_REASON.x","CCI.x","MBL.x","ACS.x","Grade123","Grade159","LetterGrade123","LetterGrade159") %in% colnames(IDF)]
    FRESH$Pass<-((FRESH$LetterGrade123 %in% Win) + 0)*((FRESH$LetterGrade159 %in% Win) + 0)
    for (i in 1:nrow(TestX)) {
      testF<-FRESH[which(!is.na(FRESH[,Q[i]])),c("Pass",Q[i])]
      if (nrow(testF)>40) { QDisc[i,"Fr"]<-mean(testF[which(testF$Pass==1),Q[i]],na.rm=T)-mean(testF[which(testF$Pass==0),Q[i]],na.rm=T) }
    }
    QS<-QDisc[order(-QDisc$Fr),c("Q","Fr")]
    FROSH<-FRESH
    title<-"123 & 159"  
    label<-paste(label,title)
  }
  
  print(QS[1:20,])
  BestQ<-paste("C",QS$Q[1:12],sep="")
  GoodQ<-paste("C",QS$Q[1:20],sep="")
  FROSH$T1<-rowMeans(FROSH[,BestQ],na.rm=TRUE)
  FROSH$T2<-rowMeans(FROSH[,GoodQ],na.rm=TRUE)
  
  #BestQ<-paste("C",c(2,6,8,68,71,73,74,83,89,99,116,124),sep="")
  #GoodQ<-paste("C",c(2,6,8,68,71,73,74,83,89,99,116,124,4,76,82,91,110,113,117,123),sep="")
  # Most significant - 73,71
  
  #logitF <- glm(Pass ~ C2+C6+C8+C68+C71+C73+C74+C83+C89+C99+C116+C124, data = test, family = "binomial")
  par(mfrow=c(2,4))
  if (classcut!=227) {
    plot(density(FROSH$CCI.x[which(FROSH$Pass==0)],na.rm=TRUE),main=paste("CCI Scores \nBlue - Passed",label,"\n Red -",label0),col=2,xlab="CCI Pretest Score",xlim=c(0,1.0))
    lines(density(FROSH$CCI.x[which(FROSH$Pass==1)],na.rm=TRUE),col=4)
  }
  else {
    plot(density(FROSH$CSEM.x[which(FROSH$Pass==0)],na.rm=TRUE),main=paste("CSEM Scores \nBlue - Passed",title,"\n Red -",label),col=2,xlab="CSEM Pretest Score",xlim=c(0,1.0))
    lines(density(FROSH$CSEM.x[which(FROSH$Pass==1)],na.rm=TRUE),col=4)
  }    
  plot(density(FROSH$FCI.x[which(FROSH$Pass==0)],na.rm=TRUE),main=paste("FCI Scores \nBlue - Passed",label,"\n Red -",label0),col=2,xlab="FCI Pretest Score",xlim=c(0,1.0))
  lines(density(FROSH$FCI.x[which(FROSH$Pass==1)],na.rm=TRUE),col=4)
  plot(density(FROSH$T1[which(FROSH$Pass==0)],na.rm=TRUE),main=paste("12 Question Test \nBlue - Passed",label,"\n Red -",label0),col=2,xlab="12Q Score",xlim=c(0,1.0))
  lines(density(FROSH$T1[which(FROSH$Pass==1)],na.rm=TRUE),col=4)
  plot(density(FROSH$T2[which(FROSH$Pass==0)],na.rm=TRUE),main=paste("20 Question Test \nBlue - Passed",label,"\n Red -",label0),col=2,xlab="20Q Score",xlim=c(0,1.0))
  lines(density(FROSH$T2[which(FROSH$Pass==1)],na.rm=TRUE),col=4)
  if (classcut!=227) {
    hist(FROSH$CCI.x[which(FROSH$Pass==0)],breaks=c(0:10*0.1),main=paste("CCI Scores \nBlue - Passed",label,"\n Red -",label0),col=2,xlab="CCI Score",xlim=c(0,1.0))
    hist(FROSH$CCI.x[which(FROSH$Pass==1)],breaks=c(0:10*0.1),add=T,col=4)
  }
  else {
    hist(FROSH$CSEM.x[which(FROSH$Pass==0)],breaks=c(0:10*0.1),main=paste("CSEM Scores \nBlue - Passed",label,"\n Red -",label0),col=2,xlab="CSEM Score",xlim=c(0,1.0))
    hist(FROSH$CSEM.x[which(FROSH$Pass==1)],add=T,col=4,breaks=c(0:10*0.1))
  }
  hist(FROSH$FCI.x[which(FROSH$Pass==0)],breaks=c(0:10*0.1),main=paste("FCI Scores \nBlue - Passed",label,"\n Red -",label0),col=2,xlab="FCI Score",xlim=c(0,1.0))
  hist(FROSH$FCI.x[which(FROSH$Pass==1)],add=T,breaks=c(0:10*0.1),col=4)
  hist(FROSH$T1[which(FROSH$Pass==0)],breaks=c(0:10*0.1),main=paste("12 Question Test \nBlue - Passed",label,"\n Red -",label0),col=2,xlab="12Q Score",xlim=c(0,1.0))
  hist(FROSH$T1[which(FROSH$Pass==1)],add=T,breaks=c(0:10*0.1),col=4)
  hist(FROSH$T2[which(FROSH$Pass==0)],breaks=c(0:10*0.1),main=paste("20 Question Test \nBlue - Passed",label,"\n Red -",label0),col=2,xlab="20Q Score",xlim=c(0,1.0))
  hist(FROSH$T2[which(FROSH$Pass==1)],add=T,breaks=c(0:10*0.1),col=4)
  
  ####################################################################
  ### First graph using Pp/NP as Correct ID
  ### "Slackers" are Type 1 Error (Pass test cut but underachieve)
  ### "Movers" are the Type 2 Error (Fail test but overachieve)
  par(mfrow=c(2,2))
  TestList<-c("BestQ","GoodQ","FCTest")
  TestLength<-c(length(BestQ),length(GoodQ),nrow(FCTest))
  testname<-c("12-Question","20-Question","FCI")
  for (ii in 1:length(TestList)) {
    test<-TestList[ii]
    if (test=="BestQ") {
      FRESH<-FROSH[which(!is.na(FROSH$T1)),]
    }
    if (test=="GoodQ") {
      FRESH<-FROSH[which(!is.na(FROSH$T2)),]
    }
    if (test=="FCTest") {
      FRESH<-FROSH[which(!is.na(FROSH$FCI.x)),]
    }    
    N<-nrow(FRESH)
    Passed<-length(which(FRESH$Pass==1))
    Lower<-length(which(FRESH$Pass==0))
    print(paste("Total % of students who got less than",label,":",format(Lower/N*100,digits=5)))
    CID<-vector(length=0,mode="numeric")
    Type1<-vector(length=0,mode="numeric")
    Type2<-vector(length=0,mode="numeric")
    Stud<-vector(length=0,mode="numeric")
    Fails<-vector(length=0,mode="numeric")
    for (i in 0:TestLength[ii]) {
      if (test=="BestQ") {
        PTest<-subset(FRESH,T1>=1*i/length(BestQ))
        FTest<-subset(FRESH,T1<1*i/length(BestQ))
      }
      if (test=="GoodQ") {
        PTest<-subset(FRESH,T2>=1/length(GoodQ)*i)
        FTest<-subset(FRESH,T2<1/length(GoodQ)*i)
      }
      if (test=="FCTest") {
        PTest<-subset(FRESH,FCI.x>=1/nrow(FCTest)*i)
        FTest<-subset(FRESH,FCI.x<1/nrow(FCTest)*i)
      }
      Stud<-c(Stud,nrow(PTest)/N) # % That pass highest cut, gets lower as cut goes up
      CID<-c(CID,length(which(PTest$Pass==1))/Passed)
      Type1<-c(Type1,length(which(PTest$Pass==0))/nrow(PTest)) 
      Type2<-c(Type2,length(which(FTest$Pass==1))/nrow(FTest)) 
      Fails<-c(Fails,1-(length(which(PTest$Pass==1))+length(which(FTest$Pass==0)))/N)
    }
    if (test=="BestQ") {
      ds0<-density(FRESH$T1[which(FRESH$Pass==0)])
      ds0$y<-ds0$y/range(ds0$y)[2]
      ds1<-density(FRESH$T1[which(FRESH$Pass==1)])
      ds1$y<-ds1$y/range(ds1$y)[2]
    }
    if (test=="GoodQ") {
      ds0<-density(FRESH$T2[which(FRESH$Pass==0)])
      ds0$y<-ds0$y/range(ds0$y)[2]
      ds1<-density(FRESH$T2[which(FRESH$Pass==1)])
      ds1$y<-ds1$y/range(ds1$y)[2]
    }
    if (test=="FCTest") {
      ds0<-density(FRESH$FCI.x[which(FRESH$Pass==0)])
      ds0$y<-ds0$y/range(ds0$y)[2]
      ds1<-density(FRESH$FCI.x[which(FRESH$Pass==1)])
      ds1$y<-ds1$y/range(ds1$y)[2]
    }
    plot(ds1, type="n",xlim=c(0,1),main=paste(label,"group ID as a function of \n",testname[ii],"cutoff Score"),xlab="Cut Grade",ylim=c(0,1.0),ylab=paste("%",label))
    polygon(ds0, col="lightpink", border="red")
    polygon(ds1, col="lightcyan", border="blue",lwd=3)
    lines(ds0, col="red",lwd=3)
    lines(c(0:(length(Type1)-1))/(length(Type1)-1),Stud,col=1,lty=1)
    lines(c(0:(length(Type1)-1))/(length(Type1)-1),CID,col=3,lty=1)
    lines(c(0:(length(Type1)-1))/(length(Type1)-1),Type1,col=2,lty=2)
    lines(c(0:(length(Type1)-1))/(length(Type1)-1),Type2,col=6,lty=2)
    lines(c(0:(length(Type1)-1))/(length(Type1)-1),Fails,col=4,lty=5)
    legend(0.75,0.95,c("PassCut",label,"Slacker","Movers","Fail"),lty=c(1,1,2,2,5),col=c(1,3,2,6,4))
  }  
  
  if (classcut!=227) {
    FRESH<-FROSH[which(!is.na(FROSH$CCI.x)),]
    testname<-"CCI"
  } else {
    FRESH<-FROSH[which(!is.na(FROSH$CSEM.x)),]    
    testname<-"CSEM"
  }
  N<-nrow(FRESH)
  Passed<-length(which(FRESH$Pass==1))
  Lower<-length(which(FRESH$Pass==0))
  print(paste("Total % of students who got less than",label,":",format(Lower/N*100,digits=5)))  
  CID<-vector(length=0,mode="numeric")
  Type1<-vector(length=0,mode="numeric")
  Type2<-vector(length=0,mode="numeric")
  Stud<-vector(length=0,mode="numeric")
  Fails<-vector(length=0,mode="numeric")
  if (classcut!=227) {
    for (i in 0:length(which(TestX$Category=="CCI"))) {
      PTest<-subset(FRESH,CCI.x>=1/length(which(TestX$Category=="CCI"))*i)
      FTest<-subset(FRESH,CCI.x<1/length(which(TestX$Category=="CCI"))*i)
      Stud<-c(Stud,nrow(PTest)/N) # % That pass highest cut, gets lower as cut goes up
      CID<-c(CID,length(which(PTest$Pass==1))/Passed)
      Type1<-c(Type1,length(which(PTest$Pass==0))/nrow(PTest)) 
      Type2<-c(Type2,length(which(FTest$Pass==1))/nrow(FTest)) 
      Fails<-c(Fails,1-(length(which(PTest$Pass==1))+length(which(FTest$Pass==0)))/N)
    }
    ds0<-density(FRESH$CCI.x[which(FRESH$Pass==0)])
    ds0$y<-ds0$y/range(ds0$y)[2]
    ds1<-density(FRESH$CCI.x[which(FRESH$Pass==1)])
    ds1$y<-ds1$y/range(ds1$y)[2]
  }  else {
    for (i in 0:length(which(TestX$Category=="CSEM"))) {
      PTest<-subset(FROSH,CSEM.x>=1/length(which(TestX$Category=="CSEM"))*i)
      FTest<-subset(FROSH,CSEM.x<1/length(which(TestX$Category=="CSEM"))*i)
      Stud<-c(Stud,nrow(PTest)/N) # % That pass highest cut, gets lower as cut goes up
      CID<-c(CID,length(which(PTest$Pass==1))/Passed)
      Type1<-c(Type1,length(which(PTest$Pass==0))/nrow(PTest)) 
      Type2<-c(Type2,length(which(FTest$Pass==1))/nrow(FTest)) 
      Fails<-c(Fails,1-(length(which(PTest$Pass==1))+length(which(FTest$Pass==0)))/N)
    }
    ds0<-density(FRESH$CSEM.x[which(FRESH$Pass==0)])
    ds0$y<-ds0$y/range(ds0$y)[2]
    ds1<-density(FRESH$CSEM.x[which(FRESH$Pass==1)])
    ds1$y<-ds1$y/range(ds1$y)[2]
  }
  plot(ds1, type="n",xlim=c(0,1),main=paste(label,"group ID as a function of \n",testname,"cutoff Score"),xlab="Cut Grade",ylim=c(0,1.0),ylab=paste("%",label))
  polygon(ds0, col="lightpink", border="red")
  polygon(ds1, col="lightcyan", border="blue",lwd=3)
  lines(ds0, col="red",lwd=3)
  lines(c(0:(length(Type1)-1))/(length(Type1)-1),Stud,col=1,lty=1)
  lines(c(0:(length(Type1)-1))/(length(Type1)-1),CID,col=3,lty=1)
  lines(c(0:(length(Type1)-1))/(length(Type1)-1),Type1,col=2,lty=2)
  lines(c(0:(length(Type1)-1))/(length(Type1)-1),Type2,col=6,lty=2)
  lines(c(0:(length(Type1)-1))/(length(Type1)-1),Fails,col=4,lty=5)
  legend(0.75,0.95,c("PassCut",label,"Slacker","Movers","Fail"),lty=c(1,1,2,2,5),col=c(1,3,2,6,4))
  
  ###########  12 Question Test Analysis  #################### 
  # SAT_MATH & SAT_TOTAL
  par(mfrow=c(2,3))
  TestNames<-c("BestQ","GoodQ")
  TestLength<-c(length(BestQ),length(GoodQ))
  testname<-c("12-Question","20-Question")
  for (j in 1:2) {    
    test<-TestNames[j]
    if (test=="BestQ") {
      FRESH<-FROSH[which(!is.na(FROSH$T1)),]
    } else {
      FRESH<-FROSH[which(!is.na(FROSH$T2)),]
    }      
    N<-nrow(FRESH)
    Passed<-length(which(FRESH$Pass==1))
    Lower<-length(which(FRESH$Pass==0))
    print(paste("Total % of students who got less than",label,":",format(Lower/N*100,digits=5)))
    CID<-vector(length=0,mode="numeric")
    Type1<-vector(length=0,mode="numeric")
    Type2<-vector(length=0,mode="numeric")
    Stud<-vector(length=0,mode="numeric")
    Fails<-vector(length=0,mode="numeric")
    for (i in 0:TestLength[j]) {
      if (test=="BestQ") {
        PTest<-subset(FRESH,T1>=1*i/length(BestQ))
        FTest<-subset(FRESH,T1<1*i/length(BestQ))
      } else {
        PTest<-subset(FRESH,T2>=1*i/length(GoodQ))
        FTest<-subset(FRESH,T2<1*i/length(GoodQ))        
      }
      Stud<-c(Stud,nrow(PTest)/N) # % That pass highest cut, gets lower as cut goes up
      CID<-c(CID,length(which(PTest$Pass==1))/Passed)
      Type1<-c(Type1,length(which(PTest$Pass==0))/nrow(PTest)) 
      Type2<-c(Type2,length(which(FTest$Pass==1))/nrow(FTest)) 
      Fails<-c(Fails,1-(length(which(PTest$Pass==1))+length(which(FTest$Pass==0)))/N)
    }
    if (test=="BestQ") {
      ds0<-density(FRESH$T1[which(FRESH$Pass==0)])
      ds0$y<-ds0$y/range(ds0$y)[2]
      ds1<-density(FRESH$T1[which(FRESH$Pass==1)])
      ds1$y<-ds1$y/range(ds1$y)[2]
    } else {
      ds0<-density(FRESH$T2[which(FRESH$Pass==0)])
      ds0$y<-ds0$y/range(ds0$y)[2]
      ds1<-density(FRESH$T2[which(FRESH$Pass==1)])
      ds1$y<-ds1$y/range(ds1$y)[2]      
    }
    plot(ds1, type="n",xlim=c(0,1),main=paste(label,"group ID as a function of \n",testname[j],"cutoff Score"),xlab="Cut Grade",ylim=c(0,1.0),ylab=paste("%",label))
    polygon(ds0, col="lightpink", border="red")
    polygon(ds1, col="lightcyan", border="blue",lwd=3)
    lines(ds0, col="red",lwd=3)
    lines(c(0:(length(Type1)-1))/(length(Type1)-1),Stud,col=1,lty=1)
    lines(c(0:(length(Type1)-1))/(length(Type1)-1),CID,col=3,lty=1)
    lines(c(0:(length(Type1)-1))/(length(Type1)-1),Type1,col=2,lty=2)
    lines(c(0:(length(Type1)-1))/(length(Type1)-1),Type2,col=6,lty=2)
    lines(c(0:(length(Type1)-1))/(length(Type1)-1),Fails,col=4,lty=5)
    legend(0.75,0.95,c("PassCut",label,"Slacker","Movers","Fail"),lty=c(1,1,2,2,5),col=c(1,3,2,6,4))
  
    if (j==1) {
      plot(FRESH$SAT_M[which(FRESH$Pass==1)],FRESH$T1[which(FRESH$Pass==1)],ylim=c(0,1),col="blue",pch=19,main=paste(testname[j],"Test vs. SAT MATH"),xlab="SAT Math Scores",ylab=paste(testname[j],"Test Score"))
      points(FRESH$SAT_M[which(FRESH$Pass==0)],FRESH$T1[which(FRESH$Pass==0)],col="red",pch=19)

      plot((FRESH$SAT_M[which(FRESH$Pass==1)]+FRESH$SAT_V[which(FRESH$Pass==1)]),FRESH$T1[which(FRESH$Pass==1)],ylim=c(0,1),col="blue",pch=19,main=paste(testname[j],"Test vs. SAT TOTAL"),xlab="SAT Total Scores",ylab=paste(testname[j],"Test Score"))
      points((FRESH$SAT_M[which(FRESH$Pass==0)]+FRESH$SAT_V[which(FRESH$Pass==0)]),FRESH$T1[which(FRESH$Pass==0)],col="red",pch=19)
    } else {
      plot(FRESH$SAT_M[which(FRESH$Pass==1)],FRESH$T2[which(FRESH$Pass==1)],ylim=c(0,1),col="blue",pch=19,main=paste(testname[j],"Test vs. SAT MATH"),xlab="SAT Math Scores",ylab=paste(testname[j],"Test Score"))
      points(FRESH$SAT_M[which(FRESH$Pass==0)],FRESH$T2[which(FRESH$Pass==0)],col="red",pch=19)
      
      plot((FRESH$SAT_M[which(FRESH$Pass==1)]+FRESH$SAT_V[which(FRESH$Pass==1)]),FRESH$T2[which(FRESH$Pass==1)],ylim=c(0,1),col="blue",pch=19,main=paste(testname[j],"Test vs. SAT TOTAL"),xlab="SAT Total Scores",ylab=paste(testname[j],"Test Score"))
      points((FRESH$SAT_M[which(FRESH$Pass==0)]+FRESH$SAT_V[which(FRESH$Pass==0)]),FRESH$T2[which(FRESH$Pass==0)],col="red",pch=19)      
    }
  }    
}

#,paste("QP",1:42,sep=""),paste("QC",1:50,sep="")