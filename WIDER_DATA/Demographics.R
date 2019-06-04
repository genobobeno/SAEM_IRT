############################################################################
# This file is for demographic questions
# Currently the functions here are
#    1. Correlations() <- runs a set of correlations on Gender, Ethnicity, SATs, MHI, etc.
#        Spits out multiple files in your current working directory.
#    2. InfoPlots() <- plots a series of distributions for information purposes

CLASSCorrelations <- function(extracut=NA) {
  if (!is.na(extracut)) {
    IDF<-IDF[eval(parse(text=extracut)),] 
  }

  Ccats<-outer(c("Over","All","PI","RW","Gen","Con","Soph","SE","CC","CL","Per"),c("F","U"),function(x,y) paste(x,y,sep=""))
  Pcats<-outer(c("Over","All","PI","RW","Gen","Con","Soph","SE","CU","AC"),c("F","U"),function(x,y) paste(x,y,sep=""))
  PostCcols<-which(colnames(IDF) %in% paste("cclass",t(Ccats),sep=""))
  PreCcols<-which(colnames(IDF) %in% paste("CCLASS",t(Ccats),sep=""))
  PostPcols<-which(colnames(IDF) %in% paste("pclass",t(Pcats),sep=""))
  PrePcols<-which(colnames(IDF) %in% paste("PCLASS",t(Pcats),sep=""))
  
  DF<-IDF[,PreCcols[1:2]]/100
  DF[,ncol(DF)+1]<-DF[,1]+DF[,2]
  DF[,ncol(DF)+1]<-DF[,1]-DF[,2]
  DF[,ncol(DF)+1]<-log(DF[,1]/DF[,2])
  colnames(DF)[ncol(DF)-2:0]<-paste("C",substr(Ccats[1,1], 1, nchar(Ccats[1,1])-1),c("Opin","Polar","Map"),sep="")
    
  for (i in 2:(length(PreCcols)/2)) {
    DF<-cbind(DF,IDF[,PreCcols[2*i-1:0]]/100)
    C1<-DF[,ncol(DF)-1]+DF[,ncol(DF)]
    C2<-DF[,ncol(DF)-1]-DF[,ncol(DF)]
    C3<-log(DF[,ncol(DF)-1]/DF[,ncol(DF)])
    DF<-cbind(DF,C1,C2,C3)
    colnames(DF)[ncol(DF)-2:0]<-paste("C",substr(Ccats[i,1], 1, nchar(Ccats[i,1])-1),c("Opin","Polar","Map"),sep="")
  }  
  
  for (i in 1:(length(PostCcols)/2)) {
    DF<-cbind(DF,IDF[,PostCcols[2*i-1:0]]/100)
    C1<-DF[,ncol(DF)-1]+DF[,ncol(DF)]
    C2<-DF[,ncol(DF)-1]-DF[,ncol(DF)]
    C3<-log(as.numeric(DF[,ncol(DF)-1]/DF[,ncol(DF)]))
    #C3<-log(DF[,ncol(DF)-1]/DF[,ncol(DF)])
    DF<-cbind(DF,C1,C2,C3)
    colnames(DF)[ncol(DF)-2:0]<-paste("c",substr(Ccats[i,1], 1, nchar(Ccats[i,1])-1),c("Opin","Polar","Map"),sep="")
  }  
  
  for (i in 1:(length(PrePcols)/2)) {
    DF<-cbind(DF,IDF[,PrePcols[2*i-1:0]]/100)
    C1<-DF[,ncol(DF)-1]+DF[,ncol(DF)]
    C2<-DF[,ncol(DF)-1]-DF[,ncol(DF)]
    C3<-log(as.numeric(DF[,ncol(DF)-1]/DF[,ncol(DF)]))
    #C3<-log(DF[,ncol(DF)-1]/DF[,ncol(DF)])
    DF<-cbind(DF,C1,C2,C3)
    colnames(DF)[ncol(DF)-2:0]<-paste("P",substr(Pcats[i,1], 1, nchar(Pcats[i,1])-1),c("Opin","Polar","Map"),sep="")
  }  
  
  for (i in 1:(length(PostPcols)/2)) {
    DF<-cbind(DF,IDF[,PostPcols[2*i-1:0]]/100)
    C1<-DF[,ncol(DF)-1]+DF[,ncol(DF)]
    C2<-DF[,ncol(DF)-1]-DF[,ncol(DF)]
    C3<-log(DF[,ncol(DF)-1]/DF[,ncol(DF)])
    DF<-cbind(DF,C1,C2,C3)
    colnames(DF)[ncol(DF)-2:0]<-paste("p",substr(Pcats[i,1], 1, nchar(Pcats[i,1])-1),c("Opin","Polar","Map"),sep="")
  }  
  #str(IDF[,-c(which(colnames(IDF)=="C1"):ncol(IDF))])[c(4:28)]
  Post<-IDF[,which(grepl(".y",colnames(IDF))==TRUE)]
  DF<-cbind(IDF[,-c(which(colnames(IDF)=="C1"):ncol(IDF))],Post[,1:4],DF)
  DF<-DF[,2:ncol(DF)]  # Removed IDs
  DF<-DF[,-c(1,2,4,5,10,13,16,17,20)] # Removed Factors
  write.csv(cor(DF,use="pairwise.complete.obs"),file = paste(OutDir,"/Corr_CLASS.csv",sep=""))  
#  for (i in 1:5) write.csv(cor(DF[which(DF$SAT_V_Quin==i),c(1,4:22)],use="pairwise.complete.obs"),file = paste(OutDir,"/Corr_SAT_V_Quin_",i,".csv",sep=""))  
#  for (i in 1:5) write.csv(cor(DF[which(DF$SAT_M_Quin==i),c(1,4:22)],use="pairwise.complete.obs"),file = paste(OutDir,"/Corr_SAT_M_Quin_",i,".csv",sep=""))  
}


Correlations <- function(extracut=NA) {
  if (!is.na(extracut)) {
    IDF<-IDF[eval(parse(text=extracut)),] 
  }
  
  DF<-IDF[,c(4:10,12,13,15,16,18:24,281:284)]
  print("Load data into a truncated data.frame using:  COR_DF<-Correlations()  ")
  print("The cuts you can use on COR_DF are: ")  
  colnames(DF)
  DF$GENDER_CD<-as.factor(DF$GENDER_CD)
  DF$ETHNIC_CD<-as.factor(DF$ETHNIC_CD)
  cats <- c(0,27300,53230,85500,132000,Inf)
  DF$MHI_Quin <- findInterval(DF$MHI,cats,all.inside = TRUE)
  DF$SAT_M_Quin<-findInterval(DF$SAT_M,unique(quantile(DF$SAT_M,seq(0,1.0,0.2),na.rm=TRUE)),all.inside=TRUE)
  DF$SAT_V_Quin<-findInterval(DF$SAT_V,unique(quantile(DF$SAT_V,seq(0,1.0,0.2),na.rm=TRUE)),all.inside=TRUE)
  DF$MHI_Quin <- as.factor(DF$MHI_Quin)
  DF$SAT_M_Quin <- as.factor(DF$SAT_M_Quin)
  DF$SAT_V_Quin <- as.factor(DF$SAT_V_Quin)
  
  for (i in 2:5) write.csv(cor(DF[which(DF$MHI_Quin==i),c(1,4:22)],use="pairwise.complete.obs"),file = paste(OutDir,"/Corr_MHI_Quin_",i,".csv",sep=""))  
  for (i in 1:5) write.csv(cor(DF[which(DF$SAT_V_Quin==i),c(1,4:22)],use="pairwise.complete.obs"),file = paste(OutDir,"/Corr_SAT_V_Quin_",i,".csv",sep=""))  
  for (i in 1:5) write.csv(cor(DF[which(DF$SAT_M_Quin==i),c(1,4:22)],use="pairwise.complete.obs"),file = paste(OutDir,"/Corr_SAT_M_Quin_",i,".csv",sep=""))  
  print("Wrote Correlation files into your Output Directory by MHI and SAT quintiles")  
  return(DF)
  # ALL_C<-cor(DF[,c(9,10,15:20)],use="pairwise.complete.obs")
  # ALL_M<-cor(IDF[which(IDF$GENDER_CD=="M"),c(9,10,15:20)],use="pairwise.complete.obs")
  # ALL_F<-cor(IDF[which(IDF$GENDER_CD=="F"),c(9,10,15:20)],use="pairwise.complete.obs")
  # ETH_3<-cor(IDF[which(IDF$ETHNIC_CD==3),c(9,10,15:20)],use="pairwise.complete.obs")
  # ETH_4<-cor(IDF[which(IDF$ETHNIC_CD==4),c(9,10,15:20)],use="pairwise.complete.obs")
  # ETH_6<-cor(IDF[which(IDF$ETHNIC_CD==6),c(9,10,15:20)],use="pairwise.complete.obs")
  # ETH_8<-cor(IDF[which(IDF$ETHNIC_CD==8),c(9,10,15:20)],use="pairwise.complete.obs")  
  
}




InfoPlots <- function(classcut=123,yearcut=17,gradesplit="C",IncludeW=TRUE) {
  GradeList<-c("A","B+","B","C+","D","F","W")
  plot(IDF$Grade123~IDF$Grade159)
  colSums(table(IDF$LetterGrade123,IDF$LetterGrade159))
  par(mfrow=c(1,2))
  plot(as.factor(IDF$LetterGrade123[which(!is.na(IDF$LetterGrade159)&!is.na(IDF$LetterGrade123))]),main="Grades in 123",ylim=c(0,200))
  plot(as.factor(IDF$LetterGrade159[which(!is.na(IDF$LetterGrade159)&!is.na(IDF$LetterGrade123))]),main="Grades in 159",ylim=c(0,200))
  
  par(mfrow=c(1,3))
  plot(as.factor(IDF$LetterGrade115[which(!is.na(IDF$LetterGrade159)&!is.na(IDF$LetterGrade115))]),main="Grades in 115",ylim=c(0,200))
  plot(as.factor(IDF$LetterGrade123[which(!is.na(IDF$LetterGrade159)&!is.na(IDF$LetterGrade123))]),main="Grades in 123",ylim=c(0,200))
  # Stacked Bar Plot with Colors and Legend
  barplot(rbind(table(IDF$LetterGrade159[which(!is.na(IDF$LetterGrade159)&!is.na(IDF$LetterGrade123))]),table(IDF$LetterGrade159[which(!is.na(IDF$LetterGrade159)&!is.na(IDF$LetterGrade115))])), main="159 Grade Distribution \n by Physics Class",
          col=c("darkblue","red"),
          legend = c("123","115"),ylim=c(0,200))   
  
  table(IDF$LetterGrade123[which(IDF$Grade123<(-0.25*IDF$Grade159+140))],IDF$LetterGrade159[which(IDF$Grade123<(-0.25*IDF$Grade159+140))])
  table(IDF$LetterGrade123[which(IDF$Grade123<67&IDF$Grade159<315)],IDF$LetterGrade159[which(IDF$Grade123<67&IDF$Grade159<315)]) # 87
  
  par(mfrow=c(2,2))
  plot(IDF$Grade123[which(IDF$Grade159>315|IDF$Grade123>67)]~IDF$Grade159[which(IDF$Grade159>315|IDF$Grade123>67)],col=4,ylab="123 Grades",xlab="159 Grades",main="Profile of Course Grades\n with 'C' Cutoff",ylim=c(0,100),xlim=c(0,480))
  points(IDF$Grade123[which(IDF$Grade123<67&IDF$Grade159<315)]~IDF$Grade159[which(IDF$Grade123<67&IDF$Grade159<315)],col=2)
  abline(v=315)
  abline(h=67)
  
  hist(IDF$MHI[which(!is.na(IDF$Grade159*IDF$Grade123))],main="Median Household Income",col=4,xlab="($)",ylab="Students")
  hist(IDF$MHI[which(IDF$Grade123<67&IDF$Grade159<315)],add=T,col=2)
  hist(IDF$SAT_M[which(!is.na(IDF$Grade159*IDF$Grade123))],col=4,main="SAT Math Scores",xlab="Score",ylab="Students")
  hist(IDF$SAT_M[which(IDF$Grade123<67&IDF$Grade159<315)],add=T,col=2)
  hist(IDF$SAT_V[which(!is.na(IDF$Grade159*IDF$Grade123))],col=4,main="SAT Verbal Scores",xlab="Score",ylab="Students")
  hist(IDF$SAT_V[which(IDF$Grade123<67&IDF$Grade159<315)],add=T,col=2)
  
  plot(IDF$ETHNIC_CD[which(!is.na(IDF$Grade159*IDF$Grade123))])
  plot(IDF$ETHNIC_CD[which(IDF$Grade123<67&IDF$Grade159<315)]) 
}