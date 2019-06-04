# Set the Working Directory here (Hint: It's the directory where you just found "TestAnalysis.R")
#   You can also set the Working Directory by clicking: Session -> Set Working Directory -> To Source File Location
#    you will then notice it run the command setwd() in the bottom left window.

source("WIDER_DATA/InitiateData.R")

print("Functions you can use:")
print("     Anonymize()")
print("      ================================================================================")
print("         Function: Run it by typing Anonymize() and hitting Enter. It will create")
print("                   a file called AnonItemData.csv in your current working directory.")
print("                   This file will have no RUIDs or names attached to the rows of data.")
print("")
print("     ItemStats(ItemNumber, classcut=int, yearcut=int, schoolcut=int, testcut=char, formcut=int)")
print("      ================================================================================")
print("         Function: Will display percent correct and standard deviation for blanks=0 and blanks not included.")
print("                   Cuts shown are an example and are not necessary as input arguments. When including them,")
print("                   you must write it out in the format shown. Examples of use are shown here...")
print("")
print("         Example:  ItemStats(35)      <- will run analysis on item 35 from Master List")
print("                   ***Note:  All cuts are optional and described in ShowCutOptions()")
print("")
print("     CutCategories()")
print("      ================================================================================")
print("         Function: Prints categories that are worth cutting on, be it demographics,")
print("                   SATs, Grades, etc.  Also shows subcategories of categorical ")
print("                   variables that can be cut. ")
print("")
print("     CLASSCategories()")
print("      ================================================================================")
print("         Function: Prints CLASS 'categories' from Excel analysis using CLASS Template.")
print("                   Gives a printout of category and the respective column names for ")
print("                   the pretest and posttest of the physics and chemistry CLASS surveys.")
print("")
print("     ScoreNorm()")
print("      ================================================================================")
print("         Function: Plots Final Grades and Standardized Final Grade Scores with Q-Q.")
print("")
print("     PlotGrades()")
print("      ================================================================================")
print("         Function: Plots Final Grade distributions.")
print("")
# Function to get item statistics

source("ItemFunctions.R")
source("Demographics.R")
source("LatentClass.R")
source("Predictions.R")
source("IRT.R")
source("Diagnostics.R")
source("Options.R")
source("CLASS.R")
source("GFA.R")


colnames(IDF)

PrePostItem(12,classcut=123)
TestX
FCIC<-paste("C",64:93,sep="")
Response<-IDF[which(!is.na(IDF$FCI.x)&(!is.na(IDF$LetterGrade123)|!is.na(IDF$LetterGrade115))),FCIC]
str(Response)
IRTFit<-ApplyIRT(RP=Response,model="ogive")
GFAFit<-ApplyGFA(RP=Response,GParam=FALSE,AParams=1,HDQs=NA,Plot=FALSE)
points(1:60,c(IRTFit$A,IRTFit$B),col=6,cex=2)

MatchPrePost<-function(col1="CSEM.x",col2="CSEM.y",course=227,itemflag=F,extracut=NA) {
  if(!is.na(extracut)) {
    IDF<-IDF[eval(parse(text=extracut))]
  }
  coursename<-paste("Final",course,sep="")
  DF<-IDF[which(!is.na(IDF[,col1]) & !is.na(IDF[,col2]) & !is.na(IDF[,coursename])),]
  print(paste("Number of students:",nrow(DF)))
  if (itemflag==TRUE) {
    Nones1<-sum(DF[,col1])
    Nones2<-sum(DF[,col2])
    NAll<-nrow(DF)
    print(paste(col1,"Mean:",format(mean(DF[,col1]),digits=4)," 95% Conf. Int.:",format(binom.confint(Nones1,NAll, methods = "wilson")$lower,digits=4),
                "-",format(binom.confint(Nones1,NAll, methods = "wilson")$upper,digits=4)))
    print(paste(col2,"Mean:",format(mean(DF[,col2]),digits=4)," 95% Conf. Int.:",format(binom.confint(Nones2,NAll, methods = "wilson")$lower,digits=4),
                "-",format(binom.confint(Nones2,NAll, methods = "wilson")$upper,digits=4)))
  } else {
    print(paste(col1,"Mean:",format(mean(DF[,col1]),digits=4),"  SD:",format(sd(DF[,col1]),digits=4)))
    print(paste(col2,"Mean:",format(mean(DF[,col2]),digits=4),"  SD:",format(sd(DF[,col2]),digits=4)))
  }
}



CLASSCategories()
par(mfrow=c(2,2))
plot(IDF$Grade123,IDF$PCLASSOverF-IDF$PCLASSOverU)
plot(IDF$Grade123,IDF$PCLASSAllF-IDF$PCLASSAllU)
plot(IDF$Grade123,IDF$PCLASSPIF-IDF$PCLASSPIU)
plot(IDF$Grade123,IDF$PCLASSRWF-IDF$PCLASSRWU)
plot(IDF$Grade123,IDF$PCLASSGenF-IDF$PCLASSGenU)
plot(IDF$Grade123,IDF$PCLASSConF-IDF$PCLASSConU)
plot(IDF$Grade123,IDF$PCLASSSophF-IDF$PCLASSSophU)
plot(IDF$Grade123,IDF$PCLASSSEF-IDF$PCLASSSEU)
plot(IDF$Grade123,IDF$PCLASSCUF-IDF$PCLASSCUU)
plot(IDF$Grade123,IDF$PCLASSACF-IDF$PCLASSACU)

cor(IDF$PCLASSOverF-IDF$PCLASSOverU,IDF$Grade123,use="pairwise.complete.obs")
cor(IDF$PCLASSAllF-IDF$PCLASSAllU,IDF$Grade123,use="pairwise.complete.obs")
cor(IDF$PCLASSPIF-IDF$PCLASSPIU,IDF$Grade123,use="pairwise.complete.obs")
cor(IDF$PCLASSRWF-IDF$PCLASSRWU,IDF$Grade123,use="pairwise.complete.obs")
cor(IDF$PCLASSGenF-IDF$PCLASSGenU,IDF$Grade123,use="pairwise.complete.obs")
cor(IDF$PCLASSConF-IDF$PCLASSConU,IDF$Grade123,use="pairwise.complete.obs")
cor(IDF$PCLASSSophF-IDF$PCLASSSophU,IDF$Grade123,use="pairwise.complete.obs")
cor(IDF$PCLASSSEF-IDF$PCLASSSEU,IDF$Grade123,use="pairwise.complete.obs")
cor(IDF$PCLASSCUF-IDF$PCLASSCUU,IDF$Grade123,use="pairwise.complete.obs")
cor(IDF$PCLASSACF-IDF$PCLASSACU,IDF$Grade123,use="pairwise.complete.obs")

colnames(IDF)
PCL<-42
Cvec<-rep(0,PCL)
for(i in 1:PCL) {
Cvec[i]<-cor(IDF[,760+i],IDF$Grade123,use="pairwise.complete.obs")
}
Cvec  # 794, 785
par(mfrow=c(3,3))
PCRows1<-which(IDF[,794]==1&IDF[,785]==1)
PCRows2<-which(IDF[,794]==1&IDF[,785]==0)
PCRows3<-which(IDF[,794]==1&IDF[,785]==-1)
PCRows4<-which(IDF[,794]==0&IDF[,785]==1)
PCRows5<-which(IDF[,794]==0&IDF[,785]==0)
PCRows6<-which(IDF[,794]==0&IDF[,785]==-1)
PCRows7<-which(IDF[,794]==-1&IDF[,785]==1)
PCRows8<-which(IDF[,794]==-1&IDF[,785]==0)
PCRows9<-which(IDF[,794]==-1&IDF[,785]==-1)
BK<-seq(0,100,length.out=11)
hist(IDF$Grade123[PCRows1],breaks=BK)
hist(IDF$Grade123[PCRows2],breaks=BK)
hist(IDF$Grade123[PCRows3],breaks=BK)
hist(IDF$Grade123[PCRows4],breaks=BK)
hist(IDF$Grade123[PCRows5],breaks=BK)
hist(IDF$Grade123[PCRows6],breaks=BK)
hist(IDF$Grade123[PCRows7],breaks=BK)
hist(IDF$Grade123[PCRows8],breaks=BK)
hist(IDF$Grade123[PCRows9],breaks=BK)


DFBad<-IDF[which(IDF$LetterGrade123 %in% c("F","D","W","C") & IDF$LetterGrade159 %in% c("F","D","W","C")),]
DFGood<-IDF[which(IDF$LetterGrade123 %in% c("A","B+") & IDF$LetterGrade159 %in% c("A","B+")),]

DFBad<-IDF[which(IDF$LetterGrade159 %in% c("F","D","W","C")),]
DFGood<-IDF[which(IDF$LetterGrade159 %in% c("A","B+")),]

nrow(DFBad)
nrow(DFGood)

DFBad<-DFBad[order(scale(DFBad$Grade123)+scale(DFBad$Grade159)),]
DFGood<-DFGood[order(-scale(DFGood$Grade123)-scale(DFGood$Grade159)),]

#DFBad<-DFBad[order(scale(DFBad$Grade159)),]
#DFGood<-DFGood[order(-scale(DFGood$Grade159)),]

CCols<-paste("C",1:129,sep="")
QPCols<-paste("QP",1:42,sep="")
QCCols<-paste("QC",1:50,sep="")
Cols<-c(CCols,QPCols,QCCols)

# par(mfrow=c(3,1))
# for (i in 5:min(nrow(DFBad),nrow(DFGood))) {
#   barplot(colSums(DFGood[1:i,Cols],na.rm=TRUE)-colSums(DFBad[1:i,Cols],na.rm=TRUE),main="Difference")
#   barplot(colSums(DFGood[1:i,Cols],na.rm=TRUE),main="Good Students")
#   barplot(colSums(DFBad[1:i,Cols],na.rm=TRUE),main="Bad Students")
#   Sys.sleep(0.2)
# }


PTestL<-40
PreMaxQs<-colSums(DFGood[1:min(nrow(DFBad),nrow(DFGood)),Cols],na.rm=TRUE)-colSums(DFBad[1:min(nrow(DFBad),nrow(DFGood)),Cols],na.rm=TRUE)
PreMaxQs<-PreMaxQs[order(-abs(PreMaxQs))]
PreTestQs<-names(PreMaxQs)[1:PTestL]

ATT<-(grepl("QP",PreTestQs)+0) + (grepl("QC",PreTestQs)+0)
SCI<-(-1*(ATT-1))
sum(SCI)

Trans<-(-1*((PreMaxQs[1:PTestL]<0)+0))+((PreMaxQs[1:PTestL]>0)+0)

DFFrosh<-IDF[which(!is.na(IDF$LetterGrade123) & !is.na(IDF$LetterGrade159)),]
DFFrosh$PSCORE<-NA
DFFrosh$ASCORE<-NA
DFFrosh$SSCORE<-NA
for(i in 1:nrow(DFFrosh)) {
  DFFrosh$PSCORE[i]<-sum(DFFrosh[i,PreTestQs]*(Trans),na.rm=TRUE)/length(which(!is.na(DFFrosh[i,PreTestQs])))
  ATrans<-ATT*Trans
  DFFrosh$ASCORE[i]<-sum(DFFrosh[i,PreTestQs]*(ATrans),na.rm=TRUE)/length(which(!is.na(DFFrosh[i,PreTestQs]*ATT)))
  STrans<-SCI*Trans
  DFFrosh$SSCORE[i]<-sum(DFFrosh[i,PreTestQs]*(STrans),na.rm=TRUE)/length(which(!is.na(DFFrosh[i,PreTestQs]*SCI)))
}

# par(mfrow=c(1,2))
# plot(DFFrosh$PSCORE,DFFrosh$Grade123,main="123 Grade vs. PreTest Score")
# plot(DFFrosh$PSCORE,DFFrosh$Grade159,main="159 Grade vs. PreTest Score")

# cor(DFFrosh$PSCORE,DFFrosh$Grade123,use="pairwise.complete.obs")
# cor(DFFrosh$PSCORE,DFFrosh$Grade159,use="pairwise.complete.obs")

m1<-lm(DFFrosh$Grade123~scale(DFFrosh$PSCORE)+scale(DFFrosh$SAT_M))
summary(m1)
anova(m1)

m1<-lm(DFFrosh$Grade123~scale(DFFrosh$ASCORE)+scale(DFFrosh$SAT_M))
summary(m1)
anova(m1)
m1<-lm(DFFrosh$Grade123~scale(DFFrosh$SSCORE)+scale(DFFrosh$SAT_M))
summary(m1)
anova(m1)
m1<-lm(DFFrosh$Grade123~scale(DFFrosh$ASCORE)+scale(DFFrosh$SSCORE)+scale(DFFrosh$SAT_M))
summary(m1)
anova(m1)
m1<-lm(DFFrosh$Grade123~scale(DFFrosh$ASCORE)+scale(DFFrosh$SSCORE))
summary(m1)
anova(m1)
m1<-lm(DFFrosh$Grade123~scale(DFFrosh$SSCORE))
summary(m1)
anova(m1)
m1<-lm(DFFrosh$Grade123~scale(DFFrosh$ASCORE))
summary(m1)
anova(m1)


plot(m1$model[,2],m1$residuals)
str(m1)


plot(density(IDF$CSEM.y[which(!is.na(IDF$Grade116))],na.rm=TRUE))
# lines(density(IDF$CSEM.y[which(is.na(IDF$Grade116))],na.rm=TRUE),col=2)
# plot(density(IDF$CSEM.x[which(!is.na(IDF$Grade116))],na.rm=TRUE))
# lines(density(IDF$CSEM.x[which(is.na(IDF$Grade116))],na.rm=TRUE),col=2)

# Code for running on each class
# IncludeW=TRUE
# gradesplit="C" # This means we're cutting at Grade X >= C & X < C {A, B, C+ are all "greater" than C, while D, F, & W are "less than"}
# GradeList<-c("A","B+","B","C+","C","D","F","W")
# if (IncludeW==FALSE) { 
#   GradeList<-GradeList[-length(GradeList)] 
# }
# GCut<-match(gradesplit,GradeList)
# Win<-GradeList[1:GCut]
# Lose<-GradeList[(GCut+1):length(GradeList)]
# IDF$Pass115<-(IDF$LetterGrade115 %in% Win) + 0
# IDF$Pass123<-(IDF$LetterGrade123 %in% Win) + 0
# IDF$Pass159<-(IDF$LetterGrade159 %in% Win) + 0
# IDF$Pass227<-(IDF$LetterGrade227 %in% Win) + 0
# IDF$Pass.123.159<-((IDF$LetterGrade123 %in% Win) + 0)*((IDF$LetterGrade159 %in% Win) + 0)



# Now you can cut on 


# DF227<-read.csv("../SuzanneCode/CLASS227PrePost.csv",header=T,stringsAsFactors=F,as.is=T,sep=",")
# CLASSCategories()
# which(colnames(IDF)==c("PCLASSOverF","pclassOverF"))
# ncol(DF227)
# nrow(IDF[which(!is.na(IDF$Grade227)&!is.na(IDF$LetterGrade116)&!is.na(IDF$PCLASSOverF)),403:422])
# length(IDF[which(!is.na(IDF$Grade227)&is.na(IDF$LetterGrade116)),804:823])
# DFN116<-DF227
# DFN116[1,2:21]<-colMeans(IDF[which(!is.na(IDF$Grade227)&is.na(IDF$LetterGrade116)),403:422],na.rm=TRUE)
# DFN116[2,2:21]<-colMeans(IDF[which(!is.na(IDF$Grade227)&is.na(IDF$LetterGrade116)),804:823],na.rm=TRUE)
# write.csv(DFN116, file = "../SuzanneCode/CLASS227N116PrePost.csv",row.names=FALSE)
# DF227
# DF116