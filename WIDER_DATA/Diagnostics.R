##########################################
# A series of exploratory operations for data analysis
# Currently being worked on... Brainstorming proper work flow...
#    1. ScoreNorm()        <- Normalizes scores for the requested test
#    2. Anonymize()        <- Spits out file without RUIDs or names
#    3. PlotGrades()       <- Plots grade distribution of each class
#    4. CutCategories()    <- Outputs the factor/categorical variables that can be used for demographics
#    5. CLASSCategories()  <- Outputs the CLASS categories you can call by name


ScoreNorm <- function() {
  par(mfrow=c(4,3))
  Z115<-qqnorm(IDF$Grade115,main="Q-Q 115 Grades")$x
  hist(Z115,main=paste("115 Final Grade (Std.), N(",format(mean(IDF$Grade115,na.rm=TRUE),digits=4),",sd=",format(sd(IDF$Grade115,na.rm=TRUE),digits=4),")",sep=""))
  plot(IDF$LetterGrade115, main = paste("Physics 115 Grades; N = ",length(IDF$LetterGrade115[which(!is.na(IDF$LetterGrade115))])))
  Z123<-qqnorm(IDF$Grade123,main="Q-Q 123 Grades")$x
  hist(Z123,main=paste("123 Final Grade (Std.), N(",format(mean(IDF$Grade123,na.rm=TRUE),digits=4),",sd=",format(sd(IDF$Grade123,na.rm=TRUE),digits=4),")",sep=""))
  plot(IDF$LetterGrade123, main = paste("Chem 123 Grades; N = ",length(IDF$LetterGrade123[which(!is.na(IDF$LetterGrade123))])))
  Z159<-qqnorm(IDF$Grade159,main="Q-Q 159 Grades")$x
  hist(Z159,main=paste("159 Final Grade (Std.), N(",format(mean(IDF$Grade159,na.rm=TRUE),digits=4),",sd=",format(sd(IDF$Grade159,na.rm=TRUE),digits=4),")",sep=""))
  plot(IDF$LetterGrade159, main = paste("Chem 159 Grades; N = ",length(IDF$LetterGrade159[which(!is.na(IDF$LetterGrade159))])))
  Z227<-qqnorm(IDF$Grade227,main="Q-Q 227 Grades")$x
  hist(Z227,main=paste("227 Final Grade (Std.), N(",format(mean(IDF$Grade227,na.rm=TRUE),digits=4),",sd=",format(sd(IDF$Grade227,na.rm=TRUE),digits=4),")",sep=""))
  plot(IDF$LetterGrade227, main = paste("Chem 227 Grades; N = ",length(IDF$LetterGrade227[which(!is.na(IDF$LetterGrade227))])))
  COLS<-paste("C",1:nrow(Test),sep="")
  IDF$Score<-rowMeans(!is.na(IDF[,SCOLS]))
}

Anonymize <- function() {
  IDF<-IDF[,c(1,4:ncol(IDF))]
  IDF$ID<-c(1:nrow(IDF))
  write.csv(IDF, file = paste(OutDir,"AnonItemData.csv",sep="/"),row.names=FALSE)
}

PlotGrades<- function() {
  par(mfrow=c(2,2))
  plot(IDF$LetterGrade115, main = paste("Physics 115 Grades; N = ",length(IDF$LetterGrade115[which(!is.na(IDF$LetterGrade115))])))
  plot(IDF$LetterGrade123, main = paste("Physics 123 Grades; N = ",length(IDF$LetterGrade123[which(!is.na(IDF$LetterGrade123))])))
  plot(IDF$LetterGrade159, main = paste("Chem 159 Grades; N = ",length(IDF$LetterGrade159[which(!is.na(IDF$LetterGrade159))])))
  plot(IDF$LetterGrade227, main = paste("Physics 227 Grades; N = ",length(IDF$LetterGrade227[which(!is.na(IDF$LetterGrade227))])))
}

CutCategories<- function() {
  print("Columns worth cutting on:")
  print(colnames(IDF)[-c(which(colnames(IDF)=="C1"):ncol(IDF))])
  #which(class(IDF[,colnames(IDF)])=="factor")
  print("Factor categories:")
  print(sapply(IDF[,colnames(IDF[,which(sapply(IDF, function(x) is.factor(x)))])],levels))
  #paste(names(LVector),LVector[1:length(LVector)])
  #which(sapply(IDF, function(x) is.numeric(x)))
  #print(str(IDF[,-c(which(colnames(IDF)=="C1"):ncol(IDF))]))
}

CLASSCategories<-function() {
  print("        PHYSICS CLASS CATEGORY     :     PRETEST     :     POSTTEST    ")
  print("==============================================================================")
  print("                  Overall - Fav    :   PCLASSOverF   :   pclassOverF")
  print("                  Overall - Unfav  :   PCLASSOverU   :   pclassOverU")
  print("                      All - Fav    :   PCLASSAllF    :   pclassAllF")
  print("                      All - Unfav  :   PCLASSAllU    :   pclassAllU")
  print("        Personal Interest - Fav    :   PCLASSPIF     :   pclassPIF")
  print("        Personal Interest - Unfav  :   PCLASSPIU     :   pclassPIU")
  print("    Real World Connection - Fav    :   PCLASSRWF     :   pclassRWF")
  print("    Real World Connection - Unfav  :   PCLASSRWU     :   pclassRWU")
  print("                  General - Fav    :   PCLASSGenF    :   pclassGenF")
  print("                  General - Unfav  :   PCLASSGenU    :   pclassGenU")
  print("               Confidence - Fav    :   PCLASSConF    :   pclassConF")
  print("               Confidence - Unfav  :   PCLASSConU    :   pclassConU")
  print("           Sophistication - Fav    :   PCLASSSophF   :   pclassSophF")
  print("           Sophistication - Unfav  :   PCLASSSophU   :   pclassSophU")
  print("      Sense Making/Effort - Fav    :   PCLASSSEF     :   pclassSEF")
  print("      Sense Making/Effort - Unfav  :   PCLASSSEU     :   pclassSEU")
  print(" Conceptual Understanding - Fav    :   PCLASSCUF     :   pclassCUF")
  print(" Conceptual Understanding - Unfav  :   PCLASSCUU     :   pclassCUU")
  print("    Applied Understanding - Fav    :   PCLASSACF     :   pclassACF")
  print("    Applied Understanding - Unfav  :   PCLASSACU     :   pclassACU")
  print(" ..............................................................................")
  print("                 FOR CHEMISTRY, CHANGE THE 'P' to a 'C'   ")
  print(" ..............................................................................")
  print("")
  print("      CHEMISTRY CLASS CATEGORY     :     PRETEST     :     POSTTEST    ")
  print("==============================================================================")
  print("                  Overall - Fav    :   CCLASSOverF   :   cclassOverF")
  print("                  Overall - Unfav  :   CCLASSOverU   :   cclassOverU")
  print("                      All - Fav    :   CCLASSAllF    :   cclassAllF")
  print("                      All - Unfav  :   CCLASSAllU    :   cclassAllU")
  print("        Personal Interest - Fav    :   CCLASSPIF     :   cclassPIF")
  print("        Personal Interest - Unfav  :   CCLASSPIU     :   cclassPIU")
  print("    Real World Connection - Fav    :   CCLASSRWF     :   cclassRWF")
  print("    Real World Connection - Unfav  :   CCLASSRWU     :   cclassRWU")
  print("                  General - Fav    :   CCLASSGenF    :   cclassGenF")
  print("                  General - Unfav  :   CCLASSGenU    :   cclassGenU")
  print("               Confidence - Fav    :   CCLASSConF    :   cclassConF")
  print("               Confidence - Unfav  :   CCLASSConU    :   cclassConU")
  print("           Sophistication - Fav    :   CCLASSSophF   :   cclassSophF")
  print("           Sophistication - Unfav  :   CCLASSSophU   :   cclassSophU")
  print("      Sense Making/Effort - Fav    :   CCLASSSEF     :   cclassSEF")
  print("      Sense Making/Effort - Unfav  :   CCLASSSEU     :   cclassSEU")
  print("   Conceptual Connections - Fav    :   CCLASSCCF     :   cclassCCF")
  print("   Conceptual Connections - Unfav  :   CCLASSCCU     :   cclassCCU")
  print("      Conceptual Learning - Fav    :   CCLASSCLF     :   cclassCLF")
  print("      Conceptual Learning - Unfav  :   CCLASSCLU     :   cclassCLU")
  print("Atom/Molecule Perspective - Fav    :   CCLASSPerF    :   cclassPerF")
  print("Atom/Molecule Perspective - Unfav  :   CCLASSPerU    :   cclassPerU")
  print("==============================================================================")
  print("    ...So to access General-Fav, use IDF$PCLASSGenF ")
}