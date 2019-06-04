###########################################################################
# This file must be sourced FIRST for analysis on Test Data to RUN.
###########################################################################

#####################################################
# Some Constant Variables about the Test for Setup  #
#####################################################
InDir = "Input"
OutDir = "Output"
PhyFormMod = 8
EMFormMod = 6
ChemFormMod = 6
PrePhyTestLength = 43
PreChemTestLength = 37
EMTestLength = 37
StatTestLength = 30
PostPhyTestLength = 37
PostChemTestLength = 8

####################################################################################################
# FILL THE ITEM DATA FRAMES # Item columns C%d:>1>0>1>0>0<1<1<<<<B%d:< <4>5>3>2>5>1<1<<<<<<
####################################################################################################

IDF<-read.csv("WIDER_DATA/IDF.csv",header = TRUE,stringsAsFactors = FALSE)
IDF$ETHNIC_CD<-as.factor(IDF$ETHNIC_CD)
IDF$GENDER_CD<-as.factor(IDF$GENDER_CD)
IDF$School<-as.factor(IDF$School)
IDF$Year<-as.factor(IDF$Year)
IDF$LetterGrade115<-factor(IDF$LetterGrade115,levels=c("A","B+","B","C+","C","D","F","W","TZ"))
IDF$LetterGrade116<-factor(IDF$LetterGrade116,levels=c("A","B+","B","C+","C","D","F","W","TZ"))
IDF$LetterGrade116.S13<-factor(IDF$LetterGrade116.S13,levels=c("A","B+","B","C+","C","D","F","W","TZ"))
IDF$LetterGrade123<-factor(IDF$LetterGrade123,levels=c("A","B+","B","C+","C","D","F","W","TZ"))
IDF$LetterGrade124<-factor(IDF$LetterGrade124,levels=c("A","B+","B","C+","C","D","F","W","TZ"))
IDF$LetterGrade159<-factor(IDF$LetterGrade159,levels=c("A","B+","B","C+","C","D","F","W","TZ"))
IDF$LetterGrade227<-factor(IDF$LetterGrade227,levels=c("A","B+","B","C+","C","D","F","W","TZ"))


#####################################################################################################
# FILL THE TEST BANKS with Items, Answers, Descriptions, Etc...
#####################################################################################################
### The Pretest
TestX<-read.csv("WIDER_DATA/MasterPreLayout.csv",header=T,stringsAsFactors=F,as.is=T,sep=",")
### The Posttest
TestY<-read.csv("WIDER_DATA/MasterPostLayout.csv",header=T,stringsAsFactors=F,as.is=T,sep=",")
TestZ<-read.csv("WIDER_DATA/MasterPostPostLayout.csv",header=T,stringsAsFactors=F,as.is=T,sep=",")
### Frame for changing numbers to letters
KeySwitch<-data.frame(LLetters=c("a","b","c","d","e"),CLetters=c("A","B","C","D","E"),Numbers=c(1:5))
### The Chemistry Pretest  ########
CHPreTest<-read.csv(paste("WIDER_DATA/",InDir,"PreChemKey.csv",sep="/"),header=F,stringsAsFactors=F,as.is=T,sep=",")
names(CHPreTest)<-c("CHQs","CHAs")
CHPreTest<-merge(CHPreTest,KeySwitch,by.x="CHAs",by.y="LLetters", all=TRUE)
CHPreTest<-arrange(CHPreTest, desc(-CHQs))
CHPreTest$CH1<-NA
CHPreTest$CH2<-NA
CHPreTest$CH3<-NA
CHPreTest$CH4<-NA
CHPreTest$CH5<-NA
CHPreTest$CH6<-NA
CHPreTest$CH1<-CHPreTest$Numbers
CHPreTest$CH2<-CHPreTest$CH1
CHPreTest$CH3<-CHPreTest$CH1
CHPreTest$CH4<-CHPreTest$CH1
CHPreTest$CH5<-CHPreTest$CH1
CHPreTest$CH6<-CHPreTest$CH1
### The Chemistry Posttest  ############
CHPostTest<-read.csv(paste("WIDER_DATA/",InDir,"PostChemKey.csv",sep="/"),header=F,stringsAsFactors=F,as.is=T,sep=",")
names(CHPostTest)<-c("CHQs","CHAs")
CHPostTest<-merge(CHPostTest,KeySwitch,by.x="CHAs",by.y="LLetters", all=TRUE)
CHPostTest<-arrange(CHPostTest, desc(-CHQs))
CHPostTest$CH1<-NA
CHPostTest$CH2<-NA
CHPostTest$CH3<-NA
CHPostTest$CH4<-NA
CHPostTest$CH5<-NA
CHPostTest$CH6<-NA
CHPostTest$CH1<-CHPostTest$Numbers
CHPostTest$CH2<-CHPostTest$CH1
CHPostTest$CH3<-CHPostTest$CH1
CHPostTest$CH4<-CHPostTest$CH1
CHPostTest$CH5<-CHPostTest$CH1
CHPostTest$CH6<-CHPostTest$CH1
#### The E&M Pretest  ##########
EMTest<-read.csv(paste("WIDER_DATA/",InDir,"PreEMKey.csv",sep="/"),header=T,stringsAsFactors=F,as.is=T,sep=",")
start = length(EMTest)+1
for (i in 1:6) {
  EMTest<-merge(EMTest,KeySwitch,by.x=paste("EMv",i,"LL",sep=""),by.y="LLetters", all.x=TRUE)
  names(EMTest)[start+i]<-paste("EM",i,sep="")
  drops <- c("CLetters")
  EMTest<-EMTest[,!(names(EMTest) %in% drops)]
}
EMTest<-arrange(EMTest, desc(-EMQs))
### The FCI  #############
FCTest<-read.csv(paste("WIDER_DATA/",InDir,"PreFCIKey.csv",sep="/"),header=F,stringsAsFactors=F,as.is=T,sep=",")
names(FCTest)<-c("FCQs","FCAs")
FCTest<-merge(FCTest,KeySwitch,by.x="FCAs",by.y="CLetters", all=TRUE)
FCTest<-arrange(FCTest, desc(-FCQs))
FCTest$FC1<-NA
FCTest$FC1<-FCTest$Numbers
### The Mathematization questions  #############
MATest<-read.csv(paste("WIDER_DATA/",InDir,"PreMathKey.csv",sep="/"),header=T,stringsAsFactors=F,as.is=T,sep=",")
MATest<-merge(MATest,KeySwitch,by.x="MAAs",by.y="LLetters", all.x=TRUE)
MATest<-arrange(MATest, desc(-MAQs))
#### The Mechanics Baseline Questions  ##########
MBTest<-read.csv(paste("WIDER_DATA/",InDir,"PreMechKey.csv",sep="/"),header=T,stringsAsFactors=F,as.is=T,sep=",")
MBTest<-merge(MBTest,KeySwitch,by.x="MBAs",by.y="LLetters", all.x=TRUE)
MBTest<-arrange(MBTest, desc(-MBQs))
### The 115 and 123 Pretest  ###########
MEPreTest<-data.frame(MEQs=c(1:PrePhyTestLength),ME1=c(MATest$Numbers,FCTest$Numbers,MBTest$Numbers))
MEPreTest$ME2<-NA
MEPreTest$ME3<-NA
MEPreTest$ME4<-NA
MEPreTest$ME5<-NA
MEPreTest$ME6<-NA
MEPreTest$ME7<-NA
MEPreTest$ME8<-NA
MEPreTest$ME2<-MEPreTest$ME1
MEPreTest$ME3<-MEPreTest$ME1
MEPreTest$ME4<-MEPreTest$ME1
MEPreTest$ME5<-MEPreTest$ME1
MEPreTest$ME6<-MEPreTest$ME1
MEPreTest$ME7<-MEPreTest$ME1
MEPreTest$ME8<-MEPreTest$ME1
MEPreTest$ME7[1]<-2
MEPreTest$ME8[1]<-2
### The 115 and 123 Posttest  ##############
MEPostTest<-data.frame(MEQs=c(1:PostPhyTestLength),ME1=c(MATest$Numbers,FCTest$Numbers))
MEPostTest$ME2<-NA
MEPostTest$ME3<-NA
MEPostTest$ME4<-NA
MEPostTest$ME5<-NA
MEPostTest$ME6<-NA
MEPostTest$ME7<-NA
MEPostTest$ME8<-NA
MEPostTest$ME2<-MEPostTest$ME1
MEPostTest$ME3<-MEPostTest$ME1
MEPostTest$ME4<-MEPostTest$ME1
MEPostTest$ME5<-MEPostTest$ME1
MEPostTest$ME6<-MEPostTest$ME1
MEPostTest$ME7<-MEPostTest$ME1
MEPostTest$ME8<-MEPostTest$ME1
MEPostTest$ME7[1]<-2
MEPostTest$ME8[1]<-2


#####################################################
### CREATE FCI Misconceptions Taxonomy
### For Analysis of FCI
######################################################
# AF 1-7, AR 1-2, CI 1-3, CF, Ob, R 1-3, G 1-5
FCI_MCstr<-c(paste("K",1:4,sep=""),
             paste("I",1:5,sep=""),
             paste("AF",1:7,sep=""),
             "AR1","AR2","CI1","CI2","CI3","CF","Ob","R1","R2","R3",
             paste("G",1:5,sep=""))
length(FCI_MCstr)
MC3D<-rep(0,length(FCI_MCstr)*length(FCTest$FC1)*6)
MC3D<-array(MC3D,c(length(FCTest$FC1),6,length(FCI_MCstr)))
MC3D[19,2:4,which(FCI_MCstr=="K1")]<-rep(1,3)  # B,C,D
MC3D[19,1,which(FCI_MCstr=="K2")]<-1   
MC3D[20,2:3,which(FCI_MCstr=="K2")]<-rep(1,2)
MC3D[9,3,which(FCI_MCstr=="K3")]<-1          
MC3D[14,1:2,which(FCI_MCstr=="K4")]<-rep(1,2) 
# 5C,D,E;11B,C;27D;30B,D,E
MC3D[5,3:5,which(FCI_MCstr=="I1")]<-rep(1,3) 
MC3D[11,2:3,which(FCI_MCstr=="I1")]<-rep(1,2) 
MC3D[27,4,which(FCI_MCstr=="I1")]<-1 
MC3D[30,c(2,4,5),which(FCI_MCstr=="I1")]<-rep(1,3) 
MC3D[7,4,which(FCI_MCstr=="I2")]<-1 
MC3D[8,c(3,5),which(FCI_MCstr=="I2")]<-rep(1,2) 
MC3D[21,1,which(FCI_MCstr=="I2")]<-1 
MC3D[23,c(1,4),which(FCI_MCstr=="I2")]<-rep(1,2) 
MC3D[12,3:4,which(FCI_MCstr=="I3")]<-rep(1,2) 
MC3D[13,1:3,which(FCI_MCstr=="I3")]<-rep(1,3) 
MC3D[14,5,which(FCI_MCstr=="I3")]<-1 
MC3D[23,4,which(FCI_MCstr=="I3")]<-1 
MC3D[24,c(3,5),which(FCI_MCstr=="I3")]<-rep(1,2) 
MC3D[27,2,which(FCI_MCstr=="I3")]<-1 
MC3D[8,4,which(FCI_MCstr=="I4")]<-1 
MC3D[10,c(2,4),which(FCI_MCstr=="I4")]<-rep(1,2) 
MC3D[21,4,which(FCI_MCstr=="I4")]<-1 
MC3D[23,5,which(FCI_MCstr=="I4")]<-1 
MC3D[26,3,which(FCI_MCstr=="I4")]<-1 
MC3D[27,5,which(FCI_MCstr=="I4")]<-1 
MC3D[5,3:5,which(FCI_MCstr=="I5")]<-rep(1,3) 
MC3D[6,1,which(FCI_MCstr=="I5")]<-1 
MC3D[7,c(1,4),which(FCI_MCstr=="I5")]<-rep(1,2) 
MC3D[18,3:4,which(FCI_MCstr=="I5")]<-rep(1,2) 
MC3D[15:16,4,which(FCI_MCstr=="AF1")]<-rep(1,2) 
MC3D[17,5,which(FCI_MCstr=="AF1")]<-1 
MC3D[18,1,which(FCI_MCstr=="AF1")]<-1 
MC3D[28,2,which(FCI_MCstr=="AF1")]<-1 
MC3D[30,1,which(FCI_MCstr=="AF1")]<-1 
MC3D[5,3:5,which(FCI_MCstr=="AF2")]<-rep(1,3) 
MC3D[27,1,which(FCI_MCstr=="AF2")]<-1 
MC3D[29,5,which(FCI_MCstr=="AF3")]<-1 
MC3D[c(22,26),1,which(FCI_MCstr=="AF4")]<-rep(1,2) 
MC3D[3,2,which(FCI_MCstr=="AF5")]<-1 
MC3D[3,1,which(FCI_MCstr=="AF6")]<-1 
MC3D[c(22,26),4,which(FCI_MCstr=="AF6")]<-rep(1,2) 
MC3D[22,c(3,5),which(FCI_MCstr=="AF7")]<-rep(1,2) 
MC3D[4,c(1,4),which(FCI_MCstr=="AR1")]<-rep(1,2) 
MC3D[c(15,16),2,which(FCI_MCstr=="AR1")]<-rep(1,2) 
MC3D[28,4,which(FCI_MCstr=="AR1")]<-1 
MC3D[c(15,16),3,which(FCI_MCstr=="AR2")]<-rep(1,2) 
MC3D[28,4,which(FCI_MCstr=="AR2")]<-1
MC3D[17,c(1,4),which(FCI_MCstr=="CI1")]<-rep(1,2)
MC3D[25,5,which(FCI_MCstr=="CI1")]<-1
MC3D[6,4,which(FCI_MCstr=="CI2")]<-1
MC3D[c(7,14,21),3,which(FCI_MCstr=="CI2")]<-rep(1,3)
MC3D[12,1,which(FCI_MCstr=="CI2")]<-1
MC3D[8,1,which(FCI_MCstr=="CI3")]<-1
MC3D[c(9,21),2,which(FCI_MCstr=="CI3")]<-rep(1,2)
MC3D[23,3,which(FCI_MCstr=="CI3")]<-1
MC3D[c(5,18),5,which(FCI_MCstr=="CF")]<-rep(1,2)
MC3D[6:7,3:5,which(FCI_MCstr=="CF")]<-rep(1,6)
MC3D[4,3,which(FCI_MCstr=="Ob")]<-1
MC3D[c(5,18,29),1,which(FCI_MCstr=="Ob")]<-rep(1,3)
MC3D[11,1:2,which(FCI_MCstr=="Ob")]<-rep(1,2)
MC3D[15:16,5,which(FCI_MCstr=="Ob")]<-rep(1,2)
MC3D[27,1:2,which(FCI_MCstr=="R1")]<-rep(1,2)
MC3D[25,c(1,2,4),which(FCI_MCstr=="R2")]<-rep(1,3)
MC3D[26,2,which(FCI_MCstr=="R2")]<-1
MC3D[26,2,which(FCI_MCstr=="R3")]<-1
MC3D[3,5,which(FCI_MCstr=="G1")]<-1
MC3D[11,1,which(FCI_MCstr=="G1")]<-1
MC3D[c(17,29),4,which(FCI_MCstr=="G1")]<-rep(1,2)
MC3D[29,3,which(FCI_MCstr=="G1")]<-1
MC3D[3,4,which(FCI_MCstr=="G2")]<-1
MC3D[c(11,13),5,which(FCI_MCstr=="G2")]<-rep(1,2)
MC3D[29,3,which(FCI_MCstr=="G2")]<-1
MC3D[1,1,which(FCI_MCstr=="G3")]<-1
MC3D[2,c(2,4),which(FCI_MCstr=="G3")]<-rep(1,2)
MC3D[c(3,13),2,which(FCI_MCstr=="G4")]<-rep(1,2)
MC3D[12,4,which(FCI_MCstr=="G5")]<-1
MC3D[13,2,which(FCI_MCstr=="G5")]<-1
MC3D[14,5,which(FCI_MCstr=="G5")]<-1
##  MC3D is 30x5x31 : length(FCI_MCstr)*length(FCTest$FC1)*5
FC_ItemMC<-apply(MC3D,1,sum)
FC_MCStats<-apply(MC3D,3,sum)
names(FC_MCStats)<-FCI_MCstr
FC_MCstats<-rep(0,length(FCI_MCstr))
names(FC_MCstats)<-FCI_MCstr
for (j in 1:length(FCI_MCstr)) {
  for (i in 1:length(FCTest$FC1)) {
    FC_MCstats[j]=FC_MCstats[j]+max(MC3D[i,,j])
  }
}

