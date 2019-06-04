
source("InitializeGeisCamilli.R")
basedir<-getwd()

data.dir<-"RealData"
if (!data.dir %in% dir()) dir.create(data.dir)

#source("WIDER_DATA/InitiateData.R")
# FCIC<-paste("C",64:93,sep="")
# Response<-IDF[which(!is.na(IDF$FCI.x)&(!is.na(IDF$LetterGrade123)|!is.na(IDF$LetterGrade115))),FCIC]
# write.csv(Response,paste0(data.dir,"/FCIPreTest.csv"),row.names=FALSE)

Response<-read.csv(paste0(data.dir,"/FCIPreTest.csv"),stringsAsFactors = FALSE,header = TRUE)
settings<-CheckParams(generate = FALSE)
settings$Adim<-3
settings$guess<-TRUE
settings$dbltrunc<-FALSE
settings$tmu<-rep(0,settings$Adim)
settings$tsigma<-diag(settings$Adim)
settings$rmethod<-NA
settings$empiricalse<-TRUE
settings$esttheta<-TRUE
settings$nesttheta<-200
settings$thinA=8
settings$thinB=5
settings$EmpIT=1000
settings$estfile<-paste0(data.dir,"/FCI_Guess_A",settings$Adim)
settings$burnin=2000
settings$eps=0.0001
settings$thetamap=FALSE
settings$record=TRUE
settings$parallel<-FALSE
FCI.Fit.3D.Guess<-AnalyzeTestData(RP = Response,settings = settings,verbose = TRUE,TargetA = NA)
FCI.Fit.3D.Guess
ROT<-RotateSlopes(fit.data = FCI.Fit.3D.Guess)
GetLikelihood(FCI.Fit.3D.Guess)
TWFitTest(fit.data = FCI.Fit.3D.Guess)


settings$Adim<-2
settings$tmu<-rep(0,settings$Adim)
settings$tsigma<-diag(settings$Adim)
settings$estfile<-paste0(data.dir,"/FCI_Guess_A",settings$Adim)
FCI.Fit.2D.Guess<-AnalyzeTestData(RP = Response,settings = settings,verbose = TRUE,TargetA = NA)
GetLikelihood(FCI.Fit.2D.Guess)
TWFitTest(fit.data = FCI.Fit.2D.Guess)
LRTest(FCI.Fit.3D.Guess,FCI.Fit.2D.Guess)

settings$Adim<-1
settings$tmu<-rep(0,settings$Adim)
settings$tsigma<-diag(settings$Adim)
settings$estfile<-paste0(data.dir,"/FCI_Guess_A",settings$Adim)
FCI.Fit.1D.Guess<-AnalyzeTestData(RP = Response,settings = settings,verbose = TRUE,TargetA = NA)
GetLikelihood(FCI.Fit.1D.Guess)
TWFitTest(fit.data = FCI.Fit.1D.Guess)
LRTest(FCI.Fit.2D.Guess,FCI.Fit.1D.Guess)

settings$Adim<-1
settings$guess<-FALSE;settings$tmu<-rep(0,settings$Adim);settings$tsigma<-diag(settings$Adim)
settings$estfile<-paste0(data.dir,"/FCI_NoGuess_A",settings$Adim)
FCI.Fit.1D.NoGuess<-AnalyzeTestData(RP = Response,settings = settings,verbose = TRUE,TargetA = NA)
GetLikelihood(FCI.Fit.1D.NoGuess)
TWFitTest(fit.data = FCI.Fit.1D.NoGuess)
LRTest(FCI.Fit.1D.Guess,FCI.Fit.1D.NoGuess)

settings$Adim<-2;settings$tmu<-rep(0,settings$Adim);settings$tsigma<-diag(settings$Adim)
settings$estfile<-paste0(data.dir,"/FCI_NoGuess_A",settings$Adim)
FCI.Fit.2D.NoGuess<-AnalyzeTestData(RP = Response,settings = settings,verbose = TRUE,TargetA = NA)
GetLikelihood(FCI.Fit.2D.NoGuess)
TWFitTest(fit.data = FCI.Fit.2D.NoGuess)
LRTest(FCI.Fit.2D.Guess,FCI.Fit.2D.NoGuess)
LRTest(FCI.Fit.2D.NoGuess,FCI.Fit.1D.NoGuess)

settings$Adim<-3;settings$tmu<-rep(0,settings$Adim);settings$tsigma<-diag(settings$Adim)
settings$estfile<-paste0(data.dir,"/FCI_NoGuess_A",settings$Adim)
FCI.Fit.3D.NoGuess<-AnalyzeTestData(RP = Response,settings = settings,verbose = TRUE,TargetA = NA)
GetLikelihood(FCI.Fit.3D.NoGuess)
TWFitTest(fit.data = FCI.Fit.3D.NoGuess,tw=0.95)
LRTest(FCI.Fit.3D.Guess,FCI.Fit.3D.NoGuess)
LRTest(FCI.Fit.3D.NoGuess,FCI.Fit.2D.NoGuess)

settings$Adim<-4;settings$tmu<-rep(0,settings$Adim);settings$tsigma<-diag(settings$Adim)
settings$estfile<-paste0(data.dir,"/FCI_NoGuess_A",settings$Adim)
FCI.Fit.4D.NoGuess<-AnalyzeTestData(RP = Response,settings = settings,verbose = TRUE,TargetA = NA)
GetLikelihood(FCI.Fit.4D.NoGuess)
TWFitTest(fit.data = FCI.Fit.4D.NoGuess)
LRTest(FCI.Fit.4D.NoGuess,FCI.Fit.3D.NoGuess)

settings$Adim<-5;settings$tmu<-rep(0,settings$Adim);settings$tsigma<-diag(settings$Adim)
settings$estfile<-paste0(data.dir,"/FCI_NoGuess_A",settings$Adim)
FCI.Fit.5D.NoGuess<-AnalyzeTestData(RP = Response,settings = settings,verbose = TRUE,TargetA = NA)
GetLikelihood(FCI.Fit.5D.NoGuess)
TWFitTest(fit.data = FCI.Fit.5D.NoGuess)
LRTest(FCI.Fit.5D.NoGuess,FCI.Fit.4D.NoGuess)
LRTest(FCI.Fit.5D.NoGuess,FCI.Fit.2D.NoGuess)

settings$Adim<-6; settings$tmu<-rep(0,settings$Adim); settings$tsigma<-diag(settings$Adim)
settings$estfile<-paste0(data.dir,"/FCI_NoGuess_A",settings$Adim)
FCI.Fit.6D.NoGuess<-AnalyzeTestData(RP = Response,settings = settings,verbose = TRUE,TargetA = NA)
GetLikelihood(FCI.Fit.6D.NoGuess)
TWFitTest(fit.data = FCI.Fit.6D.NoGuess)
LRTest(FCI.Fit.6D.NoGuess,FCI.Fit.5D.NoGuess)
LRTest(FCI.Fit.6D.NoGuess,FCI.Fit.3D.NoGuess)

settings$Adim<-7; settings$tmu<-rep(0,settings$Adim); settings$tsigma<-diag(settings$Adim)
settings$estfile<-paste0(data.dir,"/FCI_NoGuess_A",settings$Adim)
FCI.Fit.7D.NoGuess<-AnalyzeTestData(RP = Response,settings = settings,verbose = TRUE,TargetA = NA)
GetLikelihood(FCI.Fit.7D.NoGuess)
TWFitTest(fit.data = FCI.Fit.7D.NoGuess)

settings$Adim<-8; settings$tmu<-rep(0,settings$Adim); settings$tsigma<-diag(settings$Adim)
settings$estfile<-paste0(data.dir,"/FCI_NoGuess_A",settings$Adim)
FCI.Fit.8D.NoGuess<-AnalyzeTestData(RP = Response,settings = settings,verbose = TRUE,TargetA = NA)
GetLikelihood(FCI.Fit.8D.NoGuess)
TWFitTest(fit.data = FCI.Fit.8D.NoGuess)

settings$Adim<-9; settings$tmu<-rep(0,settings$Adim); settings$tsigma<-diag(settings$Adim)
settings$estfile<-paste0(data.dir,"/FCI_NoGuess_A",settings$Adim)
FCI.Fit.9D.NoGuess<-AnalyzeTestData(RP = Response,settings = settings,verbose = TRUE,TargetA = NA)
GetLikelihood(FCI.Fit.9D.NoGuess)
TWFitTest(fit.data = FCI.Fit.9D.NoGuess)

settings$Adim<-10; settings$tmu<-rep(0,settings$Adim); settings$tsigma<-diag(settings$Adim)
settings$estfile<-paste0(data.dir,"/FCI_NoGuess_A",settings$Adim)
FCI.Fit.10D.NoGuess<-AnalyzeTestData(RP = Response,settings = settings,verbose = TRUE,TargetA = NA)
GetLikelihood(FCI.Fit.10D.NoGuess)
TWFitTest(fit.data = FCI.Fit.10D.NoGuess)

settings$Adim<-11; settings$tmu<-rep(0,settings$Adim); settings$tsigma<-diag(settings$Adim)
settings$estfile<-paste0(data.dir,"/FCI_NoGuess_A",settings$Adim)
FCI.Fit.11D.NoGuess<-AnalyzeTestData(RP = Response,settings = settings,verbose = TRUE,TargetA = NA)
GetLikelihood(FCI.Fit.11D.NoGuess)
TWFitTest(fit.data = FCI.Fit.11D.NoGuess)

settings$Adim<-12; settings$tmu<-rep(0,settings$Adim); settings$tsigma<-diag(settings$Adim)
settings$estfile<-paste0(data.dir,"/FCI_NoGuess_A",settings$Adim)
FCI.Fit.12D.NoGuess<-AnalyzeTestData(RP = Response,settings = settings,verbose = TRUE,TargetA = NA)
GetLikelihood(FCI.Fit.12D.NoGuess)
TWFitTest(fit.data = FCI.Fit.12D.NoGuess)



########### CCI

#source("WIDER_DATA/InitiateData.R")
# CCIC<-paste("C",94:115,sep="")
# CCI.Response<-IDF[!is.na(IDF$CCI.x)&!is.na(IDF$LetterGrade159),CCIC]
# write.csv(CCI.Response,paste0(data.dir,"/CCIPreTest.csv"),row.names=FALSE)

Response<-read.csv(paste0(data.dir,"/CCIPreTest.csv"),stringsAsFactors = FALSE,header = TRUE)
settings<-CheckParams(generate = FALSE)
settings$Adim<-3
settings$guess<-TRUE
settings$dbltrunc<-FALSE
settings$tmu<-rep(0,settings$Adim)
settings$tsigma<-diag(settings$Adim)
settings$rmethod<-NA
settings$empiricalse<-TRUE
settings$esttheta<-TRUE
settings$nesttheta<-200
settings$thinA=8
settings$thinB=5
settings$EmpIT=1000
settings$estfile<-paste0(data.dir,"/CCI_Guess_A",settings$Adim)
settings$burnin=2000
settings$eps=0.0001
settings$thetamap=FALSE
settings$record=TRUE
settings$parallel<-FALSE
CCI.Fit.3D.Guess<-AnalyzeTestData(RP = Response,settings = settings,verbose = TRUE,TargetA = NA)
CCI.Fit.3D.Guess
ROT<-RotateSlopes(fit.data = CCI.Fit.3D.Guess)
GetLikelihood(CCI.Fit.3D.Guess)
TWFitTest(fit.data = CCI.Fit.3D.Guess)


settings$Adim<-2
settings$tmu<-rep(0,settings$Adim)
settings$tsigma<-diag(settings$Adim)
settings$estfile<-paste0(data.dir,"/CCI_Guess_A",settings$Adim)
CCI.Fit.2D.Guess<-AnalyzeTestData(RP = Response,settings = settings,verbose = TRUE,TargetA = NA)
GetLikelihood(CCI.Fit.2D.Guess)
TWFitTest(fit.data = CCI.Fit.2D.Guess)
LRTest(CCI.Fit.3D.Guess,CCI.Fit.2D.Guess)

settings$Adim<-1
settings$tmu<-rep(0,settings$Adim)
settings$tsigma<-diag(settings$Adim)
settings$estfile<-paste0(data.dir,"/CCI_Guess_A",settings$Adim)
CCI.Fit.1D.Guess<-AnalyzeTestData(RP = Response,settings = settings,verbose = TRUE,TargetA = NA)
GetLikelihood(CCI.Fit.1D.Guess)
TWFitTest(fit.data = CCI.Fit.1D.Guess)
LRTest(CCI.Fit.2D.Guess,CCI.Fit.1D.Guess)

settings$Adim<-1
settings$guess<-FALSE;settings$tmu<-rep(0,settings$Adim);settings$tsigma<-diag(settings$Adim)
settings$estfile<-paste0(data.dir,"/CCI_NoGuess_A",settings$Adim)
CCI.Fit.1D.NoGuess<-AnalyzeTestData(RP = Response,settings = settings,verbose = TRUE,TargetA = NA)
GetLikelihood(CCI.Fit.1D.NoGuess)
TWFitTest(fit.data = CCI.Fit.1D.NoGuess)
LRTest(CCI.Fit.1D.Guess,CCI.Fit.1D.NoGuess)

settings$Adim<-2;settings$tmu<-rep(0,settings$Adim);settings$tsigma<-diag(settings$Adim)
settings$estfile<-paste0(data.dir,"/CCI_NoGuess_A",settings$Adim)
CCI.Fit.2D.NoGuess<-AnalyzeTestData(RP = Response,settings = settings,verbose = TRUE,TargetA = NA)
GetLikelihood(CCI.Fit.2D.NoGuess)
TWFitTest(fit.data = CCI.Fit.2D.NoGuess)
LRTest(CCI.Fit.2D.Guess,CCI.Fit.2D.NoGuess)
LRTest(CCI.Fit.2D.NoGuess,CCI.Fit.1D.NoGuess)

settings$Adim<-3;settings$tmu<-rep(0,settings$Adim);settings$tsigma<-diag(settings$Adim)
settings$estfile<-paste0(data.dir,"/CCI_NoGuess_A",settings$Adim)
CCI.Fit.3D.NoGuess<-AnalyzeTestData(RP = Response,settings = settings,verbose = TRUE,TargetA = NA)
GetLikelihood(CCI.Fit.3D.NoGuess)
TWFitTest(fit.data = CCI.Fit.3D.NoGuess,tw=0.95)
LRTest(CCI.Fit.3D.Guess,CCI.Fit.3D.NoGuess)
LRTest(CCI.Fit.3D.NoGuess,CCI.Fit.2D.NoGuess)

settings$Adim<-4;settings$tmu<-rep(0,settings$Adim);settings$tsigma<-diag(settings$Adim)
settings$estfile<-paste0(data.dir,"/CCI_NoGuess_A",settings$Adim)
CCI.Fit.4D.NoGuess<-AnalyzeTestData(RP = Response,settings = settings,verbose = TRUE,TargetA = NA)
GetLikelihood(CCI.Fit.4D.NoGuess)
TWFitTest(fit.data = CCI.Fit.4D.NoGuess)
LRTest(CCI.Fit.4D.NoGuess,CCI.Fit.3D.NoGuess)

settings$Adim<-5;settings$tmu<-rep(0,settings$Adim);settings$tsigma<-diag(settings$Adim)
settings$estfile<-paste0(data.dir,"/CCI_NoGuess_A",settings$Adim)
CCI.Fit.5D.NoGuess<-AnalyzeTestData(RP = Response,settings = settings,verbose = TRUE,TargetA = NA)
GetLikelihood(CCI.Fit.5D.NoGuess)
TWFitTest(fit.data = CCI.Fit.5D.NoGuess)
LRTest(CCI.Fit.5D.NoGuess,CCI.Fit.4D.NoGuess)
LRTest(CCI.Fit.5D.NoGuess,CCI.Fit.2D.NoGuess)

settings$Adim<-6; settings$tmu<-rep(0,settings$Adim); settings$tsigma<-diag(settings$Adim)
settings$estfile<-paste0(data.dir,"/CCI_NoGuess_A",settings$Adim)
CCI.Fit.6D.NoGuess<-AnalyzeTestData(RP = Response,settings = settings,verbose = TRUE,TargetA = NA)
GetLikelihood(CCI.Fit.6D.NoGuess)
TWFitTest(fit.data = CCI.Fit.6D.NoGuess)
LRTest(CCI.Fit.6D.NoGuess,CCI.Fit.5D.NoGuess)
LRTest(CCI.Fit.6D.NoGuess,CCI.Fit.3D.NoGuess)

settings$Adim<-7; settings$tmu<-rep(0,settings$Adim); settings$tsigma<-diag(settings$Adim)
settings$estfile<-paste0(data.dir,"/CCI_NoGuess_A",settings$Adim)
CCI.Fit.7D.NoGuess<-AnalyzeTestData(RP = Response,settings = settings,verbose = TRUE,TargetA = NA)
GetLikelihood(CCI.Fit.7D.NoGuess)
TWFitTest(fit.data = CCI.Fit.7D.NoGuess)
LRTest(CCI.Fit.7D.NoGuess,CCI.Fit.6D.NoGuess)
LRTest(CCI.Fit.7D.NoGuess,CCI.Fit.4D.NoGuess)


settings$Adim<-10; settings$tmu<-rep(0,settings$Adim); settings$tsigma<-diag(settings$Adim)
settings$estfile<-paste0(data.dir,"/CCI_NoGuess_A",settings$Adim)
CCI.Fit.10D.NoGuess<-AnalyzeTestData(RP = Response,settings = settings,verbose = TRUE,TargetA = NA)
GetLikelihood(CCI.Fit.10D.NoGuess)
TWFitTest(fit.data = CCI.Fit.10D.NoGuess)
LRTest(CCI.Fit.10D.NoGuess,CCI.Fit.6D.NoGuess)
LRTest(CCI.Fit.10D.NoGuess,CCI.Fit.4D.NoGuess)

### QOL Data

Response = scan(paste0(data.dir,"/qol.dat"), what = "numeric",sep = "\n")
Response = data.matrix(do.call(rbind,lapply(strsplit(Response,"\\s"), function(x) as.numeric(x[-1]))))
settings<-CheckParams(generate = FALSE)
settings$ncat<-5
settings$Adim<-1
settings$dbltrunc<-TRUE
settings$tmu<-rep(0,settings$Adim)
settings$tsigma<-diag(settings$Adim)
settings$rmethod<-NA
settings$empiricalse<-TRUE
settings$esttheta<-TRUE
settings$nesttheta<-200
settings$thinA=8
settings$thinB=5
settings$EmpIT=1000
settings$estfile<-paste0(data.dir,"/QOL_A",settings$Adim)
settings$burnin=2000
settings$eps=0.0001
settings$thetamap=FALSE
settings$record=TRUE
settings$cores<-4
settings$parallel<-TRUE
QOL.Fit.1D<-AnalyzeTestData(RP = Response,settings = settings,verbose = TRUE,TargetA = NA)
GetLikelihood(QOL.Fit.1D)
TWFitTest(fit.data = QOL.Fit.1D)

settings$Adim<-2; settings$tmu<-rep(0,settings$Adim); settings$tsigma<-diag(settings$Adim)
settings$estfile<-paste0(data.dir,"/QOL_A",settings$Adim)
QOL.Fit.2D<-AnalyzeTestData(RP = Response,settings = settings,verbose = TRUE,TargetA = NA)
GetLikelihood(QOL.Fit.2D)
TWFitTest(fit.data = QOL.Fit.2D)

settings$Adim<-3; settings$tmu<-rep(0,settings$Adim); settings$tsigma<-diag(settings$Adim)
settings$estfile<-paste0(data.dir,"/QOL_A",settings$Adim)
QOL.Fit.3D<-AnalyzeTestData(RP = Response,settings = settings,verbose = TRUE,TargetA = NA)
GetLikelihood(QOL.Fit.3D)
TWFitTest(fit.data = QOL.Fit.3D)

settings$Adim<-4; settings$tmu<-rep(0,settings$Adim); settings$tsigma<-diag(settings$Adim)
settings$estfile<-paste0(data.dir,"/QOL_A",settings$Adim)
QOL.Fit.4D<-AnalyzeTestData(RP = Response,settings = settings,verbose = TRUE,TargetA = NA)
GetLikelihood(QOL.Fit.4D)
TWFitTest(fit.data = QOL.Fit.4D)

settings$Adim<-5; settings$tmu<-rep(0,settings$Adim); settings$tsigma<-diag(settings$Adim)
settings$estfile<-paste0(data.dir,"/QOL_A",settings$Adim)
QOL.Fit.5D<-AnalyzeTestData(RP = Response,settings = settings,verbose = TRUE,TargetA = NA)
GetLikelihood(QOL.Fit.5D)
TWFitTest(fit.data = QOL.Fit.5D)

settings$Adim<-6; settings$tmu<-rep(0,settings$Adim); settings$tsigma<-diag(settings$Adim)
settings$estfile<-paste0(data.dir,"/QOL_A",settings$Adim)
QOL.Fit.6D<-AnalyzeTestData(RP = Response,settings = settings,verbose = TRUE,TargetA = NA)
GetLikelihood(QOL.Fit.6D)
TWFitTest(fit.data = QOL.Fit.6D)

settings$Adim<-7; settings$tmu<-rep(0,settings$Adim); settings$tsigma<-diag(settings$Adim)
settings$estfile<-paste0(data.dir,"/QOL_A",settings$Adim)
QOL.Fit.7D<-AnalyzeTestData(RP = Response,settings = settings,verbose = TRUE,TargetA = NA)
GetLikelihood(QOL.Fit.7D)
TWFitTest(fit.data = QOL.Fit.7D)

settings$Adim<-10; settings$tmu<-rep(0,settings$Adim); settings$tsigma<-diag(settings$Adim)
settings$estfile<-paste0(data.dir,"/QOL_A",settings$Adim)
QOL.Fit.10D<-AnalyzeTestData(RP = Response,settings = settings,verbose = TRUE,TargetA = NA)
GetLikelihood(QOL.Fit.10D)
TWFitTest(fit.data = QOL.Fit.10D)

MultipleTWFitTests(fit.data.list = list(QOL.Fit.2D,
                                        QOL.Fit.3D,
                                        QOL.Fit.4D,
                                        QOL.Fit.5D,
                                        QOL.Fit.10D),title = "QOL")

source("WIDER_DATA/ItemFunctions.R")
source("WIDER_DATA/Demographics.R")
source("WIDER_DATA/LatentClass.R")
source("WIDER_DATA/Predictions.R")
source("WIDER_DATA/IRT.R")
source("WIDER_DATA/Diagnostics.R")
source("WIDER_DATA/Options.R")
source("WIDER_DATA/CLASS.R")
source("WIDER_DATA/GFA.R")

FCIC<-paste("C",64:93,sep="")
Response<-IDF[which(!is.na(IDF$FCI.x)&(!is.na(IDF$LetterGrade123)|!is.na(IDF$LetterGrade115))),FCIC]

IRTFit<-ApplyIRT(RP=Response,model="ogive")
GFAFit<-ApplyGFA(RP=Response,GParam=FALSE,AParams=1,HDQs=NA,Plot=FALSE)
points(1:60,c(IRTFit$A,IRTFit$B),col=6,cex=2)



CCols<-paste("C",1:129,sep="")
QPCols<-paste("QP",1:42,sep="")
QCCols<-paste("QC",1:50,sep="")
Cols<-c(CCols,QPCols,QCCols)

plot(density(IDF$CSEM.y[which(!is.na(IDF$Grade116))],na.rm=TRUE))
# lines(density(IDF$CSEM.y[which(is.na(IDF$Grade116))],na.rm=TRUE),col=2)
# plot(density(IDF$CSEM.x[which(!is.na(IDF$Grade116))],na.rm=TRUE))
# lines(density(IDF$CSEM.x[which(is.na(IDF$Grade116))],na.rm=TRUE),col=2)