
source("InitializeGeisCamilli.R")
basedir<-getwd()

data.dir<-"RealData"
if (!data.dir %in% dir()) dir.create(data.dir)

#source("WIDER_DATA/InitiateData.R")
# FCIC<-paste("C",64:93,sep="")
# Response<-IDF[which(!is.na(IDF$FCI.x)&(!is.na(IDF$LetterGrade123)|!is.na(IDF$LetterGrade115))),FCIC]
# write.csv(Response,paste0(data.dir,"/FCIPreTest.csv"),row.names=FALSE)

Response<-read.csv(paste0(data.dir,"/FCIPreTest.csv"),stringsAsFactors = FALSE,header = TRUE) # N=748
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

settings$Adim<-4;settings$tmu<-rep(0,settings$Adim);settings$tsigma<-diag(settings$Adim)
settings$estfile<-paste0(data.dir,"/FCI_Guess_A",settings$Adim)
FCI.Fit.4D.Guess<-AnalyzeTestData(RP = Response,settings = settings,verbose = TRUE,TargetA = NA)
TWFitTest(fit.data = FCI.Fit.4D.Guess)
LRTest(FCI.Fit.4D.Guess,FCI.Fit.3D.Guess)


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
# GetLikelihood(FCI.Fit.6D.NoGuess)
# TWFitTest(fit.data = FCI.Fit.6D.NoGuess)
LRTest(FCI.Fit.6D.NoGuess,FCI.Fit.5D.NoGuess)
LRTest(FCI.Fit.6D.NoGuess,FCI.Fit.3D.NoGuess)

settings$Adim<-7; settings$tmu<-rep(0,settings$Adim); settings$tsigma<-diag(settings$Adim)
settings$estfile<-paste0(data.dir,"/FCI_NoGuess_A",settings$Adim)
FCI.Fit.7D.NoGuess<-AnalyzeTestData(RP = Response,settings = settings,verbose = TRUE,TargetA = NA)
#GetLikelihood(FCI.Fit.7D.NoGuess)
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

###############################################

source("WIDER_DATA/InitiateData.R")
source("WIDER_DATA/ItemFunctions.R")

WIDER.data<-read.csv("WIDER_DATA/IDF.csv",header=TRUE)
ItemStats(WIDER.data,it=64)

FCIC<-paste("C",64:93,sep="")
Response<-IDF[which(!is.na(IDF$FCI.x)&(!is.na(IDF$LetterGrade123)|!is.na(IDF$LetterGrade115))),FCIC]
FCIB<-paste("B",64:93,sep="")
b.RP<-IDF[which(!is.na(IDF$FCI.x)&(!is.na(IDF$LetterGrade123)|!is.na(IDF$LetterGrade115))),FCIB]
Correct<-FCTest[,"Numbers"]
scores<-table(100*round(rowMeans(Response),digits = 3))
par(mfrow=c(1,2))
barplot(scores,main="FCI Pretest Score Distribution",names.arg = names(scores),xlab = "% Correct",ylab = "Frequency",las=2)
qqnorm(y=100*round(rowMeans(Response),digits = 3),main="FCI Pretest Q-Q Analysis")




fitlist.score<-ResponseCurves(responses = b.RP,scores = rowMeans(Response),
                              prows = 1,pcols = 3,correct = Correct,j.legend = c(1,10))
fitlist.1D<-ResponseCurves(responses = b.RP,scores = FCI.1D.NoGuess$That[,"Theta"],prows = 3,pcols = 5,correct = Correct)
fitlist.1D.G<-ResponseCurves(responses = b.RP,scores = FCI.1D.Guess$That[,"Theta"],prows = 3,pcols = 5,correct = Correct)
par(mfrow=c(1,4)) 
resid = data.frame(Mean.Score=numeric(),SD.Score=numeric(),Mean.Theta=numeric(),SD.Theta=numeric(),Mean.Theta.G=numeric(),SD.Theta.G=numeric())
for (j in 1:ncol(Response)) {
  plot(density(fitlist.score[[j]]$y-fitlist.score[[j]]$fitted.values),main="Residual Comparisons",
       xlab=expression("P(Y=1) -"~hat(P)(Y=1)),lty=2,ylab="density",
       ylim=range(c(density(fitlist.1D[[j]]$y-fitlist.1D[[j]]$fitted.values)$y,density(fitlist.score[[j]]$y-fitlist.score[[j]]$fitted.values)$y)))
    lines(density(fitlist.1D[[j]]$y-fitlist.1D[[j]]$fitted.values),lty=1,col=2)
    resid[j,]<-c(mean(fitlist.score[[j]]$y-fitlist.score[[j]]$fitted.values),
                 sd(fitlist.score[[j]]$y-fitlist.score[[j]]$fitted.values),
                 mean(fitlist.1D[[j]]$y-fitlist.1D[[j]]$fitted.values),
                 sd(fitlist.1D[[j]]$y-fitlist.1D[[j]]$fitted.values),
                 mean(fitlist.1D.G[[j]]$y-fitlist.1D.G[[j]]$fitted.values),
                 sd(fitlist.1D.G[[j]]$y-fitlist.1D.G[[j]]$fitted.values))
}

sum(resid$SD.Score>resid$SD.Theta)
sum(resid$Mean.Score>resid$Mean.Theta)
sum(resid$SD.Score>resid$SD.Theta.G)
sum(resid$Mean.Score>resid$Mean.Theta.G)
sum(resid$SD.Theta>resid$SD.Theta.G)
sum(resid$Mean.Theta>resid$Mean.Theta.G)

for (r in 1:nrow(resid)) cat(paste0("$",r,"$ & $",paste(signif(resid[r,],digits=3),collapse="$ & $"),"$ \\\\ \\hline \n"))

resid.saem = data.frame(Mean.Score=numeric(),SD.Score=numeric(),Mean.Theta=numeric(),SD.Theta=numeric(),Mean.Theta.G=numeric(),SD.Theta.G=numeric())
for (j in 1:ncol(Response)) {
  resid.saem[j,]<-c(mean(fitlist.score[[j]]$y-fitlist.score[[j]]$fitted.values),
               sd(fitlist.score[[j]]$y-fitlist.score[[j]]$fitted.values),
               mean(fitlist.1D[[j]]$y - ProbOgive(xi=FCI.1D.NoGuess$xi,
                                                  theta = fitlist.1D[[j]]$data$x.data,guess = FALSE,j = j)),
               sd(fitlist.1D[[j]]$y - ProbOgive(xi=FCI.1D.NoGuess$xi,
                                                theta = fitlist.1D[[j]]$data$x.data,guess = FALSE,j = j)),
               mean(fitlist.1D.G[[j]]$y - ProbOgive(xi=FCI.1D.Guess$xi,
                                                    theta = fitlist.1D[[j]]$data$x.data,guess = TRUE,j = j)),
               sd(fitlist.1D.G[[j]]$y - ProbOgive(xi=FCI.1D.Guess$xi,
                                                  theta = fitlist.1D[[j]]$data$x.data,guess = TRUE,j = j)))
}

sum(resid.saem$SD.Score>resid.saem$SD.Theta)
sum(resid.saem$Mean.Score>resid.saem$Mean.Theta)
sum(resid.saem$SD.Score>resid.saem$SD.Theta.G)
sum(resid.saem$Mean.Score>resid.saem$Mean.Theta.G)
sum(resid.saem$SD.Theta>resid.saem$SD.Theta.G)
sum(resid.saem$Mean.Theta>resid.saem$Mean.Theta.G)


load("RealData/FCI_NoGuess_A1.rda")
GetLikelihood(FitDATA)
FCI.1D.NoGuess<-FitDATA
#     LogL       AIC       BIC       DOF 
# -11390.13  22900.27  23177.31     60.00 
load("RealData/FCI_Guess_A1.rda")
GetLikelihood(FitDATA)
FCI.1D.Guess<-FitDATA
#     LogL       AIC       BIC       DOF 
# -11325.46  22830.93  23246.50     90.00 
LRTest(FCI.1D.Guess,FCI.1D.NoGuess)
# "Ratio = 129.337377235424 : DoF = 30 : P = 2.67563748934663e-14"

load("RealData/FCI_NoGuess_A2.rda")
GetLikelihood(FitDATA)
FCI.2D.NoGuess<-FitDATA
#     LogL       AIC       BIC       DOF 
# -10611.83  21403.66  21819.22     90.00 
LRTest(FCI.1D.NoGuess,FCI.2D.NoGuess)
# "Ratio = 1556.61022200324 : DoF = 30 : P = 0"

load("RealData/FCI_Guess_A2.rda")
GetLikelihood(FitDATA)
FCI.2D.Guess<-FitDATA
#     LogL       AIC       BIC       DOF 
# -10544.31  21328.63  21882.71    120.00 
LRTest(FCI.1D.Guess,FCI.2D.Guess)
#"Ratio = 1562.30376406312 : DoF = 30 : P = 0"
LRTest(FCI.2D.NoGuess,FCI.2D.Guess)
# "Ratio = 135.030919295306 : DoF = 30 : P = 2.77555756156289e-15"

load("RealData/FCI_NoGuess_A3.rda")
GetLikelihood(FitDATA)
FCI.3D.NoGuess<-FitDATA
#     LogL       AIC       BIC       DOF 
# -9991.786 20223.571 20777.660   120.000 
LRTest(FCI.2D.NoGuess,FCI.3D.NoGuess)
# "Ratio = 1240.08528262367 : DoF = 30 : P = 0

load("RealData/FCI_Guess_A3.rda")
GetLikelihood(FitDATA)
FCI.3D.Guess<-FitDATA
#    LogL       AIC       BIC       DOF 
# -10164.98  20629.96  21322.57    150.00 
LRTest(FCI.2D.Guess,FCI.3D.Guess)
#"Ratio = 758.662203360804 : DoF = 30 : P = 0"
LRTest(FCI.3D.NoGuess,FCI.3D.Guess)
#"Ratio = -346.392159967556 : DoF = 30 : P = 1"

load("RealData/FCI_NoGuess_A4.rda")
GetLikelihood(FitDATA)
FCI.4D.NoGuess<-FitDATA
#    LogL       AIC       BIC       DOF 
# -9547.221 19394.441 20087.052   150.000 
LRTest(FCI.3D.NoGuess,FCI.4D.NoGuess)
# "Ratio = 889.130257364344 : DoF = 30 : P = 0"

# load("RealData/FCI_Guess_A4.rda")
# GetLikelihood(FitDATA)
# FCI.4D.Guess<-FitDATA
# LRTest(FCI.3D.Guess,FCI.4D.Guess)
# LRTest(FCI.4D.NoGuess,FCI.4D.Guess)

for (r in 1:nrow(FCI.1D.NoGuess$xi)) cat(paste0("$",r,"$ & $",
                                                paste(c(round(c(1,1,1/1.7)*FCI.1D.NoGuess$xi[r,c(1,2,2)],digits=3),
                                                        round(c(1,1,1/1.7,1)*FCI.1D.Guess$xi[r,c(1,2,2,3)],digits=3)),collapse="$ & $"),
                                                "$ \\\\ \\hline \n"))

t(apply(FCI.1D.NoGuess$xi,1,function(x) signif(x*c(1,1),3)))

# WANG: A: 13, 18, 5, 25, 26, 30  Low B: 6, 1, 12, 24, 3 High B: 17, 26, 15, 25, 5  Guess: 16, 9==22,  6, 7, 27, 21 
# CHEN: A: 18, 13, 5, 25, 11, 30  Low B: 1, 3, 19, 12, 7 High B: 15, 17 ,26, 25, 4  Guess: 6, 16, 9, 8, 24 
# GEIS: A: 13, 5, 18, 2, 26, 25   Low B: 12, 6, 24, 10, 1  High B: 5, 13, 2, 18, 26  Guess: 1, 9, 19, 6, 7

########### CCI

#source("WIDER_DATA/InitiateData.R")
# CCIC<-paste("C",94:115,sep="")
# CCI.Response<-IDF[!is.na(IDF$CCI.x)&!is.na(IDF$LetterGrade159),CCIC]
# write.csv(CCI.Response,paste0(data.dir,"/CCIPreTest.csv"),row.names=FALSE)

Response<-read.csv(paste0(data.dir,"/CCIPreTest.csv"),stringsAsFactors = FALSE,header = TRUE) # N=628
settings<-CheckParams(generate = FALSE)
settings$Adim<-4
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
CCI.Fit.4D.Guess<-AnalyzeTestData(RP = Response,settings = settings,verbose = TRUE,TargetA = NA)
GetLikelihood(CCI.Fit.4D.Guess)
TWFitTest(fit.data = CCI.Fit.4D.Guess)

settings$Adim<-3
settings$tmu<-rep(0,settings$Adim)
settings$tsigma<-diag(settings$Adim)
settings$estfile<-paste0(data.dir,"/CCI_Guess_A",settings$Adim)
CCI.Fit.3D.Guess<-AnalyzeTestData(RP = Response,settings = settings,verbose = TRUE,TargetA = NA)
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

settings$Adim<-1
settings$guess<-TRUE
settings$tmu<-rep(0,settings$Adim)
settings$tsigma<-diag(settings$Adim)
settings$estfile<-paste0(data.dir,"/CCI_Guess_A",settings$Adim)
CCI.Fit.1D.Guess<-AnalyzeTestData(RP = Response,settings = settings,verbose = TRUE,TargetA = NA)
GetLikelihood(CCI.Fit.1D.Guess)
TWFitTest(fit.data = CCI.Fit.1D.Guess)

settings$Adim<-1
settings$guess<-FALSE;settings$tmu<-rep(0,settings$Adim);settings$tsigma<-diag(settings$Adim)
settings$estfile<-paste0(data.dir,"/CCI_NoGuess_A",settings$Adim)
CCI.Fit.1D.NoGuess<-AnalyzeTestData(RP = Response,settings = settings,verbose = TRUE,TargetA = NA)
GetLikelihood(CCI.Fit.1D.NoGuess)
TWFitTest(fit.data = CCI.Fit.1D.NoGuess)

settings$Adim<-2;settings$tmu<-rep(0,settings$Adim);settings$tsigma<-diag(settings$Adim)
settings$estfile<-paste0(data.dir,"/CCI_NoGuess_A",settings$Adim)
CCI.Fit.2D.NoGuess<-AnalyzeTestData(RP = Response,settings = settings,verbose = TRUE,TargetA = NA)
GetLikelihood(CCI.Fit.2D.NoGuess)
TWFitTest(fit.data = CCI.Fit.2D.NoGuess)

settings$Adim<-3;settings$tmu<-rep(0,settings$Adim);settings$tsigma<-diag(settings$Adim)
settings$estfile<-paste0(data.dir,"/CCI_NoGuess_A",settings$Adim)
CCI.Fit.3D.NoGuess<-AnalyzeTestData(RP = Response,settings = settings,verbose = TRUE,TargetA = NA)
GetLikelihood(CCI.Fit.3D.NoGuess)
TWFitTest(fit.data = CCI.Fit.3D.NoGuess,tw=0.95)

settings$Adim<-4;settings$tmu<-rep(0,settings$Adim);settings$tsigma<-diag(settings$Adim)
settings$estfile<-paste0(data.dir,"/CCI_NoGuess_A",settings$Adim)
CCI.Fit.4D.NoGuess<-AnalyzeTestData(RP = Response,settings = settings,verbose = TRUE,TargetA = NA)
GetLikelihood(CCI.Fit.4D.NoGuess)
TWFitTest(fit.data = CCI.Fit.4D.NoGuess)

settings$Adim<-5;settings$tmu<-rep(0,settings$Adim);settings$tsigma<-diag(settings$Adim)
settings$estfile<-paste0(data.dir,"/CCI_NoGuess_A",settings$Adim)
CCI.Fit.5D.NoGuess<-AnalyzeTestData(RP = Response,settings = settings,verbose = TRUE,TargetA = NA)
#GetLikelihood(CCI.Fit.5D.NoGuess)
TWFitTest(fit.data = CCI.Fit.5D.NoGuess)
LRTest(CCI.Fit.5D.NoGuess,CCI.Fit.4D.NoGuess)
LRTest(CCI.Fit.5D.NoGuess,CCI.Fit.2D.NoGuess)

settings$Adim<-6; settings$tmu<-rep(0,settings$Adim); settings$tsigma<-diag(settings$Adim)
settings$estfile<-paste0(data.dir,"/CCI_NoGuess_A",settings$Adim)
CCI.Fit.6D.NoGuess<-AnalyzeTestData(RP = Response,settings = settings,verbose = TRUE,TargetA = NA)
#GetLikelihood(CCI.Fit.6D.NoGuess)
#TWFitTest(fit.data = CCI.Fit.6D.NoGuess)
LRTest(CCI.Fit.6D.NoGuess,CCI.Fit.5D.NoGuess)
LRTest(CCI.Fit.6D.NoGuess,CCI.Fit.3D.NoGuess)

settings$Adim<-7; settings$tmu<-rep(0,settings$Adim); settings$tsigma<-diag(settings$Adim)
settings$estfile<-paste0(data.dir,"/CCI_NoGuess_A",settings$Adim)
CCI.Fit.7D.NoGuess<-AnalyzeTestData(RP = Response,settings = settings,verbose = TRUE,TargetA = NA)
#GetLikelihood(CCI.Fit.7D.NoGuess)
TWFitTest(fit.data = CCI.Fit.7D.NoGuess)
LRTest(CCI.Fit.7D.NoGuess,CCI.Fit.6D.NoGuess)
LRTest(CCI.Fit.7D.NoGuess,CCI.Fit.4D.NoGuess)


settings$Adim<-8; settings$tmu<-rep(0,settings$Adim); settings$tsigma<-diag(settings$Adim)
settings$estfile<-paste0(data.dir,"/CCI_NoGuess_A",settings$Adim)
CCI.Fit.8D.NoGuess<-AnalyzeTestData(RP = Response,settings = settings,verbose = TRUE,TargetA = NA)
#GetLikelihood(CCI.Fit.7D.NoGuess)
TWFitTest(fit.data = CCI.Fit.8D.NoGuess)
LRTest(CCI.Fit.8D.NoGuess,CCI.Fit.7D.NoGuess)
LRTest(CCI.Fit.8D.NoGuess,CCI.Fit.6D.NoGuess)

settings$Adim<-9; settings$tmu<-rep(0,settings$Adim); settings$tsigma<-diag(settings$Adim)
settings$estfile<-paste0(data.dir,"/CCI_NoGuess_A",settings$Adim)
CCI.Fit.9D.NoGuess<-AnalyzeTestData(RP = Response,settings = settings,verbose = TRUE,TargetA = NA)
#GetLikelihood(CCI.Fit.7D.NoGuess)
TWFitTest(fit.data = CCI.Fit.9D.NoGuess)
LRTest(CCI.Fit.9D.NoGuess,CCI.Fit.8D.NoGuess)
LRTest(CCI.Fit.9D.NoGuess,CCI.Fit.7D.NoGuess)

settings$Adim<-10; settings$tmu<-rep(0,settings$Adim); settings$tsigma<-diag(settings$Adim)
settings$estfile<-paste0(data.dir,"/CCI_NoGuess_A",settings$Adim)
CCI.Fit.10D.NoGuess<-AnalyzeTestData(RP = Response,settings = settings,verbose = TRUE,TargetA = NA)
#GetLikelihood(CCI.Fit.10D.NoGuess)
TWFitTest(fit.data = CCI.Fit.10D.NoGuess)
LRTest(CCI.Fit.10D.NoGuess,CCI.Fit.9D.NoGuess)
LRTest(CCI.Fit.10D.NoGuess,CCI.Fit.8D.NoGuess)

settings$Adim<-11; settings$tmu<-rep(0,settings$Adim); settings$tsigma<-diag(settings$Adim)
settings$estfile<-paste0(data.dir,"/CCI_NoGuess_A",settings$Adim)
CCI.Fit.11D.NoGuess<-AnalyzeTestData(RP = Response,settings = settings,verbose = TRUE,TargetA = NA)
#GetLikelihood(CCI.Fit.10D.NoGuess)
TWFitTest(fit.data = CCI.Fit.11D.NoGuess)
LRTest(CCI.Fit.11D.NoGuess,CCI.Fit.10D.NoGuess)
LRTest(CCI.Fit.11D.NoGuess,CCI.Fit.9D.NoGuess)

settings$Adim<-12; settings$tmu<-rep(0,settings$Adim); settings$tsigma<-diag(settings$Adim)
settings$estfile<-paste0(data.dir,"/CCI_NoGuess_A",settings$Adim)
CCI.Fit.12D.NoGuess<-AnalyzeTestData(RP = Response,settings = settings,verbose = TRUE,TargetA = NA)
#GetLikelihood(CCI.Fit.10D.NoGuess)
TWFitTest(fit.data = CCI.Fit.12D.NoGuess)
LRTest(CCI.Fit.12D.NoGuess,CCI.Fit.11D.NoGuess)
LRTest(CCI.Fit.12D.NoGuess,CCI.Fit.10D.NoGuess)

settings$Adim<-13; settings$tmu<-rep(0,settings$Adim); settings$tsigma<-diag(settings$Adim)
settings$estfile<-paste0(data.dir,"/CCI_NoGuess_A",settings$Adim)
CCI.Fit.13D.NoGuess<-AnalyzeTestData(RP = Response,settings = settings,verbose = TRUE,TargetA = NA)
#GetLikelihood(CCI.Fit.10D.NoGuess)
TWFitTest(fit.data = CCI.Fit.13D.NoGuess)
LRTest(CCI.Fit.13D.NoGuess,CCI.Fit.12D.NoGuess)

settings$Adim<-14; settings$tmu<-rep(0,settings$Adim); settings$tsigma<-diag(settings$Adim)
settings$estfile<-paste0(data.dir,"/CCI_NoGuess_A",settings$Adim)
CCI.Fit.14D.NoGuess<-AnalyzeTestData(RP = Response,settings = settings,verbose = TRUE,TargetA = NA)
#GetLikelihood(CCI.Fit.10D.NoGuess)
TWFitTest(fit.data = CCI.Fit.14D.NoGuess)
LRTest(CCI.Fit.14D.NoGuess,CCI.Fit.13D.NoGuess)

settings$Adim<-15; settings$tmu<-rep(0,settings$Adim); settings$tsigma<-diag(settings$Adim)
settings$estfile<-paste0(data.dir,"/CCI_NoGuess_A",settings$Adim)
CCI.Fit.15D.NoGuess<-AnalyzeTestData(RP = Response,settings = settings,verbose = TRUE,TargetA = NA)
#GetLikelihood(CCI.Fit.10D.NoGuess)
TWFitTest(fit.data = CCI.Fit.15D.NoGuess)
LRTest(CCI.Fit.15D.NoGuess,CCI.Fit.14D.NoGuess)

source("WIDER_DATA/InitiateData.R")
source("WIDER_DATA/ItemFunctions.R")

WIDER.data<-read.csv("WIDER_DATA/IDF.csv",header=TRUE)
ItemStats(WIDER.data,it=64)

CCIC<-paste("C",94:115,sep="")
Response<-IDF[which(!is.na(IDF$CCI.x)&!is.na(IDF$LetterGrade159)),CCIC]
CCIB<-paste("B",94:115,sep="")
b.RP<-IDF[which(!is.na(IDF$CCI.x)&!is.na(IDF$LetterGrade159)),CCIB]
Correct<-CHPreTest[8:29,"Numbers"]
scores<-table(100*round(rowMeans(Response),digits = 3))

par(mfrow=c(1,2),mar=c(5,4,3,3))
barplot(scores,main="CCI Pretest Score Distribution",names.arg = names(scores),xlab = "% Correct",ylab = "Frequency",las=2)
qqnorm(y=100*round(rowMeans(Response),digits = 3),main="CCI Pretest Q-Q Analysis")
fitlist.score<-ResponseCurves(responses = b.RP[5:22],scores = rowMeans(Response),
                              prows = 3,pcols = 3,correct = Correct[5:22],j.legend=15)
fitlist.1D<-ResponseCurves(responses = b.RP,scores = CCI.1D.NoGuess$That[,"Theta"],prows = 3,pcols = 5,correct = Correct)
fitlist.1D.G<-ResponseCurves(responses = b.RP,scores = CCI.1D.Guess$That[,"Theta"],prows = 3,pcols = 5,correct = Correct)
par(mfrow=c(1,4)) 
resid = data.frame(Mean.Score=numeric(),SD.Score=numeric(),Mean.Theta=numeric(),SD.Theta=numeric(),Mean.Theta.G=numeric(),SD.Theta.G=numeric())
for (j in 1:ncol(Response)) {
  plot(density(fitlist.score[[j]]$y-fitlist.score[[j]]$fitted.values),main="Residual Comparisons",
       xlab=expression("P(Y=1) -"~hat(P)(Y=1)),lty=2,ylab="density",
       ylim=range(c(density(fitlist.1D[[j]]$y-fitlist.1D[[j]]$fitted.values)$y,density(fitlist.score[[j]]$y-fitlist.score[[j]]$fitted.values)$y)))
  lines(density(fitlist.1D[[j]]$y-fitlist.1D[[j]]$fitted.values),lty=1,col=2)
  resid[j,]<-c(mean(fitlist.score[[j]]$y-fitlist.score[[j]]$fitted.values),
               sd(fitlist.score[[j]]$y-fitlist.score[[j]]$fitted.values),
               mean(fitlist.1D[[j]]$y-fitlist.1D[[j]]$fitted.values),
               sd(fitlist.1D[[j]]$y-fitlist.1D[[j]]$fitted.values),
               mean(fitlist.1D.G[[j]]$y-fitlist.1D.G[[j]]$fitted.values),
               sd(fitlist.1D.G[[j]]$y-fitlist.1D.G[[j]]$fitted.values))
}



### QOL Data
?Rprof
Rprof()
Response = scan(paste0(data.dir,"/qol.dat"), what = "numeric",sep = "\n") # N=753
Response = data.matrix(do.call(rbind,lapply(strsplit(Response,"\\s"), function(x) as.numeric(x[-1]))))
settings<-CheckParams(generate = FALSE)
settings$ncat<-5
settings$Adim<-5
settings$dbltrunc<-TRUE
settings$tmu<-rep(0,settings$Adim)
settings$tsigma<-diag(settings$Adim)
settings$rmethod<-NA
settings$empiricalse<-FALSE
settings$esttheta<-FALSE
settings$nesttheta<-200
settings$thinA=8
settings$thinB=5
settings$EmpIT=1000
settings$estfile<-paste0(data.dir,"/QOL_TIMED_A",settings$Adim)
settings$burnin=2000
settings$eps=0.0001
settings$thetamap=FALSE
settings$record=FALSE
settings$cores<-1
settings$parallel<-FALSE
settings$burnin=80
QOL.Fit.1D<-AnalyzeTestData(RP = Response,settings = settings,verbose = TRUE,TargetA = NA,timed=TRUE)
Rprof(NULL)
summaryRprof()
GetLikelihood(QOL.Fit.1D)
#TWFitTest(fit.data = QOL.Fit.1D)

settings$Adim<-2; settings$tmu<-rep(0,settings$Adim); settings$tsigma<-diag(settings$Adim)
settings$estfile<-paste0(data.dir,"/QOL_A",settings$Adim)
QOL.Fit.2D<-AnalyzeTestData(RP = Response,settings = settings,verbose = TRUE,TargetA = NA)
GetLikelihood(QOL.Fit.2D)
LRTest(QOL.Fit.2D,QOL.Fit.1D)
#TWFitTest(fit.data = QOL.Fit.2D)

settings$Adim<-3; settings$tmu<-rep(0,settings$Adim); settings$tsigma<-diag(settings$Adim)
settings$estfile<-paste0(data.dir,"/QOL_A",settings$Adim)
QOL.Fit.3D<-AnalyzeTestData(RP = Response,settings = settings,verbose = TRUE,TargetA = NA)
GetLikelihood(QOL.Fit.3D)
LRTest(QOL.Fit.2D,QOL.Fit.3D)
#TWFitTest(fit.data = QOL.Fit.3D)

settings$Adim<-4; settings$tmu<-rep(0,settings$Adim); settings$tsigma<-diag(settings$Adim)
settings$estfile<-paste0(data.dir,"/QOL_A",settings$Adim)
QOL.Fit.4D<-AnalyzeTestData(RP = Response,settings = settings,verbose = TRUE,TargetA = NA)
GetLikelihood(QOL.Fit.4D)
LRTest(QOL.Fit.4D,QOL.Fit.3D)
#TWFitTest(fit.data = QOL.Fit.4D)

settings$Adim<-5; settings$tmu<-rep(0,settings$Adim); settings$tsigma<-diag(settings$Adim)
settings$estfile<-paste0(data.dir,"/QOL_A",settings$Adim)
QOL.Fit.5D<-AnalyzeTestData(RP = Response,settings = settings,verbose = TRUE,TargetA = NA)
GetLikelihood(QOL.Fit.5D)
LRTest(QOL.Fit.5D,QOL.Fit.4D)
#TWFitTest(fit.data = QOL.Fit.5D)

# settings$Adim<-6; settings$tmu<-rep(0,settings$Adim); settings$tsigma<-diag(settings$Adim)
# settings$estfile<-paste0(data.dir,"/QOL_A",settings$Adim)
# QOL.Fit.6D<-AnalyzeTestData(RP = Response,settings = settings,verbose = TRUE,TargetA = NA)
# #GetLikelihood(QOL.Fit.6D)
# TWFitTest(fit.data = QOL.Fit.6D)
# 
# settings$Adim<-7; settings$tmu<-rep(0,settings$Adim); settings$tsigma<-diag(settings$Adim)
# settings$estfile<-paste0(data.dir,"/QOL_A",settings$Adim)
# QOL.Fit.7D<-AnalyzeTestData(RP = Response,settings = settings,verbose = TRUE,TargetA = NA)
# GetLikelihood(QOL.Fit.7D)
# TWFitTest(fit.data = QOL.Fit.7D)
# 
# settings$Adim<-10; settings$tmu<-rep(0,settings$Adim); settings$tsigma<-diag(settings$Adim)
# settings$estfile<-paste0(data.dir,"/QOL_A",settings$Adim)
# QOL.Fit.10D<-AnalyzeTestData(RP = Response,settings = settings,verbose = TRUE,TargetA = NA)
# GetLikelihood(QOL.Fit.10D)
# TWFitTest(fit.data = QOL.Fit.10D)



# source("WIDER_DATA/Demographics.R")
# source("WIDER_DATA/LatentClass.R")
# source("WIDER_DATA/Predictions.R")
# source("WIDER_DATA/IRT.R")
# source("WIDER_DATA/Diagnostics.R")
# source("WIDER_DATA/Options.R")
# source("WIDER_DATA/CLASS.R")
# source("WIDER_DATA/GFA.R")
# 
# FCIC<-paste("C",64:93,sep="")
# Response<-IDF[which(!is.na(IDF$FCI.x)&(!is.na(IDF$LetterGrade123)|!is.na(IDF$LetterGrade115))),FCIC]
# 
# IRTFit<-ApplyIRT(RP=Response,model="ogive")
# GFAFit<-ApplyGFA(RP=Response,GParam=FALSE,AParams=1,HDQs=NA,Plot=FALSE)
# points(1:60,c(IRTFit$A,IRTFit$B),col=6,cex=2)
# 
# 
# 
# CCols<-paste("C",1:129,sep="")
# QPCols<-paste("QP",1:42,sep="")
# QCCols<-paste("QC",1:50,sep="")
# Cols<-c(CCols,QPCols,QCCols)
# 
# plot(density(IDF$CSEM.y[which(!is.na(IDF$Grade116))],na.rm=TRUE))
# # lines(density(IDF$CSEM.y[which(is.na(IDF$Grade116))],na.rm=TRUE),col=2)
# # plot(density(IDF$CSEM.x[which(!is.na(IDF$Grade116))],na.rm=TRUE))
# # lines(density(IDF$CSEM.x[which(is.na(IDF$Grade116))],na.rm=TRUE),col=2)