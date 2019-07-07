source("InitializeGeisCamilli.R")
basedir<-getwd()
data.dir<-"RealData"
if (!data.dir %in% dir()) dir.create(data.dir)

load("RealData/FCI_NoGuess_A1.rda")
FCI.Fit.1D.NoGuess<-FitDATA
load("RealData/FCI_NoGuess_A2.rda")
FCI.Fit.2D.NoGuess<-FitDATA
load("RealData/FCI_NoGuess_A3.rda")
FCI.Fit.3D.NoGuess<-FitDATA
load("RealData/FCI_NoGuess_A4.rda")
FCI.Fit.4D.NoGuess<-FitDATA
# load("RealData/FCI_NoGuess_A10.rda")
# FCI.Fit.10D.NoGuess<-FitDATA
rm(FitDATA)
par(mfrow=c(1,2))
MultipleTWFitTests(fit.data.list = list(FCI.Fit.1D.NoGuess,
                                        FCI.Fit.2D.NoGuess,
                                        FCI.Fit.3D.NoGuess,
                                        FCI.Fit.4D.NoGuess),
                   title = "FCI (no guessing)",plot.layout = FALSE,
                   cex=0.95,ratios = FALSE)

load("RealData/FCI_Guess_A1.rda")
FCI.Fit.1D.Guess<-FitDATA
load("RealData/FCI_Guess_A2.rda")
FCI.Fit.2D.Guess<-FitDATA
load("RealData/FCI_Guess_A3.rda")
FCI.Fit.3D.Guess<-FitDATA
load("RealData/FCI_Guess_A4.rda")
FCI.Fit.4D.Guess<-FitDATA
# load("RealData/FCI_NoGuess_A10.rda")
# FCI.Fit.10D.NoGuess<-FitDATA
rm(FitDATA)
MultipleTWFitTests(fit.data.list = list(FCI.Fit.1D.Guess,
                                        FCI.Fit.2D.Guess,
                                        FCI.Fit.3D.Guess,
                                        FCI.Fit.4D.Guess),
                   title = "FCI (guessing)",plot.layout = FALSE,
                   cex=0.95,ratios = FALSE)

#########  Rotation of 2-D

A.Rot<-array(rep(0,360*length(FCI.Fit.2D.Guess$A)),
             dim = c(nrow(FCI.Fit.2D.Guess$A),
                     ncol(FCI.Fit.2D.Guess$A),360))
V.L<-matrix(rep(0,360*nrow(FCI.Fit.2D.Guess$A)),
            nrow=ncol(FCI.Fit.2D.Guess$A),
            ncol=360)
V.V<-matrix(rep(0,360*nrow(FCI.Fit.2D.Guess$A)),
            nrow=ncol(FCI.Fit.2D.Guess$A),
            ncol=360)
for (x in 1:360) {
  r<-(x/180)*pi
  A.Rot[,,x] <- FCI.Fit.2D.Guess$A %*% matrix(c(cos(r),-sin(r),sin(r),cos(r)),
                                              nrow=2,ncol=2)
  V.L[,x]<-apply(A.Rot[,,x],2,function(y) sqrt(sum(y^2)))
  V.V[,x]<-apply(A.Rot[,,x],2,function(y) var(y^2))
}
A.items<-unique(as.vector(apply(abs(A.Rot),c(2,3),which.max)))
par(mfrow=c(3,1),mar=c(1,4,3,2))
plot(1:ncol(V.L),apply(V.L,2,max),
     type="n",main="Rotated Discrimination",
     ylab="L2 Norm",xaxt="n",ylim=c(0,max(V.L)))
lines(1:ncol(V.L),V.L[1,],lty=5)
lines(1:ncol(V.L),V.L[2,],lty=3)
abline(v=43.5)
par(mar=c(3,4,2,2))
plot(1:ncol(V.V),apply(V.V,2,max),
     type="n",main=NA,
     ylab="Variance",xaxt="n",ylim=c(0,max(V.V)))
lines(1:ncol(V.V),V.V[1,],lty=5)
lines(1:ncol(V.V),V.V[2,],lty=3)
abline(v=43.5)
par(mar=c(5,4,0,2))
plot(c(1,dim(A.Rot)[3]),range(A.Rot),
     type="n",main=NA,
     ylab="Discrimination",xlab="Rotation (degrees)")
abline(h=0)
for (j in 1:length(A.items)) {
  lines(1:dim(A.Rot)[3],A.Rot[A.items[j],1,],lty=5,col=j+1)
  lines(1:dim(A.Rot)[3],A.Rot[A.items[j],2,],lty=3,col=j+1)
}
legend("top",paste("Item",A.items),col=1:length(A.items)+1,lty=1)
abline(v=43.5)

for (j in 1:30) cat(j," & ",paste0(round(A.Rot[j,,43],digits=3),collapse=" & ")," \\\\ \\hline \n")

##########################  Without Guessing
A.Rot<-array(rep(0,360*length(FCI.Fit.2D.NoGuess$A)),
             dim = c(nrow(FCI.Fit.2D.NoGuess$A),
                     ncol(FCI.Fit.2D.NoGuess$A),360))
V.L<-matrix(rep(0,360*nrow(FCI.Fit.2D.NoGuess$A)),
            nrow=ncol(FCI.Fit.2D.NoGuess$A),
            ncol=360)
V.V<-matrix(rep(0,360*nrow(FCI.Fit.2D.NoGuess$A)),
            nrow=ncol(FCI.Fit.2D.NoGuess$A),
            ncol=360)
for (x in 1:360) {
  r<-(x/180)*pi
  A.Rot[,,x] <- FCI.Fit.2D.NoGuess$A %*% matrix(c(cos(r),-sin(r),sin(r),cos(r)),
                                                nrow=2,ncol=2)
  V.L[,x]<-apply(A.Rot[,,x],2,function(y) sqrt(sum(y^2)))
  V.V[,x]<-apply(A.Rot[,,x],2,function(y) var(y^2))
}
A.items<-unique(as.vector(apply(abs(A.Rot),c(2,3),which.max)))
par(mfrow=c(3,1),mar=c(1,4,3,2))
plot(1:ncol(V.L),apply(V.L,2,max),
     type="n",main="Rotated Discrimination",
     ylab="L2 Norm",xaxt="n",ylim=c(0,max(V.L)))
lines(1:ncol(V.L),V.L[1,],lty=5)
lines(1:ncol(V.L),V.L[2,],lty=3)
abline(v=43.5)
par(mar=c(3,4,2,2))
plot(1:ncol(V.V),apply(V.V,2,max),
     type="n",main=NA,
     ylab="Variance",xaxt="n",ylim=c(0,max(V.V)))
lines(1:ncol(V.V),V.V[1,],lty=5)
lines(1:ncol(V.V),V.V[2,],lty=3)
abline(v=43.5)
par(mar=c(5,4,0,2))
plot(c(1,90),range(A.Rot[,1,]/(abs(A.Rot[,2,])+0.1)),
     type="n",main=NA,
     ylab="Discrimination",xlab="Rotation (degrees)",ylim=c(-1,15))
abline(h=0)
for (j in 1:30) {
  lines(1:90,A.Rot[j,1,1:90]/(abs(A.Rot[j,2,1:90])+0.1),lty=5,col=j+1)
  lines(1:90,A.Rot[j,2,1:90]/(abs(A.Rot[j,1,1:90])+0.1),lty=3,col=j+1)
}
legend("top",paste("Item",A.items),col=1:length(A.items)+1,lty=1)
abline(v=43.5)



load("RealData/CCI_NoGuess_A1.rda")
CCI.Fit.1D.NoGuess<-FitDATA
load("RealData/CCI_NoGuess_A2.rda")
CCI.Fit.2D.NoGuess<-FitDATA
load("RealData/CCI_NoGuess_A3.rda")
CCI.Fit.3D.NoGuess<-FitDATA
load("RealData/CCI_NoGuess_A4.rda")
CCI.Fit.4D.NoGuess<-FitDATA
load("RealData/CCI_NoGuess_A5.rda")
CCI.Fit.5D.NoGuess<-FitDATA
load("RealData/CCI_NoGuess_A10.rda")
CCI.Fit.10D.NoGuess<-FitDATA
rm(FitDATA)

MultipleTWFitTests(fit.data.list = list(CCI.Fit.1D.NoGuess,
                                        CCI.Fit.2D.NoGuess,
                                        CCI.Fit.3D.NoGuess,
                                        CCI.Fit.4D.NoGuess),title = "CCI",ratios = FALSE)

MultipleTWFitTests(fit.data.list = list(CCI.Fit.1D.Guess,
                                        CCI.Fit.2D.Guess,
                                        CCI.Fit.3D.Guess,
                                        CCI.Fit.4D.Guess),title = "CCI",ratios = FALSE)

load("RealData/CCI_Guess_A1.rda")
CCI.Fit.1D.Guess<-FitDATA
load("RealData/CCI_Guess_A2.rda")
CCI.Fit.2D.Guess<-FitDATA
load("RealData/CCI_Guess_A3.rda")
CCI.Fit.3D.Guess<-FitDATA
load("RealData/CCI_Guess_A4.rda")
CCI.Fit.4D.Guess<-FitDATA

LRTest(CCI.Fit.4D.Guess,CCI.Fit.3D.Guess)
LRTest(CCI.Fit.3D.Guess,CCI.Fit.2D.Guess)
LRTest(CCI.Fit.2D.Guess,CCI.Fit.1D.Guess)
LRTest(CCI.Fit.1D.Guess,CCI.Fit.1D.NoGuess)
LRTest(CCI.Fit.2D.Guess,CCI.Fit.2D.NoGuess)
LRTest(CCI.Fit.2D.NoGuess,CCI.Fit.1D.NoGuess)
LRTest(CCI.Fit.3D.Guess,CCI.Fit.3D.NoGuess)
LRTest(CCI.Fit.3D.NoGuess,CCI.Fit.2D.NoGuess)
LRTest(CCI.Fit.4D.NoGuess,CCI.Fit.3D.NoGuess)


A.Rot<-array(rep(0,360*length(CCI.Fit.2D.NoGuess$A)),
             dim = c(nrow(CCI.Fit.2D.NoGuess$A),
                     ncol(CCI.Fit.2D.NoGuess$A),360))
V.L<-matrix(rep(0,360*nrow(CCI.Fit.2D.NoGuess$A)),
            nrow=ncol(CCI.Fit.2D.NoGuess$A),
            ncol=360)
V.V<-matrix(rep(0,360*nrow(CCI.Fit.2D.NoGuess$A)),
            nrow=ncol(CCI.Fit.2D.NoGuess$A),
            ncol=360)
for (x in 1:360) {
  r<-(x/180)*pi
  A.Rot[,,x] <- CCI.Fit.2D.NoGuess$A %*% matrix(c(cos(r),-sin(r),sin(r),cos(r)),
                                                nrow=2,ncol=2)
  V.L[,x]<-apply(A.Rot[,,x],2,function(y) sqrt(sum(y^2)))
  V.V[,x]<-apply(A.Rot[,,x],2,function(y) var(y^2))
}
A.items<-unique(as.vector(apply(abs(A.Rot),c(2,3),function(x) (order(-x)[1:3]))))
par(mfrow=c(3,1),mar=c(1,4,3,2))
plot(1:ncol(V.L),apply(V.L,2,max),
     type="n",main="Rotated Discrimination",
     ylab="L2 Norm",xaxt="n",ylim=c(0,max(V.L)))
lines(1:ncol(V.L),V.L[1,],lty=5)
lines(1:ncol(V.L),V.L[2,],lty=3)
abline(v=10)
par(mar=c(3,4,2,2))
plot(1:ncol(V.V),apply(V.V,2,max),
     type="n",main=NA,
     ylab="Variance",xaxt="n",ylim=c(0,max(V.V)))
lines(1:ncol(V.V),V.V[1,],lty=5)
lines(1:ncol(V.V),V.V[2,],lty=3)
abline(v=10)
par(mar=c(5,4,0,2))
plot(c(1,dim(A.Rot)[3]),range(A.Rot),
     type="n",main=NA,
     ylab="Discrimination",xlab="Rotation (degrees)")
abline(h=0)
for (j in 1:length(A.items)) {
  lines(1:dim(A.Rot)[3],A.Rot[A.items[j],1,],lty=5,col=j+1)
  lines(1:dim(A.Rot)[3],A.Rot[A.items[j],2,],lty=3,col=j+1)
}
legend("top",paste("Item",A.items),col=1:length(A.items)+1,lty=1)
abline(v=10)

for (j in 1:22) cat(j," & ",paste0(round(A.Rot[j,,10],digits=3),collapse=" & ")," \\\\ \\hline \n")

ROT3<-RotateSlopes(fit.data = CCI.Fit.3D.NoGuess)
ROT3$Infomax$loadings

ROT4<-RotateSlopes(fit.data = CCI.Fit.4D.NoGuess)
ROT4$Infomax$loadings

ROT5<-RotateSlopes(fit.data = CCI.Fit.5D.NoGuess)
ROT5$Infomax$loadings

ROT6<-RotateSlopes(fit.data = CCI.Fit.6D.NoGuess)
ROT6$Infomax$loadings

ROT7<-RotateSlopes(fit.data = CCI.Fit.7D.NoGuess)
ROT7$Infomax$loadings

ROT8<-RotateSlopes(fit.data = CCI.Fit.8D.NoGuess)
round(ROT8$Infomax$loadings, digits=3)
round(ROT8$Varimax$loadings, digits=3)
round(ROT8$Bifactor$loadings, digits=3)

ROT9<-RotateSlopes(fit.data = CCI.Fit.9D.NoGuess)
round(ROT9$Infomax$loadings, digits=3)
round(ROT9$Varimax$loadings, digits=3)
round(ROT9$Bifactor$loadings, digits=3)

ROT10<-RotateSlopes(fit.data = CCI.Fit.10D.NoGuess)
round(ROT10$Infomax$loadings, digits=3)
round(ROT10$Varimax$loadings, digits=3)
round(ROT10$Bifactor$loadings, digits=3)

ROT11<-RotateSlopes(fit.data = CCI.Fit.11D.NoGuess)
ROT11$Infomax$loadings

ROT13<-RotateSlopes(fit.data = CCI.Fit.13D.NoGuess)
round(ROT13$Infomax$loadings, digits=3)

ROT14<-RotateSlopes(fit.data = CCI.Fit.14D.NoGuess)
round(ROT14$Infomax$loadings, digits=3)

apply(cbind(1:22,round(ROT10$Infomax$loadings, digits=3)),1,function(x) (paste(x,collapse=" & ")))


################# QOL

load("RealData/_A1.rda")
QOL.Fit.1D<-FitDATA
load("RealData/QOL")
QOL.Fit.2D<-FitDATA
load("RealData/CCI_NoGuess_A3.rda")
QOL.Fit.3D<-FitDATA
load("RealData/CCI_NoGuess_A4.rda")
QOL.Fit.4D<-FitDATA
load("RealData/CCI_NoGuess_A5.rda")
QOL.Fit.5D<-FitDATA
load("RealData/CCI_NoGuess_A10.rda")

par(mfrow=c(1,1),mar=c(5,4,3,2))
barplot(table((as.numeric(cut2(rowMeans(R.QOL-2,na.rm = TRUE),cuts=seq(from = -2,to = 2,by = 0.1)))-21)*(1/10)),
        main="Distribution of QOL Scores",ylab="Frequency",xlab="Average Response")


ResponseCurves(R.QOL-2,scores = rowMeans(R.QOL-2,na.rm = TRUE),
               prows = 6,pcols = 4,correct = rep(2,ncol(R.QOL)),j.legend = c(6,18),score.bins = 25)
               


MultipleTWFitTests(fit.data.list = list(QOL.Fit.1D,
                                        QOL.Fit.2D,
                                        QOL.Fit.3D,
                                        QOL.Fit.4D),title = "QOL")
mean(rowMeans(R.QOL-3,na.rm = TRUE))
sd(rowMeans(R.QOL-3,na.rm = TRUE))

QOL.TEST<-c("I wished I had more friends.",
"I want to spend more time with my family.",
"I was good at talking with adults.",
"I was good at making friends.",
"I have trouble getting along with other kids my age.",
"I had trouble getting along with my family.",
"Other kids were mean to me.",
"I got into a yelling fight with other kids.",
"I felt accepted by other kids my age.",
"I felt loved by my parents or guardians.", 
"I was afraid of other kids my age.",
"Other kids made fun of me.",
"I felt different from other kids my age.",
"I got along better with adults than other kids my age.",
"I felt comfortable with other kids my age.",
"My teachers understood me.",
"I had problems getting along with my parents or guardians.",
"I felt good about how I got along with my classmates.",
"I felt bad about how I got along with my friends.",
"I felt nervous when I was with other kids my age.", 
"I did not want to be with other kids.",
"I could talk with my friends.",
"Other kids wanted to be with me.", 
"I did things with other kids my age.")


QOL.ROT3<-RotateSlopes(QOL.Fit.3D)
for (i in 1:nrow(QOL.ROT3$Oblimin$loadings)) {
  cat(i,"&",paste0(c(round(QOL.ROT3$Oblimin$loadings[i,],digits=2),QOL.TEST[i]),collapse=" & "),"\\\\ \\hline\n")
}

QOL.Study<-list()
S<-QOL.Fit.3D$EZZ
S.gEV<-eigen(S, symmetric=TRUE)
S.cEV<-eigen(cov2cor(S), symmetric=TRUE)
gEV.sigma<-c()
cEV.sigma<-c()
Bi<-400+c(110,120,130,140,150,160,170,180,190,200)
for (b in Bi) {
  settings$burnin<-b
  settings$empiricalse<-FALSE
  settings$esttheta<-FALSE
  settings$Adim<-3; settings$tmu<-rep(0,settings$Adim); settings$tsigma<-diag(settings$Adim)
  settings$estfile<-paste0(data.dir,"/QOL_",paste0("Burn",b),"_A",settings$Adim)
  QOL.Fit.3D.B<-AnalyzeTestData(RP = Response,settings = settings,verbose = TRUE,TargetA = NA)
  QOL.Fit.3D.B$S.gEV<-eigen(QOL.Fit.3D.B$EZZ, symmetric=TRUE)
  QOL.Fit.3D.B$S.cEV<-eigen(cov2cor(QOL.Fit.3D.B$EZZ), symmetric=TRUE)
  QOL.Fit.3D.B$ROT3<-RotateSlopes(fit.data = QOL.Fit.3D.B)
  QOL.Study[[paste0("Burn",b)]]<-QOL.Fit.3D.B
  gEV.sigma<-c(gEV.sigma,sd(QOL.Fit.3D.B$S.gEV$values-S.gEV$values))
  cEV.sigma<-c(cEV.sigma,sd(QOL.Fit.3D.B$S.cEV$values-S.cEV$values))
}

# for (b in names(QOL.Study)) {
#   QOL.Study[[b]]$ROT3<-RotateSlopes(fit.data = QOL.Study[[b]])
# }

ROT3.Q<-RotateSlopes(fit.data = QOL.Fit.3D)
R.O.A.sigma<-c()
R.I.A.sigma<-c()
R.V.A.sigma<-c()
for (b in names(QOL.Study)) {
  QOL.Fit.3D.B<-QOL.Study[[b]]
  R.O.A.sigma<-c(R.O.A.sigma,sd(QOL.Fit.3D.B$ROT3$Oblimin$loadings-ROT3.Q$Oblimin$loadings))
  R.I.A.sigma<-c(R.I.A.sigma,sd(QOL.Fit.3D.B$ROT3$Infomax$loadings-ROT3.Q$Infomax$loadings))
  R.V.A.sigma<-c(R.V.A.sigma,sd(QOL.Fit.3D.B$ROT3$Varimax$loadings-ROT3.Q$Varimax$loadings))
}
par(mfrow=c(1,1))
plot(c(0,10*length(QOL.Study)),range(c(R.V.A.sigma,R.O.A.sigma,R.I.A.sigma)),
     type="n",main="Variation of Converged Rotated Loadings By Burn-In Iterations",
     xlab="Burnin Iterations",ylab=expression(sigma[A]),ylim=c(0,max(c(R.V.A.sigma,R.O.A.sigma,R.I.A.sigma))))
lines(10*1:length(QOL.Study),R.O.A.sigma,lty=2,col=2)
lines(10*1:length(QOL.Study),R.V.A.sigma,lty=3,col=3)
lines(10*1:length(QOL.Study),R.I.A.sigma,lty=4,col=4)
points(10*1:length(QOL.Study),R.O.A.sigma,pch=15,col=2)
points(10*1:length(QOL.Study),R.V.A.sigma,pch=16,col=3)
points(10*1:length(QOL.Study),R.I.A.sigma,pch=17,col=4)
legend("topright",c("Oblimin","Varimax","Infomax"),pch=15:17,col=2:4,lty=2:4)

#saveRDS(QOL.Study,"QOLFITS.rds")

# settings$burnin<-80
# settings$Adim<-5; settings$tmu<-rep(0,settings$Adim); settings$tsigma<-diag(settings$Adim)
# settings$estfile<-paste0(data.dir,"/QOL_Timed_A",settings$Adim)
# QOL.Fit.5D.T<-AnalyzeTestData(RP = Response,settings = settings,verbose = TRUE,TargetA = NA,timed = TRUE)

QOL.1D.Thresholds<-cbind(-1*((QOL.Fit.1D$tau)*1.702),QOL.Fit.1D$A*1.702)

X<-scan()  
  1 3.32 3.32 2.60 2.61 .82 .82 -.09 -.08 1.31 .10 1.31 .11
  2 3.07 3.07 2.19 2.19 .56 .56 -.59 -.59 .11 .07 .11 .07
  3 3.42 3.42 2.70 2.71 .77 .77 -.40 -.39 .73 .08 .73 .09
  4 4.53 4.53 3.72 3.72 1.88 1.89 .35 .36 1.62 .13 1.62 .13
  5 3.28 3.28 2.58 2.59 1.03 1.04 -.39 -.38 1.32 .10 1.31 .11
  6 3.42 3.42 2.46 2.46 .79 .79 -.41 -.40 .93 .09 .92 .09
  7 4.14 4.13 3.14 3.14 1.10 1.10 -.60 -.59 1.96 .13 1.94 .14
  8 3.45 3.45 2.66 2.66 1.13 1.13 .11 .12 .82 .09 .82 .09
  9 3.58 3.59 2.94 2.95 1.54 1.56 .17 .18 1.67 .12 1.67 .13
  10 4.04 4.05 3.77 3.78 2.74 2.74 1.62 1.63 .84 .12 .85 .12
  11 4.65 4.65 3.91 3.91 2.65 2.66 1.33 1.34 1.22 .13 1.23 .13
  12 4.22 4.23 3.30 3.30 1.70 1.70 .21 .21 1.71 .12 1.71 .13
  13 3.88 3.88 2.71 2.72 1.07 1.08 -.05 -.04 1.70 .12 1.70 .13
  14 2.52 2.52 1.60 1.60 -.28 -.27 -1.24 -1.24 .87 .08 .87 .09
  15 3.93 3.95 3.26 3.28 1.50 1.51 .19 .20 1.69 .13 1.70 .13
  16 3.61 3.61 2.91 2.92 1.20 1.21 .12 .13 .91 .09 .91 .09
  17 3.65 3.64 2.69 2.69 1.23 1.23 -.14 -.14 .71 .09 .70 .09
  18 4.46 4.47 3.85 3.86 1.77 1.78 .26 .27 1.88 .14 1.88 .14
  19 3.98 4.00 3.12 3.13 1.81 1.83 .47 .47 1.33 .11 1.34 .12
  20 4.72 4.72 3.61 3.62 1.74 1.75 .33 .34 1.48 .11 1.48 .12
  21 4.42 4.43 3.43 3.44 1.81 1.82 .55 .55 1.48 .12 1.49 .13
  22 4.64 4.65 4.03 4.05 2.78 2.80 1.38 1.40 1.62 .14 1.63 .15
  23 4.20 4.20 3.61 3.62 1.22 1.23 -.55 -.55 1.51 .11 1.51 .12
  24 4.17 4.17 3.20 3.21 1.60 1.61 .04 .05 1.32 .11 1.32 .11
  
LC<-matrix(X,nrow=24,ncol=13,byrow = TRUE)
LC<-LC[,c(2,4,6,8,10)]
LiCai.MH.1D<-LC

for (i in 1:ncol(Response)) {
  cat(i,"&",paste0(round(QOL.1D.Thresholds[i,],digits=2)," & ",LiCai.MH.1D[i,],collapse=" & "),"\\\\ \\hline \n")
}

Target.R<-cbind(rep(NA,24),
                c(0,0,NA,NA,0,0,0,0,NA,NA,0,0,0,0,NA,NA,0,NA,0,0,0,NA,NA,NA),
                c(0,NA,0,0,0,NA,0,0,0,NA,0,0,0,0,0,0,NA,0,0,0,0,0,0,0),
                c(0,0,NA,0,0,0,0,0,0,NA,0,0,0,NA,rep(0,10)),
                c(rep(c(0,NA),c(6,1)),rep(0,4),NA,rep(0,12)))
ROT5.T<-targetQ(QOL.Fit.5D$A,Target=Target.R)
ROT5.T$loadings<-matrix(rep(ifelse(colSums(ROT5.T$loadings)<0,-1,1),
                            nrow(QOL.Fit.5D$A)),ncol=ncol(QOL.Fit.5D$A),
                        nrow=nrow(QOL.Fit.5D$A),byrow=TRUE)*ROT5.T$loadings

QOL.Target<-round((ROT5.T$loadings)/sqrt(1+(ROT5.T$loadings)^2),digits = 2)

LiCai.MH<-rbind(c(.62 ,.14 ,-.12,.04 ,.02),
c(.00 ,.19 ,.40 ,.18 ,-.13),
c(.28,.38,.09,.32,.06),
c(.57,.47,.01,.16,-.06),
c(.61,.08,.02, -.01 ,.06),
c(.52, -.08,.55 , .02 ,.05),
c(.76,.05,-.09,.06 ,.55),
c(.44,-.04,.24,-.01,.27),
c(.55,.52,-.08,.07,.01),
c(.29,.44,.43,.34,-.03),
c(.73,-.16,-.17,.30,-.09),
c(.72,-.01,.02,.00,.51),
c(.75,-.01,.06,-.08,.13),
c(.50,.01,-.03,-.51,-.05),
c(.56,.54,-.06,.04,.01),
c(.34,.34,.28,-.01,.10),
c(.43,-.08,.68,-.05,.00),
c(.56,.56,.16,-.03,.09),
c(.65,.01,.19,-.19,-.04),
c(.75,-.03,-.03,-.01,-.15),
c(.69,.09,.04,-.15,-.08),
c(.53,.56,-.05,-.04,-.05),
c(.50,.56,-.12,-.04,.01),
c(.46,.57,-.13,-.11,-.07))      

for (i in 1:ncol(Response)) {
  cat(i,"&",paste0(ifelse(is.na(Target.R[i,]),"X",0)," & ",QOL.Target[i,]," & ",LiCai.MH[i,],collapse=" & ")," \\\\ \\hline \n")
}
