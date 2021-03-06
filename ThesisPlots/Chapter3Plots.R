
#3.1.2 Distributions
sim.list$S1
Q<-1;d<-1;mType<-"beta";K<-2
Adist=mType          # prior distribution of A's/loadings
Aparams=c(0.2,1.7)   # parameters of A's/loadings' prior distribution
Adim=Q               # 1 (univariate) or 2, 3, etc. multiple dimensions for multivariate
bdist="norm"         # distribution of B/intercept
bparams=c(0,1)       # parameters of bdist
guess=sim.list[[d]]$Guessing          # guessing ? TRUE/FALSE
ncat=K               # Ordinal Polythomous? Number of categories
taudist="norm"       # sample distribution for categorical intercepts
tauparams=c(0,1)     # parameters for taudist
cdist="unif"         # guessing parameter distribution for 3PNO or 3PL
cparams=c(0.05,0.3)  # bounds
tmu=rep(0,Q)         # Theta Prior... e.g. 0, or multivariate c(0,0) ... can be multidimensional
tsigma=if(Q==1) 1 else diag(Q)

par(mfrow=c(3,2))
j=100
#d1
A<-rbeta(j,2.5,3.0)
A<-(Aparams[2]-Aparams[1])*A+Aparams[1]
A<-cbind(A,(1.0-0.1)*rbeta(j,2.5,3.0)+0.1)
xx<-seq(-0.5,1.5,length.out = 2000)
yy<-dbeta(x = xx,shape1 = 2.5,shape2 = 3.0)
plot(density(A[,1]),main=expression("Discrimination :" ~ italic(Q[1])),
     xlab=expression(italic(a)),xlim=c(0,2),col=2,lty=3)
lines((Aparams[2]-Aparams[1])*xx+Aparams[1],1.0/(Aparams[2]-Aparams[1])*yy)
legend("topright",legend=c(expression(italic(Q[1])~"~"~ beta(2.5,3)),expression(beta(2.5,3))),
       col=c(2,1),lty=c(3,1))
#d2
plot(density(A[,2]),main=expression("Discrimination :" ~ italic(Q["2+"])),
     xlab=expression(italic(a)),xlim=c(0,2),col=4,lty=3)
lines((1-0.1)*xx+0.1,1.0/0.9*yy)
legend("topright",legend=c(expression(italic(Q["2+"])~"~"~ beta(2.5,3)),expression(italic(beta)(2.5,3))),
       col=c(4,1),lty=c(3,1))

#Difficulty
b<-rnorm(j)
xx<-seq(-3.5,3.5,length.out = 7000)
yy<-dnorm(x = xx)
plot(density(b),main=expression("Difficulty"),col=2,lty=3,xlab=expression(italic(b)))
lines(xx,yy)
legend("topright",legend=c(expression(italic(b)~"~"~ italic(N)(0,1)),expression(italic(N)(0,1))),
       col=c(4,1),lty=c(3,1))

#Poly Thresholds
ncat=4
TAU = lapply(rep(ncat-1,j),function(x) {
  y = rnorm(n = x,tauparams[1],tauparams[2])
  sort(y)
})
tau = as.matrix(do.call(rbind,TAU))
plot(density(tau[,2]),main=expression("Polytomous Thresholds :"~italic(K)~"= 4"),xlim=c(-3.5,3.5),
     col=3,lty=3,xlab=expression(italic(tau)))
lines(density(tau[,1]),col=2,lty=3)
lines(density(tau[,3]),col=4,lty=3)
legend("topright",c(expression(italic(tau)[1]),expression(italic(tau)[2]),expression(italic(tau)[3])),lty=3,col=c(2,3,4))

#Guessing
cc<-runif(j,0.05,0.3)
xx<-seq(0,0.35,length.out = 3500)
yy<-dunif(x = xx,min = 0.05,max = 0.3)
plot(density(cc),main=expression(italic(P)("Guessing")),col=2,lty=3,xlab=expression(italic(c)))
lines(xx,yy)
legend("topright",c(expression(italic(c)~"~"~ italic(U)(0.05,0.3)),expression(italic(U)(0.05,0.3))),lty=3,col=c(2,1))

#Ability
theta<-rnorm(5000)
xx<-seq(-3.5,3.5,length.out = 7000)
yy<-dnorm(x = xx)
plot(density(theta),main=expression("Ability"),col=2,lty=3,xlab=expression(italic(theta)))
lines(xx,yy)
legend("topright",legend=c(expression(italic(theta)~"~"~ italic(N)(0,1)),expression(italic(N)(0,1))),
       col=c(4,1),lty=c(3,1))


#3.1.4  Convergence

burnin=1000;eps=0.001;estgain=1
RMwindow<-ceiling(burnin*(0.2))
GC<-c(rep(1,burnin-RMwindow),
        runif(rep(1,ceiling(RMwindow/2)),
              min = 1-(1:(RMwindow/2)/RMwindow)*cos(1:(RMwindow/2))*cos(1:(RMwindow/2)), 
              max = 1.0),
        runif(rep(1,ceiling(RMwindow/2)),
              min = 0.5-(1:(RMwindow/2)/RMwindow)*cos(1:(RMwindow/2))*cos(1:(RMwindow/2)), 
              max = 1-(1:(RMwindow/2)/RMwindow)*cos(1:(RMwindow/2))*cos(1:(RMwindow/2))),
        1/(2:(ceiling(2/eps)))^estgain)
GC.u<-c(rep(1,burnin-RMwindow),
        rep(1,ceiling(RMwindow/2)),
        1-(1:(RMwindow/2)/RMwindow)*cos(1:(RMwindow/2))*cos(1:(RMwindow/2)),
        1/(2:(ceiling(2/eps)))^estgain)[1:(2*burnin)]
backTrack<-sum(c(burnin-RMwindow,ceiling(RMwindow/2),RMwindow/2))

GC.l<-c(rep(1,burnin-RMwindow),
        1-(1:(RMwindow/2)/RMwindow)*cos(1:(RMwindow/2))*cos(1:(RMwindow/2)),
        0.5-(1:(RMwindow/2)/RMwindow)*cos(1:(RMwindow/2))*cos(1:(RMwindow/2)),
        1/(2:(ceiling(2/eps)))^estgain)[1:(2*burnin)]

# plot(GC.u[700:1100])
# lines(GC.l[1100:700])
par(mfrow=c(2,1),mar=c(1,4,4,3))
plot(750:1050,c(0:300)/300,type="n",main="Psuedo Annealing Gain Space : Burn-In = 1000 iterations",
     xlab="iterations",ylab=expression(gamma[t]),xaxt="n")
polygon(x = c(750:1050,1050:750,1),c(GC.u[750:1050],GC.l[1050:750],1),col="lightblue")
#lines(750:1050,GC[750:1050],col=1)
abline(v=c(800,1000),col=2)
par(mar=c(5,4,1,3))
plot(750:1050,c(0:300)/300,type="n",main=NULL,
     xlab="iterations",ylab=expression(gamma[t]))
lines(750:1050,GC[750:1050],col=1)
abline(v=c(800,1000),col=2)




source("CreateSimulationStructure.R")
d = "S1"
r = 1
simdir<-paste0(gen.dir,"/",d)
SimList<-readRDS(paste0(simdir,"/",SFileString(sim.list[[d]],gen=TRUE),"_",r,".rds"))
fitdir<-paste0(fit.dir,"/",d)
FitList<-readRDS(paste0(fitdir,"/",SFileString(sim.list[[d]],gen=FALSE,r = r),".rds"))
load(paste0(fitdir,"/",SFileString(sim.list[[d]],gen=FALSE,r = r),".rda"))
J=12
par(mfrow=c(2,1),mar=c(5,5,4,5))
Xlim<-c(0,100*ceiling(dim(MCMCDATA$Biter)[2]/100))
Ylim<-c(0,1.05*max(MCMCDATA$Biter[J,]))
plot(1:dim(MCMCDATA$Biter)[2],MCMCDATA$Biter[J,],type="n",main=paste0("MCMC Chain: Gain Constant and Difficulty (item ",J,")"),
     xlim=Xlim,ylim=Ylim,xlab="iterations",ylab=expression(italic(b)))
lines(1:dim(MCMCDATA$Biter)[2],MCMCDATA$Biter[J,],lty=1,col=2)
axis(4,at = c(0,0.25,0.5)*max(MCMCDATA$Biter[J,]),labels = c(0,0.5,1.0))
lines(1:dim(MCMCDATA$Biter)[2],FitList$gain[1:dim(MCMCDATA$Biter)[2]])
mtext(text = "RM gain constant",side = 4,padj = 3)
abline(h=SimList$gen.xi[J,2],lty=2,col=4)
text(70,2.2,expression(italic(hat(b))[12]),col=2)
text(1300,1.55,expression(italic(b)[12]),col=4)
text(770,0.7,expression(italic(gamma)[RM]),col=1)

zoomX<-c(750,1050)
Xlim<-zoomX
Ylim<-c(1.1,1.95)
plot(zoomX[1]:zoomX[2],MCMCDATA$Biter[J,zoomX[1]:zoomX[2]],type="n",main=paste0("Zoom In: Gain Constant and Difficulty (item ",J,")"),
     xlim=Xlim,ylim=Ylim,xlab="iterations",ylab=expression(italic(b)))
lines(zoomX[1]:zoomX[2],MCMCDATA$Biter[J,zoomX[1]:zoomX[2]],lty=1,col=2)
axis(4,at = (1.1+c(0,0.5,1.0)*0.4),labels = c(0,0.5,1.0))
lines(zoomX[1]:zoomX[2],1.1+0.4*FitList$gain[zoomX[1]:zoomX[2]])
mtext(text = "RM gain constant",side = 4,padj = 3)
abline(h=SimList$gen.xi[J,2],lty=2,col=4)
text(755,1.875,expression(italic(hat(b))[12]),col=2)
text(850,1.7,expression(italic(b)[12]),col=4)
text(875,1.25,expression(italic(gamma)[RM]),col=1)
######################


BenchmarkPlots(condition = "S1",all.reps = TRUE,basedir = getwd())

RMProofPlots(condition = "S1",items=c(5,11,37,38,53,63),basedir = getwd())

ChainPlots(condition = "S1",items=c(38,37,53,63,5,11),repl=NA,bins=200,basedir = getwd(),brewpalette = "Set2",cex=1)


source("CreateSimulationStructure.R")
d = "S1";r = 3;simdir<-paste0(gen.dir,"/",d);fitdir<-paste0(fit.dir,"/",d)
SimList<-readRDS(paste0(simdir,"/",SFileString(sim.list[[d]],gen=TRUE),"_",r,".rds"))
FitList<-readRDS(paste0(fitdir,"/",SFileString(sim.list[[d]],gen=FALSE,r = r),".rds"))
load(paste0(fitdir,"/",SFileString(sim.list[[d]],gen=FALSE,r = r),".rda"))
par(mfrow=c(4,2))
AutocovMC(MCMCDATA,start = 600,end = 800,items = 1,title.append = "Burn-in : ")   # Chain during burn
AutocovMC(MCMCDATA,start = 800,end = 1000,items = 1,title.append = "Annealing : ")   # Chain during burn
AutocovMC(MCMCDATA,start = 600,end = 800,items = 2,title.append = "Burn-in : ")   # Chain during burn
AutocovMC(MCMCDATA,start = 800,end = 1000,items = 2,title.append = "Annealing : ")   # Chain during burn

BenchmarkPlots(condition = "S2",all.reps = TRUE,basedir = getwd())



source("ThesisPlots/TWFitTests.R")

TWTestAndPlot("S1",extra.dimensions = TRUE,E = 7,ratios = FALSE)

par(mfrow=c(3,2),cex=0.9)
TWSplitSC("S1",extra.dimensions = TRUE,E=8)
TWSplitSC("S2",extra.dimensions = TRUE,E=8)
TWSplitSC("S3",extra.dimensions = TRUE,E=8)


BenchmarkPlots(condition = "S2",all.reps = TRUE,basedir = getwd(),cex=1)


BenchmarkPlots(condition = "S3",all.reps = TRUE,basedir = getwd(),cex=1)

CompareConditions("S1","S3",basedir = getwd())


#   # extra.dimensions=TRUE; d="S1";r=1
#   fd<-d
#   source("CreateSimulationStructure.R")
#   d<-fd
#   fit.dir.tw <- "TWFits"
#   twd <- paste0(d,c("Plus1","Plus2"))  #paste0(d,c("Plus1","Times5"))
#   fitdir<-c(paste0(fit.dir,"/",d),paste0(fit.dir.tw,"/",twd))
#   simdir<-paste0(gen.dir,"/",d)
#   SimList<-readRDS(paste0(simdir,"/",SFileString(sim.list[[d]],gen=TRUE),"_1.rds"))
#   FitList<-readRDS(paste0(fitdir[2],"/",SFileString(sim.list[[d]],gen=FALSE,r = 1),".rds"))
#   if (d %in% c("S1","S2","S3")) {
#     if (d=="S3") {
#       S<-FitList$EZZ #-t(as.matrix(FitList$EZ))%*%as.matrix(FitList$EZ)
#     } else {
#       S<-FitList$EZZ-as.matrix(FitList$EZ)%*%t(as.matrix(FitList$EZ))
#     }
#     gEV<-eigen(S, symmetric=TRUE)
#     cEV<-eigen(cov2cor(S),symmetric = TRUE)
#   } else {
#     gEV<-eigen(FitList$EZZ,symmetric=TRUE)
#     cEV<-eigen(cov2cor(FitList$EZZ),symmetric = TRUE)
#   }
# 
# psych::cortest.bartlett(R = cov2cor(S),n = sim.list[[d]]$N,diag = TRUE)
# psych::KMO(R = cov2cor(S))
# KaiserCriterion(J=sim.list[[d]]$J,N=sim.list[[d]]$N,S = S)
layout(matrix(c(1,4,2,5,3,6),nrow=3,ncol=2,byrow = TRUE))
FactorReconstruction(condition = "S4",all.reps = TRUE,basedir = getwd())
FactorReconstruction(condition = "S5",all.reps = TRUE,basedir = getwd())

par(mfrow=c(2,2),cex=1.1)
TWSplitSC("S4",extra.dimensions = TRUE,E = 7,
          basedir="C:/Dev/MyFiles/NextCloud/Documents/Rutgers/SAEM_IRT/",
          tw = 0.95)
TWSplitSC("S5",extra.dimensions = TRUE,E = 7,
          basedir="C:/Dev/MyFiles/NextCloud/Documents/Rutgers/SAEM_IRT/")

par(mfrow=c(1,2))
PlotEVRatios("S4",extra.dimensions = TRUE,E = 9)
PlotEVRatios("S5",extra.dimensions = TRUE,E = 9)
#TWTestAndPlot("S4",extra.dimensions = TRUE,E = 8,ratios = TRUE)
#TWTestAndPlot("S5",extra.dimensions = TRUE,E = 8,ratios = TRUE)

TWTestAndPlot(c("S4","S5"),extra.dimensions = TRUE,E=8,ratios=TRUE)

BenchmarkPlots(condition = "S6",all.reps = TRUE,basedir = getwd(),cex=0.75,mar=c(5,5,2,1))
BenchmarkPlots(condition = "S7",all.reps = TRUE,basedir = getwd(),cex=0.75,mar=c(5,5,2,1))

par(mfrow=c(2,2),cex=1)
TWSplitSC("S6",extra.dimensions = TRUE,E = 9)
TWSplitSC("S7",extra.dimensions = TRUE,E = 9)

TWTestAndPlot(c("S6","S7"),extra.dimensions = TRUE,E=8,ratios=TRUE)

BenchmarkPlots(condition = "S8",all.reps = TRUE,basedir = getwd(),cex=0.75,mar=c(5,5,2,1))
BenchmarkPlots(condition = "S9",all.reps = TRUE,basedir = getwd(),cex=0.75,mar=c(5,5,2,1))

par(mfrow=c(2,2),cex=1)
TWSplitSC("S8",extra.dimensions = TRUE,E = 15)
TWSplitSC("S9",extra.dimensions = TRUE,E = 15)

TWTestAndPlot(c("S8","S9"),extra.dimensions = TRUE,E=17,ratios=TRUE)

EVTestTab<-data.frame(Dimensions = c(1,1,1,3,3,5,5,10,10), LambdaN2 = c(NA,NA,NA,1,0,1,0,2,0), 
           LambdaN1 = c(NA,NA,NA,1,0,1,0,1,0), Lambda = c(0,0,0,0,0,0,0,1,0), 
           LambdaP1 = c(0,0,1,1,1,1,1,3,1), LambdaP2 = c(1,1,2,2,2,2,2,4,2))


