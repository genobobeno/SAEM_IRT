#setwd("/home/egeis/Documents/Documents/Rutgers_GSE/Camilli/ParSAEM/SAEM_IRT")
J=40;N=20000;Q=10;K=4
structure=list(icc="ogive",          # Item Char Curve; "ogive" or "logistic" 
               Adist="test",         # prior distribution of A's/loadings
               Aparams=c(0.2,1.7),   # parameters of A's/loadings' prior distribution
               Adim=Q,               # 1 (univariate) or 2, 3, etc. multiple dimensions for multivariate
               bdist="norm",         # distribution of B/intercept
               bparams=c(0,1),    # parameters of B
               guess=FALSE,          # guessing ? TRUE/FALSE
               ncat=K,
               taudist="norm",
               tauparams=c(0,1),
               cdist="unif",         # guessing parameter distribution for 3PNO or 3PL
               cparams=c(0.05,0.3),  # bounds
               tmu=rep(0,Q),                # Theta Prior... e.g. 0, or multivariate c(0,0) # can be multidimensional
               tsigma=if(Q==1) 1 else diag(Q),
               simfile=paste0("Poly_J",J,"_N",N,"_Q",Q,"_K",K,".rda")) 
# source('GeisCamilli/R/GenerateA.R')
# source('GeisCamilli/R/GenerateB.R')
# source('GeisCamilli/R/GenerateC.R')
# source('GeisCamilli/R/GenerateTau.R')
# source('GeisCamilli/R/GenerateTheta.R')
# source('GeisCamilli/R/GenMVNorm.R')
# source('GeisCamilli/R/ProbOgive.R')
# source('GeisCamilli/R/ProbIRT.R')
# source('GeisCamilli/R/GenerateRP.R')
# source('GeisCamilli/R/GenerateTestData.R')

GEN.DATA = GenerateTestData(j=J,n=N,structure = structure)
#saveRDS(GEN.DATA,paste0("Poly_J",J,"_N",N,"_Q",Q,"_K",K,".rds"))
#GEN.DATA$ResponseData<-data.frame(cbind(GEN.DATA$THETA,GEN.DATA$RP))
#names(GEN.DATA$ResponseData)<-c(paste0("Theta",1:Q),paste0("J",1:J))
#write.csv(GEN.DATA$ResponseData,file = paste0("Poly_Q",Q,"_K",K,"_N",N,"_ResponseData.csv"),row.names=FALSE)

settings = list(Adim=Q,guess=FALSE,empiricalse=FALSE,est="rm",estgain=1,burnin=100,ncat=K,
                plots=TRUE,plotiter=5,tmu=rep(0,Q),tsigma=diag(Q),eps=1e-4,
                parallel=(parallel::detectCores()>2),thetamap=FALSE,
                estfile="FitFile_Test",record=TRUE)
settings<-CheckParams(parameters = settings,generate=FALSE)
MyFirstFit<-AnalyzeTestData(RP=GEN.DATA$RP,settings=settings,TargetA = GEN.DATA$XI[,1:Q])

load("D:/Documents/Rutgers/SAEM_IRT/Poly_J100_N1e+05_Q10_K4.rda")
rp <- gen.rp
Q<-10
settings = list(Adim=Q,guess=FALSE,empiricalse=FALSE,est="rm",estgain=1,burnin=100,ncat=K,
                plots=TRUE,plotiter=5,tmu=rep(0,Q),tsigma=diag(Q),eps=1e-4,
                parallel=(parallel::detectCores()>2),thetamap=FALSE,
                estfile="FitFile_Test",record=TRUE)
settings<-CheckParams(parameters = settings,generate=FALSE)
MyFirstFit<-AnalyzeTestData(RP=rp,settings=settings,TargetA = gen.xi[,1:Q])

parallel::stopCluster(cl)
