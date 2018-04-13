############ RUN for Presentation
MyFirstFit<-AnalyzeTestData(RP=GEN.DATA$RP,settings=settings,TargetA = GEN.DATA$XI[,1:Q])

#parallel::stopCluster(cl)

#################################
load("D:/Documents/Rutgers/SAEM_IRT/Poly_J100_N1e+05_Q10_K4.rda")
rp <- gen.rp
Q<-10
settings = list(Adim=Q,guess=FALSE,empiricalse=FALSE,est="rm",estgain=1,burnin=100,ncat=K,
                plots=TRUE,plotiter=5,tmu=rep(0,Q),tsigma=diag(Q),eps=1e-4,
                parallel=(parallel::detectCores()>2),thetamap=FALSE,
                estfile="FitFile_Test",record=TRUE)
settings<-CheckParams(parameters = settings,generate=FALSE)
MyFirstFit<-AnalyzeTestData(RP=rp,settings=settings,TargetA = gen.xi[,1:Q])