
dirName<-"ObliqueSims"
dir.create(dirName)
# A.   5 factors @ .4
# B.   3 factors @ .5
# C.   2 factors @ .6

r0=diag(10)

r1=diag(c(rep(0.6,5),rep(0.5,3),rep(0.4,2)))+
  matrix(c(rep(rep(c(0.4,0),c(5,5)),5),
           rep(rep(c(0,0.5,0),c(5,3,2)),3),
           rep(rep(c(0,0.6),c(8,2)),2)),10,10)

priorTheta = list(a=list(mu=rep(0,10),sig=r0,n=100000),
                  b=list(mu=rep(0.5,10),sig=r0,n=100000),
                  c=list(mu=rep(1,10),sig=r0,n=100000),
                  d=list(mu=rep(0,10),sig=r1,n=100000),
                  e=list(mu=rep(0.5,10),sig=r1,n=100000),
                  f=list(mu=rep(1,10),sig=r1,n=100000))
load("~/Downloads/Poly_J100_N1e+05_Q10_K4.rda")


gen.list<-list()
gen.list$gen.xi<-gen.xi
gen.list$gen.tau<-gen.tau
  # Name of saved file of generated data.
structure = CheckParams(parameters = gen.structure,generate = TRUE) # Explained in a hot minute.

for (genState in 1:length(priorTheta)) {
  # genState=1
  final.list<-priorTheta[[genState]]
  final.list$gen.xi<-gen.list$gen.xi
  final.list$gen.tau<-gen.list$gen.tau
  final.list$gen.theta<-GenMVNorm(priorTheta[[genState]]$n, 
                                  priorTheta[[genState]]$mu, 
                                  priorTheta[[genState]]$sig)
  structure=list(icc="ogive",          # Item Char Curve; "ogive" or "logistic" 
                 Adist="subscale",          # prior distribution of A's/loadings
                 Aparams=c(0.2,1.7),   # parameters of A's/loadings' prior distribution
                 Adim=ncol(final.list$gen.theta),               # 1 (univariate) or 2, 3, etc. multiple dimensions for multivariate
                 bdist="norm",         # distribution of B/intercept
                 bparams=c(0,1),       # parameters of bdist
                 guess=FALSE,          # guessing ? TRUE/FALSE
                 ncat=ncol(final.list$gen.tau)+1,               # Ordinal Polythomous? Number of categories
                 taudist="norm",       # sample distribution for categorical intercepts
                 tauparams=c(0,1),     # parameters for taudist
                 tmu=priorTheta[[genState]]$mu,         # Theta Prior... e.g. 0, or multivariate c(0,0) ... can be multidimensional
                 tsigma=priorTheta[[genState]]$sig, # Latent factor orthogonal covariance
                 simfile=paste0("Condition_",names(priorTheta)[genState],"_generatedData"))
  structure<-CheckParams(parameters = structure,generate = TRUE)
  final.list$structure<-structure
  final.list$gen.rp<-GenerateRP(xi = final.list$gen.xi,
                                theta = final.list$gen.theta,
                                tau = final.list$gen.tau,
                                structure = final.list$structure)
  saveRDS(final.list,paste0(dirName,"/",names(priorTheta)[genState],".rds"))
}


r0=diag(5)

r2=diag(5)+
  matrix(c(rep(c(0,0.5),c(1,4)),
           rep(rep(c(0.5,0),c(1,4)),4)),5,5)

priorTheta = list(g=list(mu=rep(0,5),sig=r0,n=10000),
                  h=list(mu=rep(0.5,5),sig=r0,n=10000),
                  i=list(mu=rep(1,5),sig=r0,n=10000),
                  j=list(mu=rep(0,5),sig=r2,n=10000),
                  k=list(mu=rep(0.5,5),sig=r2,n=10000),
                  l=list(mu=rep(1,5),sig=r2,n=10000))
load("~/Downloads/Poly_J100_N10000_Q5_K5.rda")                  

gen.list<-list()
gen.list$gen.xi<-gen.xi
gen.list$gen.tau<-gen.tau
# Name of saved file of generated data.
structure = CheckParams(parameters = gen.structure,generate = TRUE) # Explained in a hot minute.

for (genState in 1:length(priorTheta)) {
  # genState=1
  final.list<-priorTheta[[genState]]
  final.list$gen.xi<-gen.list$gen.xi
  final.list$gen.tau<-gen.list$gen.tau
  final.list$gen.theta<-GenMVNorm(priorTheta[[genState]]$n, 
                                  priorTheta[[genState]]$mu, 
                                  priorTheta[[genState]]$sig)
  structure=list(icc="ogive",          # Item Char Curve; "ogive" or "logistic" 
                 Adist="subscale",          # prior distribution of A's/loadings
                 Aparams=c(0.2,1.7),   # parameters of A's/loadings' prior distribution
                 Adim=ncol(final.list$gen.theta),               # 1 (univariate) or 2, 3, etc. multiple dimensions for multivariate
                 bdist="norm",         # distribution of B/intercept
                 bparams=c(0,1),       # parameters of bdist
                 guess=FALSE,          # guessing ? TRUE/FALSE
                 ncat=ncol(final.list$gen.tau)+1,               # Ordinal Polythomous? Number of categories
                 taudist="norm",       # sample distribution for categorical intercepts
                 tauparams=c(0,1),     # parameters for taudist
                 tmu=priorTheta[[genState]]$mu,         # Theta Prior... e.g. 0, or multivariate c(0,0) ... can be multidimensional
                 tsigma=priorTheta[[genState]]$sig, # Latent factor orthogonal covariance
                 simfile=paste0("Condition_",names(priorTheta)[genState],"_generatedData"))
  structure<-CheckParams(parameters = structure,generate = TRUE)
  final.list$structure<-structure
  final.list$gen.rp<-GenerateRP(xi = final.list$gen.xi,
                                theta = final.list$gen.theta,
                                tau = final.list$gen.tau,
                                structure = final.list$structure)
  saveRDS(final.list,paste0(dirName,"/",names(priorTheta)[genState],".rds"))
}

sim.list<-readRDS("ObliqueSims/g.rds")

sim.list$gen.xi - gen.list$gen.xi
settings<-CheckParams(generate = FALSE)
settings$empiricalse<-FALSE
settings$Adim<-ncol(sim.list$gen.theta)
settings$tsigma<-diag(ncol(sim.list$gen.theta))
settings$tmu<-rep(0,ncol(sim.list$gen.theta))
settings$estfile<-"ObliqueSims/G_Fit"
settings$ncat<-ncol(sim.list$gen.tau)+1
settings$thetamap<-FALSE
FIT.DATA<-AnalyzeTestData(sim.list$gen.rp,
                          settings = settings,
                          TargetA = sim.list$gen.xi[,1:ncol(sim.list$gen.theta)])

par(mfrow=c(3,2))
flip<-c(-1,-1,1,-1,1,1)
titles<-c(paste("Factor",1:5),"Intercept")
for (i in 1:ncol(sim.list$gen.xi)) plot(sim.list$gen.xi[,i],flip[i]*FIT.DATA$xi[,i],main=titles[i],
                                        ylab="Estimated",xlab="Simulated")
par(mfrow=c(2,2))
for (i in 1:ncol(sim.list$gen.tau)) plot(sim.list$gen.tau[,i],FIT.DATA$tau[,i],main=paste("Tau",i),
                                         ylab="Estimated",xlab="Simulated")

