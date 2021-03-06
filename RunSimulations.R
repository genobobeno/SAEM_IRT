if (!"knitr" %in% installed.packages()[,1]) { install.packages("knitr",dependencies = TRUE)} 
knitr::opts_chunk$set(echo=FALSE,warning = FALSE)
source("Hello.R", chdir = T)
Hello()
PACK = c(ifelse(get_os()=="windows","snow","parallel"),"foreach",
         ifelse(get_os()=="windows","doSNOW","doParallel"),
         "rlecuyer","GPArotation","mvnfast","psych","modeest","MASS","stringr",
         "Rcpp","RcppArmadillo","devtools","compiler","doRNG","abind",
         "truncdist","truncnorm","modeest","combinat","fastGHQuad","bdsmatrix","mcmc")

### Install packages not already installed ###
for (i in 1:length(PACK)){
  if(PACK[i] %in% rownames(installed.packages()) == FALSE) {
    install.packages(PACK[i],dependencies = TRUE)}
}
required<-lapply(PACK, require, character.only = TRUE)                       # Call the libraries (will deprecate)
Rcpp::sourceCpp("MVNormRand.cpp")                                  # C++ MV Norm Random (fast)
file.sources = list.files(path = "./GeisCamilli/R/",pattern="*.R") # Grab the R functions!
sourced<-sapply(paste0("./GeisCamilli/R/",file.sources),source,.GlobalEnv)  # Source them!
source("CreateSimulationStructure.R")
set.seed(789)
#CheckFiles
for (d in names(sim.list)) {
  print(paste("Checking files for:",d))
  ## Check Generated files
  fs<-SFileString(sim.list[[d]],gen=TRUE)
  simdir<-paste0(gen.dir,"/",d)
  files<-dir(simdir)
  if (length(files)==0) {
    nf<-1:sim.list[[d]]$Reps
  } else {
    nf<-str_replace(files,paste0(fs,"_"),"")
    nf<-as.integer(str_replace(nf,"\\.rda",""))
    nf<-c(1:sim.list[[d]]$Reps)[!1:sim.list[[d]]$Reps %in% nf]
  }
  print(paste("  Files to generate:",length(nf)))
  if (exists("Gen.Data")) rm(Gen.Data)
  if (length(nf)>0) {
    gf<-paste0(fs,"_",nf)
    for (sim in gf) {
      print(paste("    File:",sim))
      mType = sim.list[[d]]$mType
      J = sim.list[[d]]$J; N <- sim.list[[d]]$N; Q <- sim.list[[d]]$Q; K <-sim.list[[d]]$K
      structure=list(icc="ogive",          # Item Char Curve; "ogive" or "logistic" 
                     Adist=mType,          # prior distribution of A's/loadings
                     Aparams=c(0.2,1.7),   # parameters of A's/loadings' prior distribution
                     Adim=Q,               # 1 (univariate) or 2, 3, etc. multiple dimensions for multivariate
                     bdist="norm",         # distribution of B/intercept
                     bparams=c(0,1),       # parameters of bdist
                     guess=sim.list[[d]]$Guessing,          # guessing ? TRUE/FALSE
                     ncat=K,               # Ordinal Polythomous? Number of categories
                     taudist="norm",       # sample distribution for categorical intercepts
                     tauparams=c(0,1),     # parameters for taudist
                     cdist="unif",         # guessing parameter distribution for 3PNO or 3PL
                     cparams=c(0.05,0.3),  # bounds
                     tmu=rep(0,Q),         # Theta Prior... e.g. 0, or multivariate c(0,0) ... can be multidimensional
                     tsigma=if(Q==1) 1 else diag(Q), # Latent factor orthogonal covariance
                     simfile=paste0(simdir,"/",sim))  # Name of saved file of generated data.
      structure = CheckParams(parameters = structure,generate = TRUE) # Explained in a hot minute.
      if (exists("Gen.Data")) {
        gen.list<-Gen.Data$gen.list
      } else if (any(!1:sim.list[[d]]$Reps %in% nf)) {
        gsim<-paste0(fs,"_",c(1:(sim.list[[d]]$Reps))[!1:sim.list[[d]]$Reps %in% nf][1])
        Gen.Data<-readRDS(paste0(simdir,"/",gsim,".rds"))
        gen.list<-Gen.Data$gen.list
      } else {gen.list<-NA}
      Gen.Data<-GenerateTestData(j=J,n=N,structure=structure,gen.list=gen.list)
    }
  }
}
  # for (r in 1:sim.list[[d]]$Reps) {
#   if (grepl(,dir(paste0(gen.dir,"/",d))))
#     
#CheckSims
#BeginSimulation
