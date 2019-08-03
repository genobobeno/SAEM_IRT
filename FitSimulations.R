TEST = TRUE

##### INITIALIZE
if (!"knitr" %in% installed.packages()[,1]) { install.packages("knitr",dependencies = TRUE)} 
knitr::opts_chunk$set(echo=FALSE,warning = FALSE)
source("Hello.R", chdir = T)
Hello()
PACK = c(ifelse(get_os()=="windows","snow","parallel"),"foreach",
         ifelse(get_os()=="windows","doSNOW","doParallel"),
         "rlecuyer","GPArotation","mvnfast","psych","modeest","MASS",
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
library(stringr)
if (TEST) {
  source("TestSimulationStructure.R")
} else {
  source("CreateSimulationStructure.R")
}

## Check fit directories

if (!fit.dir %in% dir()) dir.create(fit.dir)
for (d in names(sim.list)) {
  if (length(dir(fit.dir))==0 | !d %in% dir(fit.dir)) { 
    dir.create(paste0(fit.dir,"/",d))
  } else {
    print(paste0(fit.dir,"/",d," already exists."))
  }
}
##### START FITS
# for (r in 1:sim.list[[d]]$Reps) {
#   if (grepl(,dir(paste0(gen.dir,"/",d))))
#     
#CheckSims
for (d in names(sim.list)) { 
  print(paste("Checking files for:",d))
  ## Check Generated files
  fs<-SFileString(sim.list[[d]],gen=TRUE)
  simdir<-paste0(gen.dir,"/",d)
  files<-dir(simdir)
  files<-files[grep("\\.rda",files)]
  if (length(files)==0) {
    nf<-1:sim.list[[d]]$Reps
  } else {
    nf<-str_replace(files,paste0(fs,"_"),"")
    nf<-as.integer(str_replace(nf,"\\.rda",""))
    if (sum(1:sim.list[[d]]$Reps %in% nf)==sim.list[[d]]$Reps) {
      print(paste0(d,": All Simulations Accounted for."))
    } else {
      print(paste0(d,": Simulations not adding up."))
      nf<-c(1:sim.list[[d]]$Reps)[!1:sim.list[[d]]$Reps %in% nf]
      print(paste(length(nf),"Files to generate. Files missing:"))
      print(paste(nf,collapse=","))
    }
  }
}

#BeginFits
for (d in names(sim.list)) {
  cat(paste0("Starting Sim: ",d,"\n"))
  simdir<-paste0(gen.dir,"/",d)
  fitdir<-paste0(fit.dir,"/",d)
  if (fit.dir %in% dir()) gfiles<-dir(simdir)
  gfiles<-gfiles[grep("\\.rds",gfiles)]
  R<-sort(as.numeric(str_extract_all(gfiles,paste0("(?<=_R",sim.list[[d]]$Reps,"_)[0-9]+"))))
  if (!TEST & length(dir(paste0(fit.dir,"/",d)))>0) {
    ffiles<-dir(paste0(fit.dir,"/",d))
    fR<-sort(as.numeric(unique(str_extract_all(ffiles,paste0("(?<=_R",sim.list[[d]]$Reps,"_)[0-9]+")))))
    R<-R[!R %in% fR]
  } else if (TEST) {
    R<-1
    ffiles<-dir(paste0(fit.dir,"/",d))
    fR<-sort(as.numeric(unique(str_extract_all(ffiles,paste0("(?<=_R",sim.list[[d]]$Reps,"_)[0-9]+")))))
    if (length(fR)>0 && fR==R) {
      next
    }
  } 
  fs<-SFileString(sim.list[[d]],gen=TRUE)
  cat(paste0("Gotta fit: ",paste(R,collapse=" "),"\n"))
  for (r in R) {
    gc()
    gfile<-paste0(simdir,"/",fs,"_",r)
    cat(paste0("\n\nFitting ",gfile,": ",r,"\n"))
    genlist<-readRDS(paste0(gfile,".rds"))
    if (is.na(sim.list[[d]]$K) | sim.list[[d]]$K<3) {
      settings=list(model="gifa",    # Or "irt" = Analytical EM model
                    icc="ogive",     # Probability model... logistic? ogive? 
                    Adim=sim.list[[d]]$Q,          # Multidimensional?
                    guess=sim.list[[d]]$Guessing,     # Geussing ? TRUE
                    fm="eigen",    # Factor analysis procedure: fa() methods=c("ml","minres","wls","gls") or "licai", "camilli", or "old"...etc?
                    rmethod="pstT",  # GPArotation method, currently "targetT" or "pstT"
                    empiricalse=TRUE, # Get empirical SEs by restarting sampling at converged estimates. TRUE or FALSE
                    thinA=8, # Get empirical SEs by restarting sampling at converged estimates. TRUE or FALSE
                    thinB=5, # Get empirical SEs by restarting sampling at converged estimates. TRUE or FALSE
                    EmpIT=200, # Iterations of Empirical Errors
                    est="rm",        # Estimation method? model("gifa) -> "off"=mean, "rm"=robbinsmonro, "sa"=simannealing; model("irt") -> "anr"=analytical newton-raphson, "nnr"=numerical newton-raphson #Convergence procedure
                    estgain=1,       # Constant to slow down(decrease to decimal)/speed up(increase) newton cycles, or exponent on denominator in rm-squeeze
                    burnin=5000,  #as.integer(5000000/sim.list[[d]]$N),      # MCMC Burn-In iterations... or some other convergence criteria, acf? or Euclidean distance?
                    quad="manual",   # gauss-hermite
                    nodes=15,        # nodes for quadrature
                    gridbounds=c(-3.5,3.5), #manual quadrature bounds
                    tmu=0,           # Prior mean, can be multivariate vector(Adim)
                    tsigma=1,        # Prior sigma, can be multivariate matrix(Adim x Adim)
                    eps=1e-4,        # Converged?
                    nesttheta=25,    # if esttheta="mcmc", how many random samples?
                    thetamap=FALSE,   # do an MAP estimate of Theta?
                    thetaGen=genlist$gen.theta, # Did you simulate a new set of thetas? if so, give them to me. 
                    impute=FALSE,    # Impute missing data?
                    plots=FALSE,     # Show plots for diagnostics
                    chains=5,        # How many chains to build? Diagnose convergence? Simultaneous MCMC chains?
                    initialize="best", # "best", "random"
                    record=TRUE,     # "off"
                    parallel=TRUE,  # True or false for parallel computation?
                    cores=4,
                    simfile=paste0(gfile,".rda"), # or NA
                    estfile=paste0(fitdir,"/",SFileString(sim.list[[d]],gen=FALSE,r = r))) 
      settings<-CheckParams(parameters = settings,generate = FALSE)
      # Rprof()
      Fit1D<-AnalyzeTestData(RP=genlist$gen.rp,settings=settings,timed = TRUE)
      # Rprof(NULL)
      # print(summaryRprof())
      saveRDS(Fit1D,paste0(fitdir,"/",SFileString(sim.list[[d]],gen=FALSE,r = r),".rds"))
    } else {
      settings = list(Adim=sim.list[[d]]$Q,
                      guess=FALSE,
                      empiricalse=TRUE,
                      est="rm",
                      estgain=1,
                      empiricalse=TRUE, # Get empirical SEs by restarting sampling at converged estimates. TRUE or FALSE
                      thinA=8, # Get empirical SEs by restarting sampling at converged estimates. TRUE or FALSE
                      thinB=5, # Get empirical SEs by restarting sampling at converged estimates. TRUE or FALSE
                      EmpIT=2500, # Iterations of Empirical Errors
                      burnin=2500, #as.integer(5000000/sim.list[[d]]$N),
                      ncat=sim.list[[d]]$K,
                      plots=FALSE,
                      plotiter=10,
                      tmu=rep(0,sim.list[[d]]$Q),
                      tsigma=diag(sim.list[[d]]$Q),
                      eps=1e-4,
                      thetaGen=genlist$gen.theta, # Did you simulate a new set of thetas? if so, give them to me. 
                      nesttheta=25,    # if esttheta="mcmc", how many random samples?
                      parallel=TRUE,  # True or false for parallel computation?
                      cores=4,
                      simfile=paste0(gfile,".rda"), # or NA
                      estfile=paste0(fitdir,"/",SFileString(sim.list[[d]],gen=FALSE,r = r)),
                      thetamap=FALSE,
                      record=TRUE)
      settings<-CheckParams(parameters = settings,generate=FALSE)
      Fit2D<-AnalyzeTestData(RP=genlist$gen.rp,settings=settings,TargetA = genlist$gen.xi[,1:sim.list[[d]]$Q],
                             timed=TRUE) 
      saveRDS(Fit2D,paste0(fitdir,"/",SFileString(sim.list[[d]],gen=FALSE,r = r),".rds"))
    }
  }
}
  