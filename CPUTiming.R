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
source("CreateSimulationStructure.R")


## Check fit directories
fit.dir<-"CPUTests"
if (!fit.dir %in% dir()) dir.create(fit.dir)

##### START FITS
# for (r in 1:sim.list[[d]]$Reps) {
#   if (grepl(,dir(paste0(gen.dir,"/",d))))
#     
#CheckSims
  # for (d in names(sim.list)[1:3]) {
  #   print(paste("Checking files for:",d))
  #   ## Check Generated files
  #   fs<-SFileString(sim.list[[d]],gen=TRUE)
  #   simdir<-paste0(gen.dir,"/",d)
  #   files<-dir(simdir)
  #   files<-files[grep("\\.rda",files)]
  #   if (length(files)==0) {
  #     nf<-1:sim.list[[d]]$Reps
  #   } else {
  #     nf<-str_replace(files,paste0(fs,"_"),"")
  #     nf<-as.integer(str_replace(nf,"\\.rda",""))
  #     if (sum(1:sim.list[[d]]$Reps %in% nf)==sim.list[[d]]$Reps) {
  #       print(paste0(d,": All Simulations Accounted for."))
  #     } else {
  #       print(paste0(d,": Simulations not adding up."))
  #       nf<-c(1:sim.list[[d]]$Reps)[!1:sim.list[[d]]$Reps %in% nf]
  #       print(paste(length(nf),"Files to generate. Files missing:"))
  #       print(paste(nf,collapse=","))
  #     }
  #   }
  # }

  #BeginFits
for (d in names(sim.list)[1:3]) {
  filename<-paste0(d,".rds")
  if (!filename %in% dir(fit.dir)) {
    CPU.d<-list()
    CPU.r<-data.frame(cores=c(1,2,4,6),avg.init=NA,avg.saeminit=NA,avg.time=NA,avg.iter=NA,avg.tpi=NA)
    print(paste("Sim",d))
    for (cpu in c(1,2,4,6)) {
      print(paste("CPU",cpu))
      sink(paste0(fit.dir,"/",d,"_",cpu,"_output.txt"))
      avg.init<-c()
      avg.saeminit<-c()
      avg.time<-c()
      avg.iter<-c()
      avg.tpi<-c()
      cat(paste0("Starting Sim: ",d,"\n"))
      simdir<-paste0(gen.dir,"/",d)
      gfiles<-dir(simdir)
      gfiles<-gfiles[grep("\\.rds",gfiles)]
      R<-sort(as.numeric(str_extract_all(gfiles,paste0("(?<=_R",sim.list[[d]]$Reps,"_)[0-9]+"))))
      if (max(R)>5) R<-1:10 else R<-1:5
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
                        empiricalse=FALSE, # Get empirical SEs by restarting sampling at converged estimates. TRUE or FALSE
                        thinA=7, # Get empirical SEs by restarting sampling at converged estimates. TRUE or FALSE
                        thinB=5, # Get empirical SEs by restarting sampling at converged estimates. TRUE or FALSE
                        EmpIT=2000, # Iterations of Empirical Errors
                        est="rm",        # Estimation method? model("gifa) -> "off"=mean, "rm"=robbinsmonro, "sa"=simannealing; model("irt") -> "anr"=analytical newton-raphson, "nnr"=numerical newton-raphson #Convergence procedure
                        estgain=1,       # Constant to slow down(decrease to decimal)/speed up(increase) newton cycles, or exponent on denominator in rm-squeeze
                        burnin=as.integer(5000000/sim.list[[d]]$N),      # MCMC Burn-In iterations... or some other convergence criteria, acf? or Euclidean distance?
                        quad="manual",   # gauss-hermite
                        nodes=15,        # nodes for quadrature
                        gridbounds=c(-3.5,3.5), #manual quadrature bounds
                        tmu=rep(0,sim.list[[d]]$Q),
                        tsigma=diag(sim.list[[d]]$Q),        # Prior sigma, can be multivariate matrix(Adim x Adim)
                        eps=1e-4,        # Converged?
                        nesttheta=NA,    # if esttheta="mcmc", how many random samples?
                        thetamap=FALSE,   # do an MAP estimate of Theta?
                        thetaGen=NA, #genlist$gen.theta, # Did you simulate a new set of thetas? if so, give them to me. 
                        impute=FALSE,    # Impute missing data?
                        plots=FALSE,     # Show plots for diagnostics
                        chains=5,        # How many chains to build? Diagnose convergence? Simultaneous MCMC chains?
                        initialize="best", # "best", "random"
                        record=FALSE,     # "off"
                        parallel=TRUE,  # True or false for parallel computation?
                        cores=cpu,
                        esttheta=FALSE,
                        simfile=NA, #paste0(gfile,".rda"), # or NA
                        estfile=NA) 
          settings<-CheckParams(parameters = settings,generate = FALSE)
          Fit1D<-AnalyzeTestData(RP=genlist$gen.rp,settings=settings,timed=TRUE)
          avg.init<-c(avg.init,Fit1D$settings$timed.initialize)
          avg.saeminit<-c(avg.saeminit,Fit1D$settings$timed.SAEM_Init)
          avg.time<-c(avg.time,Fit1D$settings$timed.SAEM_Cycles)
          avg.iter<-c(avg.iter,Fit1D$settings$timed.Iterations)
          avg.tpi<-c(avg.tpi,Fit1D$settings$timed.SAEM_Cycles/Fit1D$settings$timed.Iterations)
        } else {
          settings = list(Adim=sim.list[[d]]$Q,
                          guess=FALSE,
                          est="rm",
                          estgain=1,
                          empiricalse=FALSE, # Get empirical SEs by restarting sampling at converged estimates. TRUE or FALSE
                          thinA=10, # Get empirical SEs by restarting sampling at converged estimates. TRUE or FALSE
                          thinB=7, # Get empirical SEs by restarting sampling at converged estimates. TRUE or FALSE
                          EmpIT=2000, # Iterations of Empirical Errors
                          burnin=as.integer(5000000/sim.list[[d]]$N),
                          ncat=sim.list[[d]]$K,
                          plots=FALSE,
                          plotiter=10,
                          tmu=rep(0,sim.list[[d]]$Q),
                          tsigma=diag(sim.list[[d]]$Q),
                          eps=1e-4,
                          thetaGen=NA, #genlist$gen.theta, # Did you simulate a new set of thetas? if so, give them to me. 
                          esttheta=FALSE,
                          nesttheta=NA,    # if esttheta="mcmc", how many random samples?
                          parallel=TRUE,  # True or false for parallel computation?
                          cores=cpu,
                          simfile=NA, #paste0(gfile,".rda"), # or NA
                          estfile=NA,
                          thetamap=FALSE,
                          record=FALSE)
          settings<-CheckParams(parameters = settings,generate=FALSE)
          Fit2D<-AnalyzeTestData(RP=genlist$gen.rp,settings=settings,timed=TRUE) #,TargetA = genlist$gen.xi[,1:sim.list[[d]]$Q]) 
          avg.init<-c(avg.init,Fit2D$settings$timed.initialize)
          avg.saeminit<-c(avg.saeminit,Fit2D$settings$timed.SAEM_Init)
          avg.time<-c(avg.time,Fit2D$settings$timed.SAEM_Cycles)
          avg.iter<-c(avg.iter,Fit2D$settings$timed.Iterations)
          avg.tpi<-c(avg.tpi,Fit2D$settings$timed.SAEM_Cycles/Fit2D$settings$timed.Iterations)
        }
      }
      sink(NULL)
      CPU.r$avg.init[CPU.r$cores==cpu]<-mean(avg.init)
      CPU.r$avg.saeminit[CPU.r$cores==cpu]<-mean(avg.saeminit)
      CPU.r$avg.time[CPU.r$cores==cpu]<-mean(avg.time)
      CPU.r$avg.iter[CPU.r$cores==cpu]<-mean(avg.iter)
      CPU.r$avg.tpi[CPU.r$cores==cpu]<-mean(avg.tpi)
      CPU.d[[paste0(d,"_",cpu)]]<-list(avg.init=avg.init,
                                       avg.saeminit=avg.saeminit,
                                       avg.time=avg.time,
                                       avg.iter=avg.iter,
                                       avg.tpi=avg.tpi)
    }
    saveRDS(list(Averages=CPU.r,Data=CPU.d),paste0(fit.dir,"/",filename))
  } else {
    print(paste(filename,"already created."))
  }
}
# parallel::stopCluster(cl)
  