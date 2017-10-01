WD<-getwd()

##########################################################################
WorkingDirectory = "/home/egeis/Documents/Documents/Rutgers_GSE/Camilli/ParSAEM/"
setwd(WorkingDirectory)

PACK <- c("maptools","quantreg","MASS","plyr","ggplot2","Hmisc","fields",
          "truncdist","truncnorm","modeest","psych","GPArotation","abind",
          "combinat","fastGHQuad","bdsmatrix","doParallel","Rcpp","RcppArmadillo",
          "doRNG","devtools","compiler","mcmc","coda","pryr")

for (i in 1:length(PACK)){
  if(PACK[i] %in% rownames(installed.packages()) == FALSE) {install.packages(PACK[i],dependencies = TRUE)}
}

library("devtools")

if ("lineprof" %in% rownames(installed.packages()) == FALSE) devtools::install_github("hadley/lineprof")

library("maptools")  
library("quantreg")
library("MASS")
library("plyr")
library("ggplot2")
library("Hmisc")
library("fields")
library("truncdist")
library("truncnorm")
library("modeest")
library("psych")
library("GPArotation")
library("abind")
library("fastGHQuad")
library("combinat")
library("bdsmatrix")
library("doParallel")
library("Rcpp")
library("RcppArmadillo")
library("doRNG")
library("lineprof")
library("compiler")
library("mcmc")
library("coda")
library("pryr")

# improved list of objects
.ls.objects <- function (pos = 1, pattern, order.by,
                         decreasing=FALSE, head=FALSE, n=5) {
  napply <- function(names, fn) sapply(names, function(x)
    fn(get(x, pos = pos)))
  names <- ls(pos = pos, pattern = pattern)
  obj.class <- napply(names, function(x) as.character(class(x))[1])
  obj.mode <- napply(names, mode)
  obj.type <- ifelse(is.na(obj.class), obj.mode, obj.class)
  obj.prettysize <- napply(names, function(x) {
    capture.output(format(utils::object.size(x), units = "auto")) })
  obj.size <- napply(names, object.size)
  obj.dim <- t(napply(names, function(x)
    as.numeric(dim(x))[1:2]))
  vec <- is.na(obj.dim)[, 1] & (obj.type != "function")
  obj.dim[vec, 1] <- napply(names, length)[vec]
  out <- data.frame(obj.type, obj.size, obj.prettysize, obj.dim)
  names(out) <- c("Type", "Size", "PrettySize", "Rows", "Columns")
  if (!missing(order.by))
    out <- out[order(out[[order.by]], decreasing=decreasing), ]
  if (head)
    out <- head(out, n)
  out
}

# shorthand
lsos <- function(..., n=10) {
  .ls.objects(..., order.by="Size", decreasing=TRUE, head=TRUE, n=n)
}


enableJIT(1)

if (PARALLEL) {
  Sys.setenv("PKG_CXXFLAGS"="-fopenmp")
  Sys.setenv("PKG_LIBS"="-fopenmp")

  # Find out how many cores are available (if you don't already know)
  nCores<-detectCores()
  ## [1] 4
  # Create cluster with desired number of cores
  if (!exists('clust')) {
    clust <- makeCluster(4) #ceiling(nCores/2))
    # Register cluster
    registerDoParallel(clust)
    print("Made and registered cluster 'clust' as a parallel backend")
    # Find out how many are assigned
    print("NOTE: Don't forget to run 'stopCluster(clust)' upon finishing the program.")
    print(paste("You have",nCores,"Cores that can be used, we're running",getDoParWorkers()))
  } else {
    print("Cluster of cores 'cl' already created")
  }
}
  #Rcpp::sourceCpp("truncNorm.cpp")
Rcpp::sourceCpp("MVNormRand.cpp")

source("GeneratingFunctions.R")
source("Calculations.R")
source("FileManagement.R")
source("CodeDiagnostics.R")
source("Estimation.R")
source("ErrorFunctions.R")
source("GIFAFunctions.R")
source("IRTFunctions.R")
source("PlotFunctions.R")
source("PlotsOn.R")
source("EmpiricalSE.R")

# source("IRT_Response.R")
# source("RunIRT.R")
settings=list(model="gifa",    # Or "irt" = Analytical EM model
              icc="ogive",     # Probability model... logistic? ogive? 
              Adim=2,          # Multidimensional?
              Akappa=10,          # Multidimensional?
              guess=FALSE,     # Geussing ? TRUE
              fm="camilli",    # Factor analysis procedure: fa() methods=c("ml","minres","wls","gls") or "camilli", or "old", "pca",...etc?
              rmethod="pstT",    # GPArotation method, currently "targetT" or "pstT"
              empiricalse=TRUE, # Get empirical SEs by restarting sampling at converged estimates. TRUE or FALSE
              est="rm",       # Estimation method? model("gifa) -> "off"=mean, "rm"=robbinsmonro, "sa"=simannealing; model("irt") -> "anr"=analytical newton-raphson, "nnr"=numerical newton-raphson #Convergence procedure
              estgain=0.5,       # Constant to slow down(decrease to decimal)/speed up(increase) newton cycles, or exponent on denominator in rm-squeeze
              burnin=500,       # MCMC Burn-In iterations... or some other convergence criteria, acf? or Euclidean distance?
              quad="manual",   # manual, gauss-hermite
              nodes=15,        # nodes for quadrature
              gridbounds=c(-3.5,3.5), #manual quadrature bounds
              tmu=c(0,0),           # Prior mean, can be multivariate vector(Adim)
              tsigma=diag(2),        # Prior sigma, can be multivariate matrix(Adim x Adim)
              eps=1e-4,        # Converged?
              esttheta="mcmc", # "map", Must implement "eap"
              nesttheta=10,    # if esttheta="mcmc", how many random samples?
              impute=FALSE,    # Impute missing data?
              chains=1,        # How many chains to build? Diagnose convergence? Simultaneous MCMC chains?
              initialize="random", # "best", "random"
              plots=FALSE,     # Show plots for diagnostics
              record="on",     # "off"
              simfile=NA, # or NA
              estfile="GIFA_Results") 

#print(settings)
print("Use settings$estfile to declare a file to write results.")
print("If comparing to a simulation, use settings$simfile read A, b, theta.")
print("GIFA Simulation added, run using AnalyzeTestData(RP=RespPatterns,settings=settings)")
setwd(WD)
