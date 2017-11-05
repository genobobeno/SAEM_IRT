CheckParams = function(parameters,generate = TRUE) {
  parCheck = function(name,plist,dplist) {
    if (!tolower(name) %in% tolower(names(plist))) {
      plist[[name]] <- dplist[[name]]
    } 
    names(plist)[names(plist)==name] <- names(dplist)[which(tolower(names(dplist)) %in% tolower(name))]
    plist
  }
  if (generate) {
    defaultP=list( icc="ogive",          # Item Char Curve; "ogive" or "logistic" 
                   Adist="beta",         # prior distribution of A's/loadings
                   Aparams=c(0.2,1.7),   # parameters of A's/loadings' prior distribution
                   Adim=2,               # 1 (univariate) or 2, 3, etc. multiple dimensions for multivariate
                   bdist="norm",         # distribution of B/intercept
                   bparams=c(0,1),       # parameters of bdist
                   guess=FALSE,          # guessing ? TRUE/FALSE
                   ncat=NA,              # Ordinal Polythomous? Number of categories
                   taudist="norm",       # sample distribution for categorical intercepts
                   tauparams=c(0,1),     # parameters for taudist
                   cdist="unif",         # guessing parameter distribution for 3PNO or 3PL
                   cparams=c(0.05,0.3),  # bounds
                   tmu=rep(0,2),         # Theta Prior... e.g. 0, or multivariate c(0,0) ... can be multidimensional
                   tsigma=diag(2), # Latent factor orthogonal covariance
                   simfile=paste0("SimFile_",sample(1:99999,1),".rda"))  # Name of saved file of generated data.
  } else {
    defaultP=list(model="gifa",    # Or "irt" = Analytical EM model
                icc="ogive",     # Probability model... logistic? ogive? 
                Adim=2,          # Multidimensional?
                guess=FALSE,     # Geussing ? TRUE
                fm="eigen",    # Factor analysis procedure: fa() methods=c("ml","minres","wls","gls") or "licai", "camilli", or "old"...etc?
                rmethod="pstT",  # GPArotation method, currently "targetT" or "pstT"
                empiricalse=TRUE, # Get empirical SEs by restarting sampling at converged estimates. TRUE or FALSE
                thinA=5, # Get empirical SEs by restarting sampling at converged estimates. TRUE or FALSE
                thinB=2, # Get empirical SEs by restarting sampling at converged estimates. TRUE or FALSE
                EmpIT=1000, # Iterations of Empirical Errors
                est="rm",        # Estimation method? model("gifa) -> "off"=mean, "rm"=robbinsmonro, "sa"=simannealing; model("irt") -> "anr"=analytical newton-raphson, "nnr"=numerical newton-raphson #Convergence procedure
                estgain=1,       # Constant to slow down(decrease to decimal)/speed up(increase) newton cycles, or exponent on denominator in rm-squeeze
                burnin=500,      # MCMC Burn-In iterations... or some other convergence criteria, acf? or Euclidean distance?
                quad="manual",   # gauss-hermite
                nodes=15,        # nodes for quadrature
                gridbounds=c(-3.5,3.5), #manual quadrature bounds
                tmu=c(0,0),           # Prior mean, can be multivariate vector(Adim)
                tsigma=matrix(c(1,0,0,1),nrow=2,ncol=2),        # Prior sigma, can be multivariate matrix(Adim x Adim)
                eps=1e-4,        # Converged?
                nesttheta=10,    # if esttheta="mcmc", how many random samples?
                thetamap=TRUE,   # do an MAP estimate of Theta?
                thetaGen=Gen2D$THETA, # Did you simulate a new set of thetas? if so, give them to me. 
                impute=FALSE,    # Impute missing data?
                plots=PLOT,     # Show plots for diagnostics
                chains=5,        # How many chains to build? Diagnose convergence? Simultaneous MCMC chains?
                initialize="best", # "best", "random"
                record="on",     # "off"
                parallel=TRUE) #,  # True or false for parallel computation?
                # simfile=structure$simfile, # or NA
                defaultP$estfile=paste("FitFile_Test") 
  }
  if (is.na(parameters)) {
    parameters<-defaultP
  } else {
    for (i in names(defaultP)) {
      parameters<-parCheck(i,parameters,defaultP)
    }
  }
  parameters
}