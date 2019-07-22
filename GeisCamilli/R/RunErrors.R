RunErrors<-function(fitObject=NA,condition=NA,fitRDAfile=NA,CLT.start=0.15,CLT.end=0.8,
                    RMSE=NA,errorFilename="DATA",guessing=FALSE) {
  if (!is.na(condition)) { #condition<-"S3"; CLT.start=0.15; CLT.end=0.8
    source("CreateSimulationStructure.R")
    fit.dir<-"ConvergedModelFits"
    d<-condition  #<-"S1"
    SIM = 1
    sims=sim.list[[d]]$Reps
    # items=c(30,40,50)
    # examinees=c(2000,3000,4000)
    # burn = 100
    #b = 1; n = 1; j = 1; i = 1 # debug
    items=sim.list[[d]]$J
    examinees=sim.list[[d]]$N #,10000)
    
    simdir<-paste0(gen.dir,"/",d)
    SimList<-readRDS(paste0(simdir,"/",SFileString(sim.list[[d]],gen=TRUE),"_",SIM,".rds"))
    fitdir<-paste0(fit.dir,"/",d)
    FitDATA<-readRDS(paste0(fitdir,"/",SFileString(sim.list[[d]],gen=FALSE,r = SIM),".rds"))
    #load(paste0(fitdir,"/",SFileString(sim.list[[d]],gen=FALSE,r = SIM),".rda"))
    
    if (!guessing) {
      FitDATA$settings$guess<-FALSE
    }
    
    burninCLT = c(ceiling(CLT.start*FitDATA$settings$burnin),ceiling(CLT.end*FitDATA$settings$burnin))#,FitDATA$settings$burnin)
    burninMC = c(ceiling(0.6*FitDATA$settings$burnin),ceiling(0.8*FitDATA$settings$burnin))#,FitDATA$settings$burnin)
    burn = FitDATA$settings$burnin
    
    XIiter<-array(0, dim=c(items,1+FitDATA$settings$Adim+(FitDATA$settings$guess+0),sims))        
    ERRiter<-array(0, dim=c(items,1+FitDATA$settings$Adim+(FitDATA$settings$guess+0),sims)) # LMI iterative - 2PNO
    iERR<-array(0, dim=c(items,1+FitDATA$settings$Adim+(FitDATA$settings$guess+0),sims))    # LMI simple - 2PL
    oERR<-array(0, dim=c(items,1+FitDATA$settings$Adim+(FitDATA$settings$guess+0),sims))    # LMI simple - 2PNO
    MCERR.Burn<-array(0, dim=c(items,1+FitDATA$settings$Adim+(FitDATA$settings$guess+0),sims))  # Chain during burn
    # MCERR.Ann<-array(0, dim=c(items,1+FitDATA$settings$Adim+(FitDATA$settings$guess+0),sims))   # Chain during Ann window
    MCCLT.Burn<-array(0, dim=c(items,1+FitDATA$settings$Adim+(FitDATA$settings$guess+0),sims))  # CLT during burn
    # MCCLT.Ann<-array(0, dim=c(items,1+FitDATA$settings$Adim+(FitDATA$settings$guess+0),sims))   # CLT during Ann window
    # EMCERR<-array(0, dim=c(items,1+FitDATA$settings$Adim,sims))  # Chain empirical
    # EMCCLT<-array(0, dim=c(items,1+FitDATA$settings$Adim,sims))   # MCMC CLT empirical
    EmpLMIi<-array(0, dim=c(items,1+FitDATA$settings$Adim+(FitDATA$settings$guess+0),sims)) # LMI iterative - 2PL
    EmpLMIo<-array(0, dim=c(items,1+FitDATA$settings$Adim+(FitDATA$settings$guess+0),sims)) # LMI iterative - 2PNO
    for (i in 1:sims) {
      cat(":",i)
      # estfile=paste0(fitdir,"/",SFileString(sim.list[[d]],gen=FALSE,r = i),".rds")
      # load(estfile)
      FitDATA<-readRDS(paste0(fitdir,"/",SFileString(sim.list[[d]],gen=FALSE,r = i),".rds"))
      if (!guessing) {
        FitDATA$settings$guess<-FALSE
      }
      XIiter[,,i]<-as.matrix(FitDATA$xi[,1:(1+FitDATA$settings$Adim+(FitDATA$settings$guess+0))])   # converged
      ERRiter[,,i]<-as.matrix(FitDATA$xiError[,1:(1+FitDATA$settings$Adim+(FitDATA$settings$guess+0))])   # LMI iterative
      iERR[,,i]<-as.matrix(FitDATA$iError[,1:(1+FitDATA$settings$Adim+(FitDATA$settings$guess+0))])   # LMI at convergence, one calculation - 2pl
      oERR[,,i]<-as.matrix(FitDATA$oError[,1:(1+FitDATA$settings$Adim+(FitDATA$settings$guess+0))]) # LMI at convergence, one calculation - 2pno
      MCERR.Burn[,,i]<-as.matrix(mcmcTrimmedError(FitDATA$Iterations,end = burninMC[2],start = burninMC[1],settings=FitDATA$settings))   # Chain during burn
      #MCERR.Ann[,,i]<-as.matrix(mcmcTrimmedError(MCMCDATA,end = burninMC[3],start = burninMC[2]))    # Chain during Ann window
      MCCLT.Burn[,,i]<-as.matrix(mcmcCLTError(FitDATA$Iterations,end = burninCLT[2],start = burninCLT[1],settings=FitDATA$settings))  # CLT during burn
      #MCCLT.Ann[,,i]<-as.matrix(mcmcCLTError(MCMCDATA,end = burninCLT[3],start = burninCLT[2]))   # CLT during Ann window
      # EMCERR[,,i]<-as.matrix(cbind(FitDATA$EmpSE$MCSA,FitDATA$EmpSE$MCSB)) # Chain empirical
      # EMCCLT[,,i]<-as.matrix(cbind(FitDATA$EmpSE$SEA,FitDATA$EmpSE$SEB))   # MCMC CLT empirical
      EmpLMIi[,,i]<-as.matrix(FitDATA$EmpSE$VARLMIi[,1:(1+FitDATA$settings$Adim+(FitDATA$settings$guess+0))])  # LMI empirical - iterative - 2pl
      EmpLMIo[,,i]<-as.matrix(FitDATA$EmpSE$VARLMIo[,1:(1+FitDATA$settings$Adim+(FitDATA$settings$guess+0))])  # LMI empirical - iterative - 2pno
    }
    
    #### Means of the variance estimations
    mXI<-  apply(XIiter,c(1,2),mean)  # converged   :::   parameters
    mLMI<- apply(ERRiter,c(1,2),function(x) mean(x[x>0])) # LMI iterative   :::  variance
    miERR<-apply(iERR,c(1,2),function(x) mean(x[x>0]))    # LMI at convergence, one calculation - 2pl   ::: variance
    moERR<-apply(oERR,c(1,2),function(x) mean(x[x>0]))    # LMI at convergence, one calculation - 2pno  ::: variance
    mLMIi<-apply(EmpLMIi,c(1,2),function(x) mean(x[x>0]))  # LMI empirical - iterative - 2pl  ::: variance
    mLMIo<-apply(EmpLMIo,c(1,2),function(x) mean(x[x>0]))  # LMI empirical - iterative - 2pno ::: variance
    
    print("Iterative LMI Negative Variances:")
    tt<-table(unlist(lapply(apply(ERRiter,3,function(x) which(x<0,arr.ind = TRUE)),function(y) (y[,1]))))
    if (length(tt)>0) print(tt) else cat("0\n")
    print("Simple LMI 2PL  Negative Variances:")
    tt<-table(unlist(lapply(apply(iERR,3,function(x) which(x<0,arr.ind = TRUE)),function(y) (y[,1]))))
    if (length(tt)>0) print(tt) else cat("0\n")
    print("Simple LMI 2PNO Negative Variances:")
    tt<-table(unlist(lapply(apply(oERR,3,function(x) which(x<0,arr.ind = TRUE)),function(y) (y[,1]))))
    if (length(tt)>0) print(tt) else cat("0\n")
    print("Empirical LMI 2PL  Negative Variances:")
    tt<-table(unlist(lapply(apply(EmpLMIi,3,function(x) which(x<0,arr.ind = TRUE)),function(y) (y[,1]))))
    if (length(tt)>0) print(tt) else cat("0\n")
    print("Empirical LMI 2PNO Negative Variances:")
    tt<-table(unlist(lapply(apply(EmpLMIo,3,function(x) which(x<0,arr.ind = TRUE)),function(y) (y[,1]))))
    if (length(tt)>0) print(tt) else cat("0\n")
    
    print("MCMC CLT Burnin NaN Variances:")
    tt<-table(unlist(lapply(apply(MCCLT.Burn,3,function(x) which(is.nan(x),arr.ind = TRUE)),function(y) (y[,1]))))
    if (length(tt)>0) print(tt) else cat("0\n")
    # print("MCMC CLT Annealing NaN Variances:")
    # tt<-table(unlist(lapply(apply(MCCLT.Ann,3,function(x) which(is.nan(x),arr.ind = TRUE)),function(y) (y[,1]))))
    # if (length(tt)>0) print(tt) else cat("0\n")
    # print("Empirical MCMC CLT Burnin NaN Variances:")
    # tt<-table(unlist(lapply(apply(EMCCLT,3,function(x) which(is.nan(x),arr.ind = TRUE)),function(y) (y[,1]))))
    # if (length(tt)>0) print(tt) else cat("0\n")
    
    #### All SDs
    sXI<-apply(XIiter,c(1,2),sd)  # RMSE
    eLMI<-matrix(sqrt(mLMI),nrow=items,ncol=1+FitDATA$settings$Adim+(FitDATA$settings$guess+0))  # LMI iterative ::: sd
    eiERR<-matrix(sqrt(miERR),nrow=items,ncol=1+FitDATA$settings$Adim+(FitDATA$settings$guess+0)) # LMI at convergence, one calculation - 2pl  ::: sd
    eoERR<-matrix(sqrt(moERR),nrow=items,ncol=1+FitDATA$settings$Adim+(FitDATA$settings$guess+0)) # LMI at convergence, one calculation - 2pno ::: sd
    mMCE.Burn<- apply(MCERR.Burn,c(1,2),mean)   # Chain during burn         ::: sd
    # mMCE.Ann<-  apply(MCERR.Ann,c(1,2),mean)    # Chain during Ann window   ::: sd
    mCLT.Burn<- apply(MCCLT.Burn,c(1,2),function(x) mean(x[!is.nan(x)]))   # CLT during burn           ::: sd
    # mCLT.Ann<-  apply(MCCLT.Ann,c(1,2),function(x) mean(x[!is.nan(x)]))    # CLT during Ann window     ::: sd
    # mEMC<- apply(EMCERR,c(1,2),mean)   # Chain empirical    ::: sd
    # mECLT<-apply(EMCCLT,c(1,2),mean)   # MCMC CLT empirical ::: sd
    eLMIi<-matrix(sqrt(mLMIi),nrow=items,ncol=1+FitDATA$settings$Adim+(FitDATA$settings$guess+0)) # LMI empirical - iterative - 2pl  ::: sd
    eLMIo<-matrix(sqrt(mLMIo),nrow=items,ncol=1+FitDATA$settings$Adim+(FitDATA$settings$guess+0)) # LMI empirical - iterative - 2pno ::: sd
    
    gen.xi<-SimList$gen.xi[,1:(1+FitDATA$settings$Adim+(FitDATA$settings$guess+0))]
    bias<-mXI-gen.xi
    #######  SD of the variances
    sLMI <-apply(ERRiter,c(1,2),function(x) sd(sqrt(x[x>0])))
    siERR<-apply(iERR,c(1,2),function(x) sd(sqrt(x[x>0])))
    soERR<-apply(oERR,c(1,2),function(x) sd(sqrt(x[x>0])))
    # sMCE.Burn<-apply(MCERR.Burn,c(1,2),sd)
    # sMCE.Ann <-apply(MCERR.Ann,c(1,2),sd)
    # sCLT.Burn<-apply(MCCLT.Burn,c(1,2),sd)
    # sCLT.Ann <-apply(MCCLT.Ann,c(1,2),sd)
    # sEMC <-apply(EMCERR,c(1,2),sd)
    # sECLT<-apply(EMCCLT,c(1,2),sd)
    sLMIi<-apply(EmpLMIi,c(1,2),function(x) sd(sqrt(x[x>0])))
    sLMIo<-apply(EmpLMIo,c(1,2),function(x) sd(sqrt(x[x>0])))
    
    ErrorStats<-data.frame(N=rep(examinees,items),J=rep(items,items),Burn=rep(burn,items))
    for(p in 1:(FitDATA$settings$Adim+1+(FitDATA$settings$guess+0))) {
      
      pt<-ifelse(p<=FitDATA$settings$Adim,paste("A",p,sep=""),
                 ifelse(p==(FitDATA$settings$Adim+1),"B","C"))
      ErrorStats<-cbind(ErrorStats,gen.xi[,p])
      colnames(ErrorStats)[ncol(ErrorStats)]<-pt
      ErrorStats<-cbind(ErrorStats,mXI[,p])
      colnames(ErrorStats)[ncol(ErrorStats)]<-paste(pt,"est",sep="_")
      ErrorStats<-cbind(ErrorStats,bias[,p])
      colnames(ErrorStats)[ncol(ErrorStats)]<-paste(pt,"bias",sep="_")
      ErrorStats<-cbind(ErrorStats,sXI[,p])
      colnames(ErrorStats)[ncol(ErrorStats)]<-paste(pt,"RMSE",sep="_")
      ErrorStats<-cbind(ErrorStats,mLMI[,p])
      colnames(ErrorStats)[ncol(ErrorStats)]<-paste(pt,"VAR_LMI_Iter",sep="_")
      ErrorStats<-cbind(ErrorStats,eLMI[,p])
      colnames(ErrorStats)[ncol(ErrorStats)]<-paste(pt,"ERR_LMI_Iter",sep="_")
      ErrorStats<-cbind(ErrorStats,sLMI[,p])
      colnames(ErrorStats)[ncol(ErrorStats)]<-paste(pt,"sd_ERR_LMI_Iter",sep="_")
      ErrorStats<-cbind(ErrorStats,miERR[,p])
      colnames(ErrorStats)[ncol(ErrorStats)]<-paste(pt,"VAR_Simple_2PL",sep="_")
      ErrorStats<-cbind(ErrorStats,eiERR[,p])
      colnames(ErrorStats)[ncol(ErrorStats)]<-paste(pt,"ERR_Simple_2PL",sep="_")
      ErrorStats<-cbind(ErrorStats,siERR[,p])
      colnames(ErrorStats)[ncol(ErrorStats)]<-paste(pt,"sd_ERR_Simple_2PL",sep="_")
      ErrorStats<-cbind(ErrorStats,moERR[,p])
      colnames(ErrorStats)[ncol(ErrorStats)]<-paste(pt,"VAR_Simple_2PNO",sep="_")
      ErrorStats<-cbind(ErrorStats,eoERR[,p])
      colnames(ErrorStats)[ncol(ErrorStats)]<-paste(pt,"ERR_Simple_2PNO",sep="_")
      ErrorStats<-cbind(ErrorStats,soERR[,p])
      colnames(ErrorStats)[ncol(ErrorStats)]<-paste(pt,"sd_ERR_Simple_2PNO",sep="_")
      ErrorStats<-cbind(ErrorStats,mMCE.Burn[,p])
      colnames(ErrorStats)[ncol(ErrorStats)]<-paste(pt,"ERR_MCMC_Burn",sep="_")
      ErrorStats<-cbind(ErrorStats,sMCE.Burn[,p])
      colnames(ErrorStats)[ncol(ErrorStats)]<-paste(pt,"sd_ERR_MCMC_Burn",sep="_")
      # ErrorStats<-cbind(ErrorStats,mMCE.Ann[,p])
      # colnames(ErrorStats)[ncol(ErrorStats)]<-paste(pt,"ERR_MCMC_Ann",sep="_")
      # ErrorStats<-cbind(ErrorStats,sMCE.Ann[,p])
      # colnames(ErrorStats)[ncol(ErrorStats)]<-paste(pt,"sd_ERR_MCMC_Ann",sep="_")
      ErrorStats<-cbind(ErrorStats,mCLT.Burn[,p])
      colnames(ErrorStats)[ncol(ErrorStats)]<-paste(pt,"ERR_CLT_Burn",sep="_")
      ErrorStats<-cbind(ErrorStats,sCLT.Burn[,p])
      colnames(ErrorStats)[ncol(ErrorStats)]<-paste(pt,"sd_ERR_CLT_Burn",sep="_")
      # ErrorStats<-cbind(ErrorStats,mCLT.Ann[,p])
      # colnames(ErrorStats)[ncol(ErrorStats)]<-paste(pt,"ERR_CLT_Ann",sep="_")
      # ErrorStats<-cbind(ErrorStats,sCLT.Ann[,p])
      # colnames(ErrorStats)[ncol(ErrorStats)]<-paste(pt,"sd_ERR_CLT_Ann",sep="_")
      # ErrorStats<-cbind(ErrorStats,mEMC[,p])
      # colnames(ErrorStats)[ncol(ErrorStats)]<-paste(pt,"ERR_MCMC_Emp",sep="_")
      # ErrorStats<-cbind(ErrorStats,sEMC[,p])
      # colnames(ErrorStats)[ncol(ErrorStats)]<-paste(pt,"sd_ERR_MCMC_Emp",sep="_")
      # ErrorStats<-cbind(ErrorStats,mECLT[,p])
      # colnames(ErrorStats)[ncol(ErrorStats)]<-paste(pt,"ERR_CLT_Emp",sep="_")
      # ErrorStats<-cbind(ErrorStats,sECLT[,p])
      # colnames(ErrorStats)[ncol(ErrorStats)]<-paste(pt,"sd_ERR_CLT_Emp",sep="_")
      ErrorStats<-cbind(ErrorStats,mLMIi[,p])
      colnames(ErrorStats)[ncol(ErrorStats)]<-paste(pt,"VAR_LMI_Iter_Emp_2PL",sep="_")
      ErrorStats<-cbind(ErrorStats,eLMIi[,p])
      colnames(ErrorStats)[ncol(ErrorStats)]<-paste(pt,"ERR_LMI_Iter_Emp_2PL",sep="_")
      ErrorStats<-cbind(ErrorStats,sLMIi[,p])
      colnames(ErrorStats)[ncol(ErrorStats)]<-paste(pt,"sd_ERR_LMI_Iter_Emp_2PL",sep="_")
      ErrorStats<-cbind(ErrorStats,mLMIo[,p])
      colnames(ErrorStats)[ncol(ErrorStats)]<-paste(pt,"VAR_LMI_Iter_Emp_2PNO",sep="_")
      ErrorStats<-cbind(ErrorStats,eLMIo[,p])
      colnames(ErrorStats)[ncol(ErrorStats)]<-paste(pt,"ERR_LMI_Iter_Emp_2PNO",sep="_")
      ErrorStats<-cbind(ErrorStats,sLMIo[,p])
      colnames(ErrorStats)[ncol(ErrorStats)]<-paste(pt,"sd_ERR_LMI_Iter_Emp_2PNO",sep="_")
    } # else {
    write.csv(ErrorStats,paste0("ErrorStats_",d,".csv"))
  } else if (!is.na(fitObject) | !is.na(fitRDAfile)) {
    if (!is.na(fitRDAfile)) {
      load(fitRDAfile)
      if (!exists("FitDATA")) {
        print("Don't know if you've got the right rda file: NULL for YOU!")
        return(NULL)
      }
    } else {
      FitDATA<-fitObject
    }
    items=ncol(FitDATA$RP)
    examinees=nrow(FitDATA$RP)
    if (FitDATA$settings$record) {
      #MCMCDATA<-fitObject$MCMCDATA
      burninCLT = c(ceiling(CLT.start*FitDATA$settings$burnin),ceiling(CLT.end*FitDATA$settings$burnin))#,FitDATA$settings$burnin)
      burninMC = c(ceiling(0.6*FitDATA$settings$burnin),ceiling(0.8*FitDATA$settings$burnin))#,FitDATA$settings$burnin)
      burn = FitDATA$settings$burnin
      
      mMCE.Burn<-as.matrix(mcmcTrimmedError(FitDATA$Iterations,
                                             end = burninMC[2],start = burninMC[1],
                                             settings=FitDATA$settings))  # Chain during burn
      mCLT.Burn<-as.matrix(mcmcCLTError(FitDATA$Iterations,
                                         end = burninCLT[2],start = burninCLT[1],
                                         settings=FitDATA$settings))  # CLT during burn
      print("MCMC CLT Burnin NaN Variances:")
      tt<-table(unlist(lapply(apply(MCCLT.Burn,3,function(x) which(is.nan(x),arr.ind = TRUE)),function(y) (y[,1]))))
      if (length(tt)>0) print(tt) else cat("0\n")
    }
    
    # XIiter<-array(0, dim=c(items,1+FitDATA$settings$Adim+(FitDATA$settings$guess+0),sims))        
    # ERRiter<-array(0, dim=c(items,1+FitDATA$settings$Adim+(FitDATA$settings$guess+0),sims)) # LMI iterative - 2PNO
    # iERR<-array(0, dim=c(items,1+FitDATA$settings$Adim+(FitDATA$settings$guess+0),sims))    # LMI simple - 2PL
    # oERR<-array(0, dim=c(items,1+FitDATA$settings$Adim+(FitDATA$settings$guess+0),sims))    # LMI simple - 2PNO
    # EmpLMIi<-array(0, dim=c(items,1+FitDATA$settings$Adim+(FitDATA$settings$guess+0),sims)) # LMI iterative - 2PL
    # EmpLMIo<-array(0, dim=c(items,1+FitDATA$settings$Adim+(FitDATA$settings$guess+0),sims)) # LMI iterative - 2PNO
    
    # for (i in 1:sims) {
    #   cat(":",i)
    #   estfile=paste0(fitdir,"/",SFileString(sim.list[[d]],gen=FALSE,r = i),".rda")
    #   load(estfile)
    #   XIiter[,,i]<-as.matrix(FitDATA$xi)   # converged
    #   ERRiter[,,i]<-as.matrix(FitDATA$xiError)   # LMI iterative
    #   iERR[,,i]<-as.matrix(FitDATA$iError)   # LMI at convergence, one calculation - 2pl
    #   oERR[,,i]<-as.matrix(FitDATA$oError) # LMI at convergence, one calculation - 2pno
    #   MCERR.Burn[,,i]<-as.matrix(mcmcTrimmedError(MCMCDATA,end = burninMC[2],start = burninMC[1],settings=FitDATA$settings))   # Chain during burn
    #   #MCERR.Ann[,,i]<-as.matrix(mcmcTrimmedError(MCMCDATA,end = burninMC[3],start = burninMC[2]))    # Chain during Ann window
    #   MCCLT.Burn[,,i]<-as.matrix(mcmcCLTError(MCMCDATA,end = burninCLT[2],start = burninCLT[1],settings=FitDATA$settings))  # CLT during burn
    #   #MCCLT.Ann[,,i]<-as.matrix(mcmcCLTError(MCMCDATA,end = burninCLT[3],start = burninCLT[2]))   # CLT during Ann window
    #   # EMCERR[,,i]<-as.matrix(cbind(FitDATA$EmpSE$MCSA,FitDATA$EmpSE$MCSB)) # Chain empirical
    #   # EMCCLT[,,i]<-as.matrix(cbind(FitDATA$EmpSE$SEA,FitDATA$EmpSE$SEB))   # MCMC CLT empirical
    #   EmpLMIi[,,i]<-as.matrix(FitDATA$EmpSE$VARLMIi)  # LMI empirical - iterative - 2pl
    #   EmpLMIo[,,i]<-as.matrix(FitDATA$EmpSE$VARLMIo)  # LMI empirical - iterative - 2pno
    # }
    
    #### Means of the variance estimations
    mXI<-   FitDATA$xi #apply(XIiter,c(1,2),mean)  # converged   :::   parameters
    mLMI<-  FitDATA$xiError #apply(ERRiter,c(1,2),function(x) mean(x[x>0])) # LMI iterative   :::  variance
    miERR<- FitDATA$iError    #apply(iERR,c(1,2),function(x) mean(x[x>0]))    # LMI at convergence, one calculation - 2pl   ::: variance
    moERR<- FitDATA$oError #apply(oERR,c(1,2),function(x) mean(x[x>0]))    # LMI at convergence, one calculation - 2pno  ::: variance
    print("Iterative LMI Negative Variances:")
    tt<-which(FitDATA$xiError<0,arr.ind = TRUE)
    if (length(tt)>0) print(tt) else cat("0\n")
    print("Simple LMI 2PL  Negative Variances:")
    tt<-which(FitDATA$iError<0,arr.ind = TRUE)
    if (length(tt)>0) print(tt) else cat("0\n")
    print("Simple LMI 2PNO Negative Variances:")
    tt<-which(FitDATA$oError<0,arr.ind = TRUE)
    if (length(tt)>0) print(tt) else cat("0\n")
    if (FitDATA$settings$empiricalse) {
      EMCERR<-as.matrix(mcmcTrimmedError(FitDATA$EmpSE,
                                         end = FitDATA$settings$EmpIT,start = 1,
                                         settings=FitDATA$settings))  # Chain empirical
      EMCCLT<-as.matrix(mcmcCLTError(FitDATA$EmpSE,
                                     end = FitDATA$settings$EmpIT,start = 1,
                                     settings=FitDATA$settings))   # MCMC CLT empirical
      mLMIi<- FitDATA$EmpSE$VARLMIi #apply(EmpLMIi,c(1,2),function(x) mean(x[x>0]))  # LMI empirical - iterative - 2pl  ::: variance
      mLMIo<- FitDATA$EmpSE$VARLMIo #apply(EmpLMIo,c(1,2),function(x) mean(x[x>0]))  # LMI empirical - iterative - 2pno ::: variance
      print("Empirical LMI 2PL  Negative Variances:")
      tt<-which(mLMIi<0,arr.ind = TRUE)
      if (length(tt)>0) print(tt) else cat("0\n")
      print("Empirical LMI 2PNO Negative Variances:")
      tt<-which(mLMIo<0,arr.ind = TRUE)
      if (length(tt)>0) print(tt) else cat("0\n")
      eLMIi<-sqrt(mLMIi) # LMI empirical - iterative - 2pl  ::: sd
      eLMIo<-sqrt(mLMIo) # LMI empirical - iterative - 2pno ::: sd
    }
    
    
    #### All SDs
    eLMI<-sqrt(mLMI) # LMI iterative ::: sd
    eiERR<-sqrt(miERR)# LMI at convergence, one calculation - 2pl  ::: sd
    eoERR<-sqrt(moERR) # LMI at convergence, one calculation - 2pno ::: sd
  
    
    #######  SD of the variances
    # sLMI <-apply(ERRiter,c(1,2),function(x) sd(sqrt(x[x>0])))
    # siERR<-apply(iERR,c(1,2),function(x) sd(sqrt(x[x>0])))
    # soERR<-apply(oERR,c(1,2),function(x) sd(sqrt(x[x>0])))
    # sMCE.Burn<-apply(MCERR.Burn,c(1,2),sd)
    # sMCE.Ann <-apply(MCERR.Ann,c(1,2),sd)
    # sCLT.Burn<-apply(MCCLT.Burn,c(1,2),sd)
    # sCLT.Ann <-apply(MCCLT.Ann,c(1,2),sd)
    # sEMC <-apply(EMCERR,c(1,2),sd)
    # sECLT<-apply(EMCCLT,c(1,2),sd)
    # sLMIi<-apply(EmpLMIi,c(1,2),function(x) sd(sqrt(x[x>0])))
    # sLMIo<-apply(EmpLMIo,c(1,2),function(x) sd(sqrt(x[x>0])))
    if (!is.na(RMSE)[1]) {
      if (ncol(RMSE)<ncol(FitDATA$xi) | nrow(RMSE)!=nrow(FitDATA$xi)) {
        print("Check your RMSEs: NULL for you.")
        return(NULL)
      } else if (ncol(RMSE)>ncol(FitDATA$xi)) {
        print("Check your RMSEs")
      }
    }
    ErrorStats<-data.frame(N=rep(examinees,items),J=rep(items,items),Burn=rep(burn,items))
    for(p in 1:ncol(FitDATA$xi)) {
      pt<-ifelse(p<=FitDATA$settings$Adim,paste("A",p,sep=""),
                 ifelse(p==(FitDATA$settings$Adim+1),"B","C"))
      ErrorStats<-cbind(ErrorStats,gen.xi[,p])
      colnames(ErrorStats)[ncol(ErrorStats)]<-pt
      ErrorStats<-cbind(ErrorStats,mXI[,p])
      colnames(ErrorStats)[ncol(ErrorStats)]<-paste(pt,"est",sep="_")
      ErrorStats<-cbind(ErrorStats,NA)
      colnames(ErrorStats)[ncol(ErrorStats)]<-paste(pt,"bias",sep="_")
      if (!is.na(RMSE)[1]) {
        ErrorStats<-cbind(ErrorStats,RMSE[,p])
        colnames(ErrorStats)[ncol(ErrorStats)]<-paste(pt,"RMSE",sep="_")
      } else {
        ErrorStats<-cbind(ErrorStats,NA)
        colnames(ErrorStats)[ncol(ErrorStats)]<-paste(pt,"RMSE",sep="_")
      }
      ErrorStats<-cbind(ErrorStats,mLMI[,p])
      colnames(ErrorStats)[ncol(ErrorStats)]<-paste(pt,"VAR_LMI_Iter",sep="_")
      ErrorStats<-cbind(ErrorStats,eLMI[,p])
      colnames(ErrorStats)[ncol(ErrorStats)]<-paste(pt,"ERR_LMI_Iter",sep="_")
      ErrorStats<-cbind(ErrorStats,NA)
      colnames(ErrorStats)[ncol(ErrorStats)]<-paste(pt,"sd_ERR_LMI_Iter",sep="_")
      ErrorStats<-cbind(ErrorStats,miERR[,p])
      colnames(ErrorStats)[ncol(ErrorStats)]<-paste(pt,"VAR_Simple_2PL",sep="_")
      ErrorStats<-cbind(ErrorStats,eiERR[,p])
      colnames(ErrorStats)[ncol(ErrorStats)]<-paste(pt,"ERR_Simple_2PL",sep="_")
      ErrorStats<-cbind(ErrorStats,NA)
      colnames(ErrorStats)[ncol(ErrorStats)]<-paste(pt,"sd_ERR_Simple_2PL",sep="_")
      ErrorStats<-cbind(ErrorStats,moERR[,p])
      colnames(ErrorStats)[ncol(ErrorStats)]<-paste(pt,"VAR_Simple_2PNO",sep="_")
      ErrorStats<-cbind(ErrorStats,eoERR[,p])
      colnames(ErrorStats)[ncol(ErrorStats)]<-paste(pt,"ERR_Simple_2PNO",sep="_")
      ErrorStats<-cbind(ErrorStats,NA)
      colnames(ErrorStats)[ncol(ErrorStats)]<-paste(pt,"sd_ERR_Simple_2PNO",sep="_")
      ErrorStats<-cbind(ErrorStats,mMCE.Burn[,p])
      colnames(ErrorStats)[ncol(ErrorStats)]<-paste(pt,"ERR_MCMC_Burn",sep="_")
      ErrorStats<-cbind(ErrorStats,NA)
      colnames(ErrorStats)[ncol(ErrorStats)]<-paste(pt,"sd_ERR_MCMC_Burn",sep="_")
      # ErrorStats<-cbind(ErrorStats,mMCE.Ann[,p])
      # colnames(ErrorStats)[ncol(ErrorStats)]<-paste(pt,"ERR_MCMC_Ann",sep="_")
      # ErrorStats<-cbind(ErrorStats,sMCE.Ann[,p])
      # colnames(ErrorStats)[ncol(ErrorStats)]<-paste(pt,"sd_ERR_MCMC_Ann",sep="_")
      ErrorStats<-cbind(ErrorStats,mCLT.Burn[,p])
      colnames(ErrorStats)[ncol(ErrorStats)]<-paste(pt,"ERR_CLT_Burn",sep="_")
      ErrorStats<-cbind(ErrorStats,NA)
      colnames(ErrorStats)[ncol(ErrorStats)]<-paste(pt,"sd_ERR_CLT_Burn",sep="_")
      # ErrorStats<-cbind(ErrorStats,mCLT.Ann[,p])
      # colnames(ErrorStats)[ncol(ErrorStats)]<-paste(pt,"ERR_CLT_Ann",sep="_")
      # ErrorStats<-cbind(ErrorStats,sCLT.Ann[,p])
      # colnames(ErrorStats)[ncol(ErrorStats)]<-paste(pt,"sd_ERR_CLT_Ann",sep="_")
      if (FitDATA$settings$empiricalse) {
        ErrorStats<-cbind(ErrorStats,mEMC[,p])
        colnames(ErrorStats)[ncol(ErrorStats)]<-paste(pt,"ERR_MCMC_Emp",sep="_")
        ErrorStats<-cbind(ErrorStats,sEMC[,p])
        colnames(ErrorStats)[ncol(ErrorStats)]<-paste(pt,"sd_ERR_MCMC_Emp",sep="_")
        ErrorStats<-cbind(ErrorStats,mECLT[,p])
        colnames(ErrorStats)[ncol(ErrorStats)]<-paste(pt,"ERR_CLT_Emp",sep="_")
        ErrorStats<-cbind(ErrorStats,sECLT[,p])
        colnames(ErrorStats)[ncol(ErrorStats)]<-paste(pt,"sd_ERR_CLT_Emp",sep="_")
        ErrorStats<-cbind(ErrorStats,mLMIi[,p])
        colnames(ErrorStats)[ncol(ErrorStats)]<-paste(pt,"VAR_LMI_Iter_Emp_2PL",sep="_")
        ErrorStats<-cbind(ErrorStats,eLMIi[,p])
        colnames(ErrorStats)[ncol(ErrorStats)]<-paste(pt,"ERR_LMI_Iter_Emp_2PL",sep="_")
        ErrorStats<-cbind(ErrorStats,NA)
        colnames(ErrorStats)[ncol(ErrorStats)]<-paste(pt,"sd_ERR_LMI_Iter_Emp_2PL",sep="_")
        ErrorStats<-cbind(ErrorStats,mLMIo[,p])
        colnames(ErrorStats)[ncol(ErrorStats)]<-paste(pt,"VAR_LMI_Iter_Emp_2PNO",sep="_")
        ErrorStats<-cbind(ErrorStats,eLMIo[,p])
        colnames(ErrorStats)[ncol(ErrorStats)]<-paste(pt,"ERR_LMI_Iter_Emp_2PNO",sep="_")
        ErrorStats<-cbind(ErrorStats,NA)
        colnames(ErrorStats)[ncol(ErrorStats)]<-paste(pt,"sd_ERR_LMI_Iter_Emp_2PNO",sep="_")
      } else {
        ErrorStats<-cbind(ErrorStats,NA)
        colnames(ErrorStats)[ncol(ErrorStats)]<-paste(pt,"ERR_MCMC_Emp",sep="_")
        ErrorStats<-cbind(ErrorStats,NA)
        colnames(ErrorStats)[ncol(ErrorStats)]<-paste(pt,"sd_ERR_MCMC_Emp",sep="_")
        ErrorStats<-cbind(ErrorStats,NA)
        colnames(ErrorStats)[ncol(ErrorStats)]<-paste(pt,"ERR_CLT_Emp",sep="_")
        ErrorStats<-cbind(ErrorStats,NA)
        colnames(ErrorStats)[ncol(ErrorStats)]<-paste(pt,"sd_ERR_CLT_Emp",sep="_")
        ErrorStats<-cbind(ErrorStats,NA)
        colnames(ErrorStats)[ncol(ErrorStats)]<-paste(pt,"VAR_LMI_Iter_Emp_2PL",sep="_")
        ErrorStats<-cbind(ErrorStats,NA)
        colnames(ErrorStats)[ncol(ErrorStats)]<-paste(pt,"ERR_LMI_Iter_Emp_2PL",sep="_")
        ErrorStats<-cbind(ErrorStats,NA)
        colnames(ErrorStats)[ncol(ErrorStats)]<-paste(pt,"sd_ERR_LMI_Iter_Emp_2PL",sep="_")
        ErrorStats<-cbind(ErrorStats,NA)
        colnames(ErrorStats)[ncol(ErrorStats)]<-paste(pt,"VAR_LMI_Iter_Emp_2PNO",sep="_")
        ErrorStats<-cbind(ErrorStats,NA)
        colnames(ErrorStats)[ncol(ErrorStats)]<-paste(pt,"ERR_LMI_Iter_Emp_2PNO",sep="_")
        ErrorStats<-cbind(ErrorStats,NA)
        colnames(ErrorStats)[ncol(ErrorStats)]<-paste(pt,"sd_ERR_LMI_Iter_Emp_2PNO",sep="_")
      }
    } 
    write.csv(ErrorStats,paste0("ErrorStats_",errorFilename,".csv"))
  } else {
    print("Don't know what to do: NULL for you!")
    return(NULL)
  }
  ErrorStats
}