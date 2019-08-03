RunErrors<-function(fitObject=NA,condition=NA,fitFile=NA,CLT.start=0.15,CLT.end=0.8,
                    RMSE=NA,errorFilename="DATA",guessing=FALSE,item.ind=FALSE) {
  if ((is.na(fitObject) & is.na(fitFile)) & !is.na(condition)) { #condition<-"S3"; CLT.start=0.3; CLT.end=0.8;item.ind=TRUE
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
    
    if ("tau" %in% names(FitDATA)) {
      TAUiter<-array(0, dim=c(items,FitDATA$settings$ncat-1,sims))
    }
    
    if (!is.na(FitDATA$xiError)[1]) {
      ERRiter<-array(0, dim=c(items,1+FitDATA$settings$Adim+(FitDATA$settings$guess+0),sims)) # LMI iterative - 2PNO
      iERR<-array(0, dim=c(items,1+FitDATA$settings$Adim+(FitDATA$settings$guess+0),sims))    # LMI simple - 2PL
      oERR<-array(0, dim=c(items,1+FitDATA$settings$Adim+(FitDATA$settings$guess+0),sims))    # LMI simple - 2PNO
    }
    MCERR.Burn<-array(0, dim=c(items,1+FitDATA$settings$Adim+(FitDATA$settings$guess+0),sims))  # Chain during burn
    # MCERR.Ann<-array(0, dim=c(items,1+FitDATA$settings$Adim+(FitDATA$settings$guess+0),sims))   # Chain during Ann window
    MCCLT.Burn<-array(0, dim=c(items,1+FitDATA$settings$Adim+(FitDATA$settings$guess+0),sims))  # CLT during burn
    # MCCLT.Ann<-array(0, dim=c(items,1+FitDATA$settings$Adim+(FitDATA$settings$guess+0),sims))   # CLT during Ann window
    
    if ("tau" %in% names(FitDATA) && "Diter" %in% names(FitDATA$Iterations)) {
      MCERR.Burn.d<-array(0, dim=c(items,FitDATA$settings$ncat-1,sims))  # Chain during burn
      # MCERR.Ann<-array(0, dim=c(items,1+FitDATA$settings$Adim+(FitDATA$settings$guess+0),sims))   # Chain during Ann window
      MCCLT.Burn.d<-array(0, dim=c(items,FitDATA$settings$ncat-1,sims))  # CLT during burn
      # MCCLT.Ann<-array(0, dim=c(items,1+FitDATA$settings$Adim+(FitDATA$settings$guess+0),sims))   # CLT during Ann window
    }
    
    
    # EMCERR<-array(0, dim=c(items,1+FitDATA$settings$Adim,sims))  # Chain empirical
    # EMCCLT<-array(0, dim=c(items,1+FitDATA$settings$Adim,sims))   # MCMC CLT empirical
    if (!is.null(FitDATA$EmpSE$VARLMIi)[1] && !is.na(FitDATA$EmpSE$VARLMIi)[1]) {
      EmpLMIi<-array(0, dim=c(items,1+FitDATA$settings$Adim+(FitDATA$settings$guess+0),sims)) # LMI iterative - 2PL
      EmpLMIo<-array(0, dim=c(items,1+FitDATA$settings$Adim+(FitDATA$settings$guess+0),sims)) # LMI iterative - 2PNO
    }
    for (i in 1:sims) {#i<-2
      cat(":",i)
      # estfile=paste0(fitdir,"/",SFileString(sim.list[[d]],gen=FALSE,r = i),".rds")
      # load(estfile)
      FitDATA<-readRDS(paste0(fitdir,"/",SFileString(sim.list[[d]],gen=FALSE,r = i),".rds"))
      if (!guessing) {
        FitDATA$settings$guess<-FALSE
      }
      XIiter[,,i]<-as.matrix(FitDATA$xi[,1:(1+FitDATA$settings$Adim+(FitDATA$settings$guess+0))])   # converged
      if ("tau" %in% names(FitDATA)) {
        TAUiter[,,i]<-as.matrix(FitDATA$tau)
      }
      if (!is.na(FitDATA$xiError)[1]) {
        ERRiter[,,i]<-as.matrix(FitDATA$xiError[,1:(1+FitDATA$settings$Adim+(FitDATA$settings$guess+0))])   # LMI iterative
        iERR[,,i]<-as.matrix(FitDATA$iError[,1:(1+FitDATA$settings$Adim+(FitDATA$settings$guess+0))])   # LMI at convergence, one calculation - 2pl
        oERR[,,i]<-as.matrix(FitDATA$oError[,1:(1+FitDATA$settings$Adim+(FitDATA$settings$guess+0))]) # LMI at convergence, one calculation - 2pno
      }
      if (dim(FitDATA$Iterations$Aiter)[3]>burninMC[2]) {
        MCERR.Burn[,,i]<-as.matrix(mcmcTrimmedError(FitDATA$Iterations,end = burninMC[2],start = burninMC[1],settings=FitDATA$settings))   # Chain during burn
      
        #MCERR.Ann[,,i]<-as.matrix(mcmcTrimmedError(MCMCDATA,end = burninMC[3],start = burninMC[2]))    # Chain during Ann window
        MCCLT.Burn[,,i]<-as.matrix(mcmcCLTError(FitDATA$Iterations,item.ind=item.ind,
                                                end = burninCLT[2],
                                                start = burninCLT[1],
                                                settings=FitDATA$settings))  # CLT during burn
        if ("tau" %in% names(FitDATA) && "Diter" %in% names(FitDATA$Iterations)) {
          MCERR.Burn.d[,,i]<-as.matrix(mcmcTrimmedError(FitDATA$Iterations,
                                                        end = burninMC[2],
                                                        start = burninMC[1],
                                                        settings=FitDATA$settings,tau=TRUE))   # Chain during burn
          # MCERR.Ann<-array(0, dim=c(items,1+FitDATA$settings$Adim+(FitDATA$settings$guess+0),sims))   # Chain during Ann window
          MCCLT.Burn.d[,,i]<-as.matrix(mcmcCLTError(FitDATA$Iterations,item.ind=item.ind,
                                                    end = burninCLT[2],
                                                    start = burninCLT[1],
                                                    settings=FitDATA$settings,tau=TRUE))  # CLT during burn
          # MCCLT.Ann<-array(0, dim=c(items,1+FitDATA$settings$Adim+(FitDATA$settings$guess+0),sims))   # CLT during Ann window
        }
      } else {
        MCERR.Burn[,,i]<-NA
        MCCLT.Burn[,,i]<-NA
        if ("tau" %in% names(FitDATA) && "Diter" %in% names(FitDATA$Iterations)) {
          MCERR.Burn.d[,,i]<-NA
          # MCERR.Ann<-NA
          MCCLT.Burn.d[,,i]<-NA
          # MCCLT.Ann<-NA
        }
      }
      
      #MCCLT.Ann[,,i]<-as.matrix(mcmcCLTError(MCMCDATA,end = burninCLT[3],start = burninCLT[2]))   # CLT during Ann window
      # EMCERR[,,i]<-as.matrix(cbind(FitDATA$EmpSE$MCSA,FitDATA$EmpSE$MCSB)) # Chain empirical
      # EMCCLT[,,i]<-as.matrix(cbind(FitDATA$EmpSE$SEA,FitDATA$EmpSE$SEB))   # MCMC CLT empirical
      if (!is.null(FitDATA$EmpSE$VARLMIi)[1] && !is.na(FitDATA$EmpSE$VARLMIi)[1]) {
        EmpLMIi[,,i]<-as.matrix(FitDATA$EmpSE$VARLMIi[,1:(1+FitDATA$settings$Adim+(FitDATA$settings$guess+0))])  # LMI empirical - iterative - 2pl
        EmpLMIo[,,i]<-as.matrix(FitDATA$EmpSE$VARLMIo[,1:(1+FitDATA$settings$Adim+(FitDATA$settings$guess+0))])  # LMI empirical - iterative - 2pno
      }
    }
    dim(MCERR.Burn.d)
    MCCLT.Burn.d[1,1,]
    #### Means of the variance estimations
    mXI<-  apply(XIiter,c(1,2),mean)  # converged   :::   parameters
    if ("tau" %in% names(FitDATA)) {
      mTAU<-  apply(TAUiter,c(1,2),mean)  # converged   :::   parameters
    }  
    if (!is.na(FitDATA$xiError)[1]) {
      mLMI<- apply(ERRiter,c(1,2),function(x) mean(x[x>0])) # LMI iterative   :::  variance
      miERR<-apply(iERR,c(1,2),function(x) mean(x[x>0]))    # LMI at convergence, one calculation - 2pl   ::: variance
      moERR<-apply(oERR,c(1,2),function(x) mean(x[x>0]))    # LMI at convergence, one calculation - 2pno  ::: variance
    } else {
      mLMI<- NA
      miERR<-NA
      moERR<-NA
    } 
    if (!is.null(FitDATA$EmpSE$VARLMIi)[1] && !is.na(FitDATA$EmpSE$VARLMIi)[1]) {
      mLMIi<-apply(EmpLMIi,c(1,2),function(x) mean(x[x>0]))  # LMI empirical - iterative - 2pl  ::: variance
      mLMIo<-apply(EmpLMIo,c(1,2),function(x) mean(x[x>0]))  # LMI empirical - iterative - 2pno ::: variance
    } else {
      mLMIi<-NA
      mLMIo<-NA
    } 
    
    if (!is.na(FitDATA$xiError)[1]) {
      print("Iterative LMI Negative Variances:")
      tt<-table(unlist(lapply(apply(ERRiter,3,function(x) which(x<0,arr.ind = TRUE)),function(y) (y[,1]))))
      if (length(tt)>0) print(tt) else cat("0\n")
      print("Simple LMI 2PL  Negative Variances:")
      tt<-table(unlist(lapply(apply(iERR,3,function(x) which(x<0,arr.ind = TRUE)),function(y) (y[,1]))))
      if (length(tt)>0) print(tt) else cat("0\n")
      print("Simple LMI 2PNO Negative Variances:")
      tt<-table(unlist(lapply(apply(oERR,3,function(x) which(x<0,arr.ind = TRUE)),function(y) (y[,1]))))
      if (length(tt)>0) print(tt) else cat("0\n")
    }
    if (!is.null(FitDATA$EmpSE$VARLMIi)[1] && !is.na(FitDATA$EmpSE$VARLMIi)[1]) {
      print("Empirical LMI 2PL  Negative Variances:")
      tt<-table(unlist(lapply(apply(EmpLMIi,3,function(x) which(x<0,arr.ind = TRUE)),function(y) (y[,1]))))
      if (length(tt)>0) print(tt) else cat("0\n")
      print("Empirical LMI 2PNO Negative Variances:")
      tt<-table(unlist(lapply(apply(EmpLMIo,3,function(x) which(x<0,arr.ind = TRUE)),function(y) (y[,1]))))
      if (length(tt)>0) print(tt) else cat("0\n")
    }
    
    print("MCMC CLT Burnin NaN Variances:")
    tt<-table(unlist(lapply(apply(MCCLT.Burn,3,function(x) which(is.nan(x),arr.ind = TRUE)),function(y) (y[,1]))))
    if (length(tt)>0) print(tt) else cat("0\n")
    # print("MCMC CLT Annealing NaN Variances:")
    # tt<-table(unlist(lapply(apply(MCCLT.Ann,3,function(x) which(is.nan(x),arr.ind = TRUE)),function(y) (y[,1]))))
    # if (length(tt)>0) print(tt) else cat("0\n")
    # print("Empirical MCMC CLT Burnin NaN Variances:")
    # tt<-table(unlist(lapply(apply(EMCCLT,3,function(x) which(is.nan(x),arr.ind = TRUE)),function(y) (y[,1]))))
    # if (length(tt)>0) print(tt) else cat("0\n")
    
    gen.xi<-SimList$gen.xi[,1:(1+FitDATA$settings$Adim+(FitDATA$settings$guess+0))]
    bias<-mXI-gen.xi
    sXI<-apply(XIiter,c(1,2),sd)  # RMSE
    
    if ("tau" %in% names(FitDATA)) {
      gen.tau<-SimList$gen.tau
      tau.bias<-mTAU-gen.tau
      sTAU<-apply(TAUiter,c(1,2),sd)  # RMSE
    }
    
    #### All SDs
    if (!is.na(mLMI)[1]) {
      eLMI<-matrix(sqrt(mLMI),nrow=items,ncol=1+FitDATA$settings$Adim+(FitDATA$settings$guess+0))  # LMI iterative ::: sd
      sLMI <-apply(ERRiter,c(1,2),function(x) sd(sqrt(x[x>0])))
    } else {
      eLMI<-NA
      sLMI <-NA
    }
    if (!is.na(miERR)[1]) {
      eiERR<-matrix(sqrt(miERR),nrow=items,ncol=1+FitDATA$settings$Adim+(FitDATA$settings$guess+0)) # LMI at convergence, one calculation - 2pl  ::: sd
      eoERR<-matrix(sqrt(moERR),nrow=items,ncol=1+FitDATA$settings$Adim+(FitDATA$settings$guess+0)) # LMI at convergence, one calculation - 2pno ::: sd
      siERR<-apply(iERR,c(1,2),function(x) sd(sqrt(x[x>0])))
      soERR<-apply(oERR,c(1,2),function(x) sd(sqrt(x[x>0])))
    } else {
      eiERR<-NA
      eoERR<-NA
      siERR<-NA
      soERR<-NA
    } 
    mMCE.Burn<- apply(MCERR.Burn,c(1,2),mean,na.rm=TRUE)   # Chain during burn         ::: sd
    # mMCE.Ann<-  apply(MCERR.Ann,c(1,2),mean)    # Chain during Ann window   ::: sd
    mCLT.Burn<- apply(MCCLT.Burn,c(1,2),function(x) mean(x[!is.nan(x)],na.rm=TRUE))   # CLT during burn           ::: sd
    # mCLT.Ann<-  apply(MCCLT.Ann,c(1,2),function(x) mean(x[!is.nan(x)]))    # CLT during Ann window     ::: sd
    if ("tau" %in% names(FitDATA) && "Diter" %in% names(FitDATA$Iterations)) {
      mMCE.Burn.d<- apply(MCERR.Burn.d,c(1,2),mean,na.rm=TRUE)   # Chain during burn         ::: sd
      mCLT.Burn.d<- apply(MCCLT.Burn.d,c(1,2),function(x) mean(x[!is.nan(x)],na.rm=TRUE))   # CLT during burn           ::: sd
    }
    # mEMC<- apply(EMCERR,c(1,2),mean)   # Chain empirical    ::: sd
    # mECLT<-apply(EMCCLT,c(1,2),mean)   # MCMC CLT empirical ::: sd
    if (!is.na(mLMIi)[1]) {
      eLMIi<-matrix(sqrt(mLMIi),nrow=items,ncol=1+FitDATA$settings$Adim+(FitDATA$settings$guess+0)) # LMI empirical - iterative - 2pl  ::: sd
      eLMIo<-matrix(sqrt(mLMIo),nrow=items,ncol=1+FitDATA$settings$Adim+(FitDATA$settings$guess+0)) # LMI empirical - iterative - 2pno ::: sd
      sLMIi<-apply(EmpLMIi,c(1,2),function(x) sd(sqrt(x[x>0])))
      sLMIo<-apply(EmpLMIo,c(1,2),function(x) sd(sqrt(x[x>0])))
    } else {
      eLMIi<-NA
      eLMIo<-NA
      sLMIi<-NA
      sLMIo<-NA
    }
    
    
    #######  SD of the variances
    # sMCE.Burn<-apply(MCERR.Burn,c(1,2),sd)
    # sMCE.Ann <-apply(MCERR.Ann,c(1,2),sd)
    # sCLT.Burn<-apply(MCCLT.Burn,c(1,2),sd)
    # sCLT.Ann <-apply(MCCLT.Ann,c(1,2),sd)
    # sEMC <-apply(EMCERR,c(1,2),sd)
    # sECLT<-apply(EMCCLT,c(1,2),sd)
    
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
      if (!is.na(mLMI)[1]) {
        ErrorStats<-cbind(ErrorStats,mLMI[,p])
        colnames(ErrorStats)[ncol(ErrorStats)]<-paste(pt,"VAR_LMI_Iter",sep="_")
        ErrorStats<-cbind(ErrorStats,eLMI[,p])
        colnames(ErrorStats)[ncol(ErrorStats)]<-paste(pt,"ERR_LMI_Iter",sep="_")
        ErrorStats<-cbind(ErrorStats,sLMI[,p])
        colnames(ErrorStats)[ncol(ErrorStats)]<-paste(pt,"sd_ERR_LMI_Iter",sep="_")
      } else {
        ErrorStats<-cbind(ErrorStats,NA)
        colnames(ErrorStats)[ncol(ErrorStats)]<-paste(pt,"VAR_LMI_Iter",sep="_")
        ErrorStats<-cbind(ErrorStats,NA)
        colnames(ErrorStats)[ncol(ErrorStats)]<-paste(pt,"ERR_LMI_Iter",sep="_")
        ErrorStats<-cbind(ErrorStats,NA)
        colnames(ErrorStats)[ncol(ErrorStats)]<-paste(pt,"sd_ERR_LMI_Iter",sep="_")
      }
      if (!is.na(miERR)[1]) {
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
      } else {
        ErrorStats<-cbind(ErrorStats,NA)
        colnames(ErrorStats)[ncol(ErrorStats)]<-paste(pt,"VAR_Simple_2PL",sep="_")
        ErrorStats<-cbind(ErrorStats,NA)
        colnames(ErrorStats)[ncol(ErrorStats)]<-paste(pt,"ERR_Simple_2PL",sep="_")
        ErrorStats<-cbind(ErrorStats,NA)
        colnames(ErrorStats)[ncol(ErrorStats)]<-paste(pt,"sd_ERR_Simple_2PL",sep="_")
        ErrorStats<-cbind(ErrorStats,NA)
        colnames(ErrorStats)[ncol(ErrorStats)]<-paste(pt,"VAR_Simple_2PNO",sep="_")
        ErrorStats<-cbind(ErrorStats,NA)
        colnames(ErrorStats)[ncol(ErrorStats)]<-paste(pt,"ERR_Simple_2PNO",sep="_")
        ErrorStats<-cbind(ErrorStats,NA)
        colnames(ErrorStats)[ncol(ErrorStats)]<-paste(pt,"sd_ERR_Simple_2PNO",sep="_")
      } 
      ErrorStats<-cbind(ErrorStats,mMCE.Burn[,p])
      colnames(ErrorStats)[ncol(ErrorStats)]<-paste(pt,"ERR_MCMC_Burn",sep="_")
      # ErrorStats<-cbind(ErrorStats,sMCE.Burn[,p])
      # colnames(ErrorStats)[ncol(ErrorStats)]<-paste(pt,"sd_ERR_MCMC_Burn",sep="_")
      # ErrorStats<-cbind(ErrorStats,mMCE.Ann[,p])
      # colnames(ErrorStats)[ncol(ErrorStats)]<-paste(pt,"ERR_MCMC_Ann",sep="_")
      # ErrorStats<-cbind(ErrorStats,sMCE.Ann[,p])
      # colnames(ErrorStats)[ncol(ErrorStats)]<-paste(pt,"sd_ERR_MCMC_Ann",sep="_")
      ErrorStats<-cbind(ErrorStats,mCLT.Burn[,p])
      colnames(ErrorStats)[ncol(ErrorStats)]<-paste(pt,"ERR_CLT_Burn",sep="_")
      # ErrorStats<-cbind(ErrorStats,sCLT.Burn[,p])
      # colnames(ErrorStats)[ncol(ErrorStats)]<-paste(pt,"sd_ERR_CLT_Burn",sep="_")
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
      if (!is.na(mLMIi)[1]) {
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
      } else { 
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
    } # else {
    
    if ("tau" %in% names(FitDATA)) {
      for (p in 1:(FitDATA$settings$ncat-1)) {
      pt<-paste0("Tau",p)
      ErrorStats<-cbind(ErrorStats,gen.tau[,p])
      colnames(ErrorStats)[ncol(ErrorStats)]<-pt
      ErrorStats<-cbind(ErrorStats,mTAU[,p])
      colnames(ErrorStats)[ncol(ErrorStats)]<-paste(pt,"est",sep="_")
      ErrorStats<-cbind(ErrorStats,tau.bias[,p])
      colnames(ErrorStats)[ncol(ErrorStats)]<-paste(pt,"bias",sep="_")
      ErrorStats<-cbind(ErrorStats,sTAU[,p])
      colnames(ErrorStats)[ncol(ErrorStats)]<-paste(pt,"RMSE",sep="_")
      ErrorStats<-cbind(ErrorStats,mMCE.Burn.d[,p])
      colnames(ErrorStats)[ncol(ErrorStats)]<-paste(pt,"ERR_MCMC_Burn",sep="_")
      # ErrorStats<-cbind(ErrorStats,sMCE.Burn[,p])
      # colnames(ErrorStats)[ncol(ErrorStats)]<-paste(pt,"sd_ERR_MCMC_Burn",sep="_")
      # ErrorStats<-cbind(ErrorStats,mMCE.Ann[,p])
      # colnames(ErrorStats)[ncol(ErrorStats)]<-paste(pt,"ERR_MCMC_Ann",sep="_")
      # ErrorStats<-cbind(ErrorStats,sMCE.Ann[,p])
      # colnames(ErrorStats)[ncol(ErrorStats)]<-paste(pt,"sd_ERR_MCMC_Ann",sep="_")
      ErrorStats<-cbind(ErrorStats,mCLT.Burn.d[,p])
      colnames(ErrorStats)[ncol(ErrorStats)]<-paste(pt,"ERR_CLT_Burn",sep="_")
      }
    }
    write.csv(ErrorStats,paste0("ErrorStats_",d,".csv"))
  } else if ((!is.na(fitObject) | !is.na(fitFile)) & !is.na(condition)) {
    source("CreateSimulationStructure.R")
    d<-condition
    simdir<-paste0(gen.dir,"/",d)
    SimList<-readRDS(paste0(simdir,"/",SFileString(sim.list[[d]],gen=TRUE),"_1.rds"))
    gen.xi<-SimList$gen.xi
    fitdir<-paste0(fit.dir,"/",d)
    sims=sim.list[[d]]$Reps
    items=sim.list[[d]]$J
    examinees=sim.list[[d]]$N #,10000)
    XIiter<-array(0, dim=c(items,1+SimList$gen.structure$Adim+(SimList$gen.structure$guess+0),sims))
    if ("tau" %in% names(FitDATA)) {
      TAUiter<-array(0, dim=c(items,FitDATA$settings$ncat-1,sims))
    }
    
    for (i in 1:sims) {
      cat(":",i)
      FitDATA<-readRDS(paste0(fitdir,"/",SFileString(sim.list[[d]],gen=FALSE,r = i),".rds"))
      if (!guessing) {
        FitDATA$settings$guess<-FALSE
      }
      XIiter[,,i]<-as.matrix(FitDATA$xi[,1:(1+FitDATA$settings$Adim+(FitDATA$settings$guess+0))])   # converged
      if ("tau" %in% names(FitDATA)) {
        TAUiter[,,i]<-as.matrix(FitDATA$tau)
      }
    }
    RMSE<-sXI<-apply(XIiter,c(1,2),sd)  # RMSE
    if ("tau" %in% names(FitDATA)) {
      gen.tau<-SimList$gen.tau
      TAU_RMSE<-sTAU<-apply(TAUiter,c(1,2),sd)  # RMSE
    }
    
    if (!is.na(fitFile) & grepl("\\.rda",fitFile)) {
      load(fitFile)
      if (!exists("FitDATA")) {
        print("Don't know if you've got the right rda file: NULL for YOU!")
        return(NULL)
      }
    } else if (!is.na(fitFile) & grepl("\\.rds",fitFile)) {
      FitDATA<-readRDS(fitFile)
    } else if (!is.na(fitFile)) {
      print("Not sure how to read your fitFile... NULL for YOU!")
      return(NULL)
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
                                         settings=FitDATA$settings,item.ind=item.ind))  # CLT during burn
      if ("tau" %in% names(FitDATA) && "Diter" %in% names(FitDATA$Iterations)) {
        mMCE.Burn.d<- as.matrix(mcmcTrimmedError(FitDATA$Iterations,
                                                 end = burninMC[2],start = burninMC[1],
                                                 settings=FitDATA$settings,tau = TRUE))  # Chain during burn
        mCLT.Burn.d<- as.matrix(mcmcCLTError(FitDATA$Iterations,
                                             end = burninCLT[2],start = burninCLT[1],
                                             settings=FitDATA$settings,tau = TRUE,item.ind=item.ind))  # CLT during burn         ::: sd
      }
      
      print("MCMC CLT Burnin NaN Variances:")
      tt<-table(unlist(lapply(apply(MCCLT.Burn,3,function(x) which(is.nan(x),arr.ind = TRUE)),function(y) (y[,1]))))
      if (length(tt)>0) print(tt) else cat("0\n")
    }
    #### Means of the variance estimations
    mXI<-   FitDATA$xi #apply(XIiter,c(1,2),mean)  # converged   :::   parameters
    if ("tau" %in% names(FitDATA)) {
      mTAU<-  apply(TAUiter,c(1,2),mean)  # converged   :::   parameters
    }  
    
    if (!is.na(FitDATA$xiError)[1]) {
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
      #### All SDs
      eLMI<-sqrt(mLMI) # LMI iterative ::: sd
      eiERR<-sqrt(miERR)# LMI at convergence, one calculation - 2pl  ::: sd
      eoERR<-sqrt(moERR) # LMI at convergence, one calculation - 2pno ::: sd
    } else {
      mLMI<- NA 
      miERR<-NA
      moERR<-NA
      eLMI<-NA # LMI iterative ::: sd
      eiERR<-NA# LMI at convergence, one calculation - 2pl  ::: sd
      eoERR<-NA # LMI at convergence, one calculation - 2pno ::: sd
    }
    if (FitDATA$settings$empiricalse) {
      if ("Aiter" %in% names(FitDATA$EmpSE)) {
        mEMC<-as.matrix(mcmcTrimmedError(FitDATA$EmpSE,
                                           end = FitDATA$settings$EmpIT,start = 1,
                                           settings=FitDATA$settings))  # Chain empirical
        mECLT<-as.matrix(mcmcCLTError(FitDATA$EmpSE,
                                       end = FitDATA$settings$EmpIT,start = 1,
                                       settings=FitDATA$settings))   # MCMC CLT empirical
      } else {
        mEMC<-NA
        mECLT<-NA
      }
      if ("VARLMIi" %in% names(FitDATA$EmpSE)) {
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
        # mEMC<- apply(EMCERR,c(1,2),mean)   # Chain empirical    ::: sd
        # mECLT<-apply(EMCCLT,c(1,2),mean)   # MCMC CLT empirical ::: sd
      } else {
        mLMIi<-NA
        mLMIo<-NA
        eLMIi<-NA
        eLMIo<-NA
      }
    }

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
      ErrorStats<-cbind(ErrorStats,RMSE[,p])
      colnames(ErrorStats)[ncol(ErrorStats)]<-paste(pt,"RMSE",sep="_")
      if (!is.na(mLMI)[1] && ncol(mLMI)==ncol(gen.xi)) {
        ErrorStats<-cbind(ErrorStats,mLMI[,p])
        colnames(ErrorStats)[ncol(ErrorStats)]<-paste(pt,"VAR_LMI_Iter",sep="_")
        ErrorStats<-cbind(ErrorStats,eLMI[,p])
        colnames(ErrorStats)[ncol(ErrorStats)]<-paste(pt,"ERR_LMI_Iter",sep="_")
      } else {
        ErrorStats<-cbind(ErrorStats,NA)
        colnames(ErrorStats)[ncol(ErrorStats)]<-paste(pt,"VAR_LMI_Iter",sep="_")
        ErrorStats<-cbind(ErrorStats,NA)
        colnames(ErrorStats)[ncol(ErrorStats)]<-paste(pt,"ERR_LMI_Iter",sep="_")
      }
      ErrorStats<-cbind(ErrorStats,NA)
      colnames(ErrorStats)[ncol(ErrorStats)]<-paste(pt,"sd_ERR_LMI_Iter",sep="_")
      if (!is.na(miERR)[1]) {
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
      } else {
        ErrorStats<-cbind(ErrorStats,NA)
        colnames(ErrorStats)[ncol(ErrorStats)]<-paste(pt,"VAR_Simple_2PL",sep="_")
        ErrorStats<-cbind(ErrorStats,NA)
        colnames(ErrorStats)[ncol(ErrorStats)]<-paste(pt,"ERR_Simple_2PL",sep="_")
        ErrorStats<-cbind(ErrorStats,NA)
        colnames(ErrorStats)[ncol(ErrorStats)]<-paste(pt,"sd_ERR_Simple_2PL",sep="_")
        ErrorStats<-cbind(ErrorStats,NA)
        colnames(ErrorStats)[ncol(ErrorStats)]<-paste(pt,"VAR_Simple_2PNO",sep="_")
        ErrorStats<-cbind(ErrorStats,NA)
        colnames(ErrorStats)[ncol(ErrorStats)]<-paste(pt,"ERR_Simple_2PNO",sep="_")
        ErrorStats<-cbind(ErrorStats,NA)
        colnames(ErrorStats)[ncol(ErrorStats)]<-paste(pt,"sd_ERR_Simple_2PNO",sep="_")
      } 
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
      if ("VARLMIi" %in% names(FitDATA$EmpSE) && FitDATA$settings$empiricalse) {
        ErrorStats<-cbind(ErrorStats,mEMC[,p])
        colnames(ErrorStats)[ncol(ErrorStats)]<-paste(pt,"ERR_MCMC_Emp",sep="_")
        # ErrorStats<-cbind(ErrorStats,sEMC[,p])
        # colnames(ErrorStats)[ncol(ErrorStats)]<-paste(pt,"sd_ERR_MCMC_Emp",sep="_")
        ErrorStats<-cbind(ErrorStats,mECLT[,p])
        colnames(ErrorStats)[ncol(ErrorStats)]<-paste(pt,"ERR_CLT_Emp",sep="_")
        # ErrorStats<-cbind(ErrorStats,sECLT[,p])
        # colnames(ErrorStats)[ncol(ErrorStats)]<-paste(pt,"sd_ERR_CLT_Emp",sep="_")
        if (!is.na(mLMIi)[1] && ncol(mLMIi)==ncol(gen.xi)) {
          ErrorStats<-cbind(ErrorStats,mLMIi[,p])
          colnames(ErrorStats)[ncol(ErrorStats)]<-paste(pt,"VAR_LMI_Iter_Emp_2PL",sep="_")
          ErrorStats<-cbind(ErrorStats,eLMIi[,p])
          colnames(ErrorStats)[ncol(ErrorStats)]<-paste(pt,"ERR_LMI_Iter_Emp_2PL",sep="_")
        } else {
          ErrorStats<-cbind(ErrorStats,NA)
          colnames(ErrorStats)[ncol(ErrorStats)]<-paste(pt,"VAR_LMI_Iter_Emp_2PL",sep="_")
          ErrorStats<-cbind(ErrorStats,NA)
          colnames(ErrorStats)[ncol(ErrorStats)]<-paste(pt,"ERR_LMI_Iter_Emp_2PL",sep="_")
        }
        ErrorStats<-cbind(ErrorStats,NA)
        colnames(ErrorStats)[ncol(ErrorStats)]<-paste(pt,"sd_ERR_LMI_Iter_Emp_2PL",sep="_")
        if (!is.na(mLMIo)[1] && ncol(mLMIo)==ncol(gen.xi)) {
          ErrorStats<-cbind(ErrorStats,mLMIo[,p])
          colnames(ErrorStats)[ncol(ErrorStats)]<-paste(pt,"VAR_LMI_Iter_Emp_2PNO",sep="_")
          ErrorStats<-cbind(ErrorStats,eLMIo[,p])
          colnames(ErrorStats)[ncol(ErrorStats)]<-paste(pt,"ERR_LMI_Iter_Emp_2PNO",sep="_")
        } else {
          ErrorStats<-cbind(ErrorStats,NA)
          colnames(ErrorStats)[ncol(ErrorStats)]<-paste(pt,"VAR_LMI_Iter_Emp_2PNO",sep="_")
          ErrorStats<-cbind(ErrorStats,NA)
          colnames(ErrorStats)[ncol(ErrorStats)]<-paste(pt,"ERR_LMI_Iter_Emp_2PNO",sep="_")
        }
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
    if ("tau" %in% names(FitDATA)) {
      for (p in 1:(FitDATA$settings$ncat-1)) {
        pt<-paste0("Tau",p)
        ErrorStats<-cbind(ErrorStats,gen.tau[,p])
        colnames(ErrorStats)[ncol(ErrorStats)]<-pt
        ErrorStats<-cbind(ErrorStats,mTAU[,p])
        colnames(ErrorStats)[ncol(ErrorStats)]<-paste(pt,"est",sep="_")
        ErrorStats<-cbind(ErrorStats,NA)
        colnames(ErrorStats)[ncol(ErrorStats)]<-paste(pt,"bias",sep="_")
        ErrorStats<-cbind(ErrorStats,sTAU[,p])
        colnames(ErrorStats)[ncol(ErrorStats)]<-paste(pt,"RMSE",sep="_")
        ErrorStats<-cbind(ErrorStats,mMCE.Burn.d[,p])
        colnames(ErrorStats)[ncol(ErrorStats)]<-paste(pt,"ERR_MCMC_Burn",sep="_")
        # ErrorStats<-cbind(ErrorStats,sMCE.Burn[,p])
        # colnames(ErrorStats)[ncol(ErrorStats)]<-paste(pt,"sd_ERR_MCMC_Burn",sep="_")
        # ErrorStats<-cbind(ErrorStats,mMCE.Ann[,p])
        # colnames(ErrorStats)[ncol(ErrorStats)]<-paste(pt,"ERR_MCMC_Ann",sep="_")
        # ErrorStats<-cbind(ErrorStats,sMCE.Ann[,p])
        # colnames(ErrorStats)[ncol(ErrorStats)]<-paste(pt,"sd_ERR_MCMC_Ann",sep="_")
        ErrorStats<-cbind(ErrorStats,mCLT.Burn.d[,p])
        colnames(ErrorStats)[ncol(ErrorStats)]<-paste(pt,"ERR_CLT_Burn",sep="_")
      }
    }
    
    write.csv(ErrorStats,paste0("ErrorStats_",errorFilename,".csv"))
  } else {
    if (!is.na(fitFile) & grepl("\\.rda",fitFile)) {
      load(fitFile)
      if (!exists("FitDATA")) {
        print("Don't know if you've got the right rda file: NULL for YOU!")
        return(NULL)
      }
    } else if (!is.na(fitFile) & grepl("\\.rds",fitFile)) {
      FitDATA<-readRDS(fitFile)
    } else if (!is.na(fitFile)) {
      print("Not sure how to read your fitFile... NULL for YOU!")
      return(NULL)
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
      # mEMC<- apply(EMCERR,c(1,2),mean)   # Chain empirical    ::: sd
      # mECLT<-apply(EMCCLT,c(1,2),mean)   # MCMC CLT empirical ::: sd
      mEMC<-as.matrix(mcmcTrimmedError(FitDATA$EmpSE,
                                         end = FitDATA$settings$EmpIT,start = 1,
                                         settings=FitDATA$settings))  # Chain empirical
      mECLT<-as.matrix(mcmcCLTError(FitDATA$EmpSE,
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
      #ErrorStats<-cbind(ErrorStats,gen.xi[,p])
      #colnames(ErrorStats)[ncol(ErrorStats)]<-pt
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
        # ErrorStats<-cbind(ErrorStats,sEMC[,p])
        # colnames(ErrorStats)[ncol(ErrorStats)]<-paste(pt,"sd_ERR_MCMC_Emp",sep="_")
        ErrorStats<-cbind(ErrorStats,mECLT[,p])
        colnames(ErrorStats)[ncol(ErrorStats)]<-paste(pt,"ERR_CLT_Emp",sep="_")
        # ErrorStats<-cbind(ErrorStats,sECLT[,p])
        # colnames(ErrorStats)[ncol(ErrorStats)]<-paste(pt,"sd_ERR_CLT_Emp",sep="_")
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
  }
  ErrorStats
}