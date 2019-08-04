CalculateMCMCErrors<-function(Achain,Bchain,settings,Cchain=NA,
                              Dchain=NA,item.ind=TRUE,start=100,end=NA,
                              condition=NA) {
  if (!is.na(condition)) {
    l.RMSE<-GetRMSE(condition)
    RMSE<-l.RMSE$RMSE
    gen.xi<-l.RMSE$gen.xi
    settings<-l.RMSE$settings
  } else {
    gen.xi<-RMSE<-NA
  }
  if (is.na(Cchain)[1]) settings$guess<-FALSE
  if (is.na(end)[1]) end<-dim(Bchain)[2]
  if (!is.na(Dchain)[1]) {
    MClist<-list(Aiter=Achain,Biter=Bchain,Diter=Dchain)
  } else if (!is.na(Cchain)[1]) {
    MClist<-list(Aiter=Achain,Biter=Bchain,Citer=Cchain)
  } else {
    MClist<-list(Aiter=Achain,Biter=Bchain)
  }
  MCERR.Burn<-as.matrix(mcmcTrimmedError(MClist,
                                              end = end,
                                              start = start,
                                            settings=settings))   # Chain during burn
  MCCLT.Burn<-as.matrix(mcmcCLTError(MClist,
                                          end = end,
                                          start = start,
                                          settings=settings,item.ind=item.ind))  # CLT during burn
  if (settings$guess) {
    colnames(MCERR.Burn)<-c(paste0("A",1:settings$Adim,"_MCERR"),"B_MCERR","C_MCERR")
    colnames(MCCLT.Burn)<-c(paste0("A",1:settings$Adim,"_CLTERR"),
                            "B_CLTERR","C_CLTERR")
  } else {
    colnames(MCERR.Burn)<-c(paste0("A",1:settings$Adim,"_MCERR"),"B_MCERR")
    colnames(MCCLT.Burn)<-c(paste0("A",1:settings$Adim,"_CLTERR"),"B_CLTERR")
  }
  if (!is.na(Dchain)[1])  {
    MCERR.Burn.D<-as.matrix(mcmcTrimmedError(MClist,
                                           end = end,
                                           start = start,
                                           settings=settings,tau = TRUE))   # Chain during burn
    colnames(MCERR.Burn.D)<-paste0("D",1:(settings$ncat-1),"_MCERR")
    MCCLT.Burn.D<-as.matrix(mcmcCLTError(MClist,
                                       end = end,
                                       start = start,
                                       settings=settings,item.ind=item.ind,tau=TRUE))  # CLT during burn
    colnames(MCCLT.Burn.D)<-paste0("D",1:(settings$ncat-1),"_CLTERR")
  }
  if (!is.na(RMSE)[1]) {
    BayesErrors<-list()
    for (q in 1:settings$Adim) {# q=1;j=1
      BayesErrors[[paste0("A",q)]]<-sapply(1:nrow(Bchain),function(j) sum(Achain[j,q,start:end]>(gen.xi[j,q]-2*RMSE[j,paste0("A",q,"_RMSE")]) & 
                                                 Achain[j,q,start:end]<(gen.xi[j,q]+2*RMSE[j,paste0("A",q,"_RMSE")]))/(1+end-start))
      print(BayesErrors[[paste0("A",q)]])
    }
    BayesErrors[["B"]]<-sapply(1:nrow(Bchain),function(j) sum(Bchain[j,start:end]>(gen.xi[j,ncol(gen.xi)-settings$guess]-2*RMSE[j,"B_RMSE"]) & 
                                                             Bchain[j,start:end]<(gen.xi[j,ncol(gen.xi)-settings$guess]+2*RMSE[j,"B_RMSE"]))/(1+end-start))
    print(BayesErrors[["B"]])
    if (settings$guess) {
      BayesErrors[["C"]]<-sapply(1:nrow(Bchain),function(j) sum(Cchain[j,start:end]>(gen.xi[j,ncol(gen.xi)]-2*RMSE[j,"C_RMSE"]) & 
                                                               Cchain[j,start:end]<(gen.xi[j,ncol(gen.xi)]+2*RMSE[j,"C_RMSE"]))/(1+end-start))
      print(BayesErrors[["C"]])
    }
    if (!is.na(settings$ncat) && settings$ncat>2) {
      for (k in 1:(settings$ncat-1)) {
        gen.tau<-l.RMSE$gen.tau
        b<-gen.xi[,ncol(gen.xi)-settings$guess]
        BayesErrors[[paste0("D",k)]]<-sapply(1:nrow(Bchain),function(j) sum(Dchain[j,k,start:end]+b[j]>(gen.tau[j,k]-2*RMSE[j,paste0("Tau",k,"_RMSE")]) & 
                                                                           Dchain[j,k,start:end]+b[j]<(gen.tau[j,k]+2*RMSE[j,paste0("Tau",k,"_RMSE")]))/(1+end-start))
        print(BayesErrors[[paste0("D",k)]])
      }
      list(MCERR=cbind(RMSE,MCERR.Burn,MCERR.Burn.D),
           MCCLT=cbind(RMSE,MCCLT.Burn,MCCLT.Burn.D),BayesError=BayesErrors)
    } else {
      list(MCERR=cbind(RMSE,MCERR.Burn),
           MCCLT=cbind(RMSE,MCCLT.Burn),BayesError=BayesErrors)
    }
  } else {
    list(MCERR=cbind(MCERR.Burn,MCERR.Burn.D),
         MCCLT=cbind(MCCLT.Burn,MCCLT.Burn.D))  #,BayesError=BayesErrors Need a definition of different Bayes for Real Data
  }
}