FixMAPs<-function(fitFile) {
  FitList<-readRDS(fitFile)
  if ("tau" %in% names(FitList)) {
    Tmap<-ThetaMAP(aa=FitList$A,bb=FitList$B,cc=FitList$C,rp=FitList$RP,
                   settings=FitList$settings,TAU=FitList$tau,That=FitList$That)
  } else {
    Tmap<-ThetaMAP(aa=FitList$A,bb=FitList$B,cc=FitList$C,rp=FitList$RP,
                   settings=FitList$settings,TAU=NA,That=FitList$That)
  }
    FitList$Tmap<-Tmap
  saveRDS(FitList,fitFile)
}
