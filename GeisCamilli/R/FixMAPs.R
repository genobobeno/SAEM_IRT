FixMAPs<-function(fitFile) {
  FitList<-readRDS(fitFile)
  Tmap<-ThetaMAP(aa=FitList$A,bb=FitList$B,cc=FitList$C,rp=FitList$RP,
                 settings=FitList$settings,TAU=FitList$tau,That=FitList$That)
  FitList$Tmap<-Tmap
  saveRDS(FitList,fitFile)
}
