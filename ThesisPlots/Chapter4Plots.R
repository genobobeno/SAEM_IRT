source("InitializeGeisCamilli.R")
basedir<-getwd()
data.dir<-"RealData"
if (!data.dir %in% dir()) dir.create(data.dir)

load("RealData/FCI_NoGuess_A1.rda")
FCI.Fit.1D.NoGuess<-FitDATA
load("RealData/FCI_NoGuess_A2.rda")
FCI.Fit.2D.NoGuess<-FitDATA
load("RealData/FCI_NoGuess_A3.rda")
FCI.Fit.3D.NoGuess<-FitDATA
load("RealData/FCI_NoGuess_A4.rda")
FCI.Fit.4D.NoGuess<-FitDATA
load("RealData/FCI_NoGuess_A10.rda")
FCI.Fit.10D.NoGuess<-FitDATA
rm(FitDATA)
MultipleTWFitTests(fit.data.list = list(FCI.Fit.1D.NoGuess,
                                        FCI.Fit.2D.NoGuess,
                                        FCI.Fit.3D.NoGuess,
                                        FCI.Fit.4D.NoGuess,
                                        FCI.Fit.10D.NoGuess),title = "FCI")

load("RealData/CCI_NoGuess_A2.rda")
CCI.Fit.2D.NoGuess<-FitDATA
load("RealData/FCI_NoGuess_A3.rda")
CCI.Fit.3D.NoGuess<-FitDATA
load("RealData/FCI_NoGuess_A4.rda")
CCI.Fit.4D.NoGuess<-FitDATA
load("RealData/FCI_NoGuess_A5.rda")
CCI.Fit.5D.NoGuess<-FitDATA
load("RealData/FCI_NoGuess_A10.rda")
CCI.Fit.10D.NoGuess<-FitDATA
rm(FitDATA)

MultipleTWFitTests(fit.data.list = list(CCI.Fit.2D.NoGuess,
                                        CCI.Fit.3D.NoGuess,
                                        CCI.Fit.4D.NoGuess,
                                        CCI.Fit.5D.NoGuess,
                                        CCI.Fit.10D.NoGuess),title = "CCI")

