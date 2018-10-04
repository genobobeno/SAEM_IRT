setwd("~/ParSAEM/SAEM_IRT")
source("CreateSimulationStructure.R")

d = "S3"
r = 10

simdir<-paste0(gen.dir,"/",d)
SimList<-readRDS(paste0(simdir,"/",SFileString(sim.list[[d]],gen=TRUE),"_",r,".rds"))
print(SimList$gen.list$gen.xi)

fitdir<-paste0(fit.dir,"/",d)
FitList<-readRDS(paste0(fitdir,"/",SFileString(sim.list[[d]],gen=FALSE,r = r),".rds"))
# if (d=="S2") {
#   print(FitList$XI - FitList$xi)
#   plot(density(FitList$XI[,1] - FitList$xi[,1]))
#   for ( i in 2:ncol(FitList$XI))) lines(density(FitList$XI[,i] - FitList$xi[,i]),col=i)
# } else {
  print(FitList$XI - FitList$xi)
  plot(density(FitList$XI[,1] - FitList$xi[,1]))
  for ( i in 2:ncol(FitList$XI)) lines(density(FitList$XI[,i] - FitList$xi[,i]),col=i)
  if (!is.na(FitList$TAU[1])) {
    print(FitList$TAU - FitList$tau)
    plot(density(as.vector(FitList$TAU) - as.vector(FitList$tau)))
  }


plot(FitList$Iterations$Ait[1,1,],type="n")
lines(FitList$Iterations$Ait[1,1,])
acf(FitList$Iterations$Bit[1,])
eigen(FitList$EZZ)
