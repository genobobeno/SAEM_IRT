WorkingDirectory = "."
setwd(WorkingDirectory)
PARALLEL=FALSE
options(digits = 5)
# if (!exists("jj")) {
#  jj=0
# }
#"ConvergedModelFits/"
library(mcmc)
AutocovMC<-function(x,burnin=800,start=600,items=1:5,title.append=NA) {
  par(mfrow=c(5,2),mar=c(3,3,3,1))
  if (is.na(items)[1]) items<-dim(x$Aiter)[1]
  for (j in items) {
    acf(x$Aiter[j,1,start:burnin],main=paste0(ifelse(is.na(title.append),"",title.append),"Slope, Item ",j),xlab=NA)
    acf(x$Biter[j,start:burnin],main=paste0(ifelse(is.na(title.append),"",title.append),"Intercept, Item ",j),xlab=NA)
  }
}
AutocovMC(MCMCDATA,burnin = burnin[2],start = burnin[1],title.append = "Iteration 600-800: ")   # Chain during burn
AutocovMC(MCMCDATA,burnin = burnin[3],start = burnin[2],title.append = "Iteration 800-1000: ")   # Chain during Annealing window

MCErrorFunction<-function(x,burnin=800,start=600,trimA=8,trimB=5) {
  for (iii in 1:dim(x$Aiter)[2]) {
    Atrim = start+trimA*1:(floor((burnin-start)/trimA))
    if (iii==1) {
      AE<-apply(x$Aiter[,1,Atrim],1,sd)
    } else {
      AE<-cbind(AE,apply(x$Aiter[,1,Atrim],1,sd))
    }
  }
  Btrim = start+trimB*1:(floor((burnin-start)/trimB))
  BE<-apply(x$Biter[,Btrim],1,sd)
  cbind(AE,BE)
}

CLTErrorFunction<-function(x,burnin=800,start=600,trimA=8,trimB=5) {
  for (iii in 1:dim(x$Aiter)[2]) {
    Atrim = start+trimA*1:(floor((burnin-start)/trimA))
    if (iii==1) {
      AE<-apply(x$Aiter[,1,Atrim],1,function(y) (sqrt(initseq(y)$var.pos)))
    } else {
      AE<-cbind(AE,apply(x$Aiter[,1,Atrim],1,function(y) (sqrt(initseq(y)$var.pos))))
    }
  }
  Btrim = start+trimB*1:(floor((burnin-start)/trimB))
  BE<-apply(x$Biter[,Btrim],1,function(y) (sqrt(initseq(y)$var.pos)))
  cbind(AE,BE)
}

negsqrt<-function(x) {
  if (x<0) return(0)
  else return(sqrt(x))
}


###############
# RUN 1D Test #
###############
source("CreateSimulationStructure.R")
d<-condition<-"S1"
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
FitList<-readRDS(paste0(fitdir,"/",SFileString(sim.list[[d]],gen=FALSE,r = SIM),".rds"))
load(paste0(fitdir,"/",SFileString(sim.list[[d]],gen=FALSE,r = SIM),".rda"))
burnin = c(ceiling(0.6*FitList$settings$burnin),ceiling(0.8*FitList$settings$burnin),FitList$settings$burnin)
burn = burnin[2]

# for (b in 1:length(burn)) {
#   simfile=paste("Error",burn[b],"/SimError1D-",items,"-",examinees,"-s",SIM,".rda",sep="")    # Save to an ".rda" file with this name
#   load(simfile)
#   estfile=paste("Error",burn[b],"/FitError1D-",items,"-",examinees,"-s",SIM,"-r1.rda",sep="") 
#       load(estfile)
XIiter<-array(0, dim=c(items,1+FitDATA$settings$Adim,sims))        
ERRiter<-array(0, dim=c(items,1+FitDATA$settings$Adim,sims)) # LMI iterative - 2PNO
iERR<-array(0, dim=c(items,1+FitDATA$settings$Adim,sims))    # LMI simple - 2PL
oERR<-array(0, dim=c(items,1+FitDATA$settings$Adim,sims))    # LMI simple - 2PNO
MCERR.Burn<-array(0, dim=c(items,1+FitDATA$settings$Adim,sims))  # Chain during burn
MCERR.Ann<-array(0, dim=c(items,1+FitDATA$settings$Adim,sims))   # Chain during Ann window
MCCLT.Burn<-array(0, dim=c(items,1+FitDATA$settings$Adim,sims))  # CLT during burn
MCCLT.Ann<-array(0, dim=c(items,1+FitDATA$settings$Adim,sims))   # CLT during Ann window
EMCERR<-array(0, dim=c(items,1+FitDATA$settings$Adim,sims))  # Chain empirical
EMCCLT<-array(0, dim=c(items,1+FitDATA$settings$Adim,sims))   # MCMC CLT empirical
EmpLMIi<-array(0, dim=c(items,1+FitDATA$settings$Adim,sims)) # LMI iterative - 2PL
EmpLMIo<-array(0, dim=c(items,1+FitDATA$settings$Adim,sims)) # LMI iterative - 2PNO
for (i in 1:sims) {
  estfile=paste0(fitdir,"/",SFileString(sim.list[[d]],gen=FALSE,r = i),".rda")
  load(estfile)
  XIiter[,,i]<-as.matrix(FitDATA$xi)   # converged
  ERRiter[,,i]<-as.matrix(FitDATA$xiError)   # LMI iterative
  iERR[,,i]<-as.matrix(FitDATA$iError)   # LMI at convergence, one calculation - 2pl
  oERR[,,i]<-as.matrix(FitDATA$oError) # LMI at convergence, one calculation - 2pno
  MCERR.Burn[,,i]<-as.matrix(MCErrorFunction(MCMCDATA,burnin = burnin[2],start = burnin[1]))   # Chain during burn
  MCERR.Ann[,,i]<-as.matrix(MCErrorFunction(MCMCDATA,burnin = burnin[3],start = burnin[2]))    # Chain during Ann window
  MCCLT.Burn[,,i]<-as.matrix(CLTErrorFunction(MCMCDATA,burnin = burnin[2],start = burnin[1]))  # CLT during burn
  MCCLT.Ann[,,i]<-as.matrix(CLTErrorFunction(MCMCDATA,burnin = burnin[3],start = burnin[2]))   # CLT during Ann window
  EMCERR[,,i]<-as.matrix(cbind(FitDATA$EmpSE$MCSA,FitDATA$EmpSE$MCSB)) # Chain empirical
  EMCCLT[,,i]<-as.matrix(cbind(FitDATA$EmpSE$SEA,FitDATA$EmpSE$SEB))   # MCMC CLT empirical
  EmpLMIi[,,i]<-as.matrix(FitDATA$EmpSE$VARLMIi)  # LMI empirical - iterative - 2pl
  EmpLMIo[,,i]<-as.matrix(FitDATA$EmpSE$VARLMIo)  # LMI empirical - iterative - 2pno
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
print("MCMC CLT Annealing NaN Variances:")
tt<-table(unlist(lapply(apply(MCCLT.Ann,3,function(x) which(is.nan(x),arr.ind = TRUE)),function(y) (y[,1]))))
if (length(tt)>0) print(tt) else cat("0\n")
print("Empirical MCMC CLT Burnin NaN Variances:")
tt<-table(unlist(lapply(apply(EMCCLT,3,function(x) which(is.nan(x),arr.ind = TRUE)),function(y) (y[,1]))))
if (length(tt)>0) print(tt) else cat("0\n")


#### All SDs
sXI<-apply(XIiter,c(1,2),sd)  # RMSE
eLMI<-matrix(sqrt(mLMI),nrow=items,ncol=1+FitDATA$settings$Adim)  # LMI iterative ::: sd
eiERR<-matrix(sqrt(miERR),nrow=items,ncol=1+FitDATA$settings$Adim) # LMI at convergence, one calculation - 2pl  ::: sd
eoERR<-matrix(sqrt(moERR),nrow=items,ncol=1+FitDATA$settings$Adim) # LMI at convergence, one calculation - 2pno ::: sd
mMCE.Burn<- apply(MCERR.Burn,c(1,2),mean)   # Chain during burn         ::: sd
mMCE.Ann<-  apply(MCERR.Ann,c(1,2),mean)    # Chain during Ann window   ::: sd
mCLT.Burn<- apply(MCCLT.Burn,c(1,2),function(x) mean(x[!is.nan(x)]))   # CLT during burn           ::: sd
mCLT.Ann<-  apply(MCCLT.Ann,c(1,2),function(x) mean(x[!is.nan(x)]))    # CLT during Ann window     ::: sd
mEMC<- apply(EMCERR,c(1,2),mean)   # Chain empirical    ::: sd
mECLT<-apply(EMCCLT,c(1,2),mean)   # MCMC CLT empirical ::: sd
eLMIi<-matrix(sqrt(mLMIi),nrow=items,ncol=1+FitDATA$settings$Adim) # LMI empirical - iterative - 2pl  ::: sd
eLMIo<-matrix(sqrt(mLMIo),nrow=items,ncol=1+FitDATA$settings$Adim) # LMI empirical - iterative - 2pno ::: sd


gen.xi<-SimList$gen.xi
bias<-mXI-gen.xi
#######  SD of the variances
sLMI <-apply(ERRiter,c(1,2),function(x) sd(sqrt(x[x>0])))
siERR<-apply(iERR,c(1,2),function(x) sd(sqrt(x[x>0])))
soERR<-apply(oERR,c(1,2),function(x) sd(sqrt(x[x>0])))
sMCE.Burn<-apply(MCERR.Burn,c(1,2),sd)
sMCE.Ann <-apply(MCERR.Ann,c(1,2),sd)
sCLT.Burn<-apply(MCCLT.Burn,c(1,2),sd)
sCLT.Ann <-apply(MCCLT.Ann,c(1,2),sd)
sEMC <-apply(EMCERR,c(1,2),sd)
sECLT<-apply(EMCCLT,c(1,2),sd)
sLMIi<-apply(EmpLMIi,c(1,2),function(x) sd(sqrt(x[x>0])))
sLMIo<-apply(EmpLMIo,c(1,2),function(x) sd(sqrt(x[x>0])))

stats1D<-data.frame(N=rep(examinees,items),J=rep(items,items),Burn=rep(burn,items))
for(p in 1:ncol(gen.xi)) {
  ifelse(p==ncol(gen.xi),pt<-"B",pt<-paste("A",p,sep=""))
  stats1D<-cbind(stats1D,gen.xi[,p])
  colnames(stats1D)[ncol(stats1D)]<-pt
  stats1D<-cbind(stats1D,mXI[,p])
  colnames(stats1D)[ncol(stats1D)]<-paste(pt,"est",sep="_")
  stats1D<-cbind(stats1D,bias[,p])
  colnames(stats1D)[ncol(stats1D)]<-paste(pt,"bias",sep="_")
  stats1D<-cbind(stats1D,sXI[,p])
  colnames(stats1D)[ncol(stats1D)]<-paste(pt,"RMSE",sep="_")
  stats1D<-cbind(stats1D,mLMI[,p])
  colnames(stats1D)[ncol(stats1D)]<-paste(pt,"VAR_LMI_Iter",sep="_")
  stats1D<-cbind(stats1D,eLMI[,p])
  colnames(stats1D)[ncol(stats1D)]<-paste(pt,"ERR_LMI_Iter",sep="_")
  stats1D<-cbind(stats1D,sLMI[,p])
  colnames(stats1D)[ncol(stats1D)]<-paste(pt,"sd_ERR_LMI_Iter",sep="_")
  stats1D<-cbind(stats1D,miERR[,p])
  colnames(stats1D)[ncol(stats1D)]<-paste(pt,"VAR_Simple_2PL",sep="_")
  stats1D<-cbind(stats1D,eiERR[,p])
  colnames(stats1D)[ncol(stats1D)]<-paste(pt,"ERR_Simple_2PL",sep="_")
  stats1D<-cbind(stats1D,siERR[,p])
  colnames(stats1D)[ncol(stats1D)]<-paste(pt,"sd_ERR_Simple_2PL",sep="_")
  stats1D<-cbind(stats1D,moERR[,p])
  colnames(stats1D)[ncol(stats1D)]<-paste(pt,"VAR_Simple_2PNO",sep="_")
  stats1D<-cbind(stats1D,eoERR[,p])
  colnames(stats1D)[ncol(stats1D)]<-paste(pt,"ERR_Simple_2PNO",sep="_")
  stats1D<-cbind(stats1D,soERR[,p])
  colnames(stats1D)[ncol(stats1D)]<-paste(pt,"sd_ERR_Simple_2PNO",sep="_")
  stats1D<-cbind(stats1D,mMCE.Burn[,p])
  colnames(stats1D)[ncol(stats1D)]<-paste(pt,"ERR_MCMC_Burn",sep="_")
  stats1D<-cbind(stats1D,sMCE.Burn[,p])
  colnames(stats1D)[ncol(stats1D)]<-paste(pt,"sd_ERR_MCMC_Burn",sep="_")
  stats1D<-cbind(stats1D,mMCE.Ann[,p])
  colnames(stats1D)[ncol(stats1D)]<-paste(pt,"ERR_MCMC_Ann",sep="_")
  stats1D<-cbind(stats1D,sMCE.Ann[,p])
  colnames(stats1D)[ncol(stats1D)]<-paste(pt,"sd_ERR_MCMC_Ann",sep="_")
  stats1D<-cbind(stats1D,mCLT.Burn[,p])
  colnames(stats1D)[ncol(stats1D)]<-paste(pt,"ERR_CLT_Burn",sep="_")
  stats1D<-cbind(stats1D,sCLT.Burn[,p])
  colnames(stats1D)[ncol(stats1D)]<-paste(pt,"sd_ERR_CLT_Burn",sep="_")
  stats1D<-cbind(stats1D,mCLT.Ann[,p])
  colnames(stats1D)[ncol(stats1D)]<-paste(pt,"ERR_CLT_Ann",sep="_")
  stats1D<-cbind(stats1D,sCLT.Ann[,p])
  colnames(stats1D)[ncol(stats1D)]<-paste(pt,"sd_ERR_CLT_Ann",sep="_")
  stats1D<-cbind(stats1D,mEMC[,p])
  colnames(stats1D)[ncol(stats1D)]<-paste(pt,"ERR_MCMC_Emp",sep="_")
  stats1D<-cbind(stats1D,sEMC[,p])
  colnames(stats1D)[ncol(stats1D)]<-paste(pt,"sd_ERR_MCMC_Emp",sep="_")
  stats1D<-cbind(stats1D,mECLT[,p])
  colnames(stats1D)[ncol(stats1D)]<-paste(pt,"ERR_CLT_Emp",sep="_")
  stats1D<-cbind(stats1D,sECLT[,p])
  colnames(stats1D)[ncol(stats1D)]<-paste(pt,"sd_ERR_CLT_Emp",sep="_")
  stats1D<-cbind(stats1D,mLMIi[,p])
  colnames(stats1D)[ncol(stats1D)]<-paste(pt,"VAR_LMI_Iter_Emp_2PL",sep="_")
  stats1D<-cbind(stats1D,eLMIi[,p])
  colnames(stats1D)[ncol(stats1D)]<-paste(pt,"ERR_LMI_Iter_Emp_2PL",sep="_")
  stats1D<-cbind(stats1D,sLMIi[,p])
  colnames(stats1D)[ncol(stats1D)]<-paste(pt,"sd_ERR_LMI_Iter_Emp_2PL",sep="_")
  stats1D<-cbind(stats1D,mLMIo[,p])
  colnames(stats1D)[ncol(stats1D)]<-paste(pt,"VAR_LMI_Iter_Emp_2PNO",sep="_")
  stats1D<-cbind(stats1D,eLMIo[,p])
  colnames(stats1D)[ncol(stats1D)]<-paste(pt,"ERR_LMI_Iter_Emp_2PNO",sep="_")
  stats1D<-cbind(stats1D,sLMIo[,p])
  colnames(stats1D)[ncol(stats1D)]<-paste(pt,"sd_ERR_LMI_Iter_Emp_2PNO",sep="_")
} # else {
#   st<-cbind(rep(examinees,items),rep(items,items),rep(burn[b],items))
#   for(p in 1:ncol(gen.xi)) {
#     st<-cbind(st,gen.xi[,p],bias[,p],sXI[,p],mLMI[,p],sqrt(mLMI[,p]),sLMI[,p],miERR[,p],
#               sqrt(miERR[,p]),siERR[,p],moERR[,p],sqrt(moERR[,p]),soERR[,p],mMCE[,p],sMCE[,p],
#               mLMIi[,p],eLMIi[,p],sLMIi[,p],mLMIo[,p],eLMIo[,p],sLMIo[,p])
#   }
#   colnames(st)<-colnames(stats1D)
#   stats1D<-rbind(stats1D,st)
# }

write.csv(stats1D,"Stats1D.csv")


###### PLOTS
#nf <- layout(matrix(c(1,2,3,4,6,7,8,9,5,5,5,5,10,10,10,10),8,2), c(16,16), c(4,3,3,5,4,3,3,5), TRUE)
nf <- layout(matrix(c(1,2,3,4,5,6,7,8,9,10,11,12,12,12,12,12,13,13,13,13,13,13),11,2), c(14,20), 
             c(4,3,3,3,3,3,3,3,3,3,5,4,3,3,3,3,3,3,3,3,3,5), TRUE)
#layout.show(nf)

YR<-range(c(stats1D$A1_ERR_CLT_Ann,
            stats1D$A1_ERR_CLT_Burn,
            stats1D$A1_ERR_CLT_Emp,
            stats1D$A1_ERR_LMI_Iter,
            stats1D$A1_ERR_LMI_Iter_Emp_2PL,
            stats1D$A1_ERR_LMI_Iter_Emp_2PNO,
            stats1D$A1_ERR_Simple_2PL,
            stats1D$A1_ERR_Simple_2PNO,
            stats1D$A1_ERR_MCMC_Ann,
            stats1D$A1_ERR_MCMC_Burn,
            stats1D$A1_ERR_MCMC_Emp)/stats1D$A1_RMSE)

txt.labels<-c(expression("CLT"["Ann"]),
              expression("CLT"["Burn"]),
              expression("CLT"["Emp"]),
              expression("MCMC"["Ann"]),
              expression("MCMC"["Burn"]),
              expression("MCMC"["Emp"]),
              expression("ICE"),
              expression("IPCE"["2PL"]),
              expression("IPCE"["2PNO"]),
              expression("SPCE"["2PL"]),
              expression("SPCE"["2PNO"]))

ticks <- pretty(YR)[-c(length(pretty(YR))-1:0)]
labels <- format(ticks, big.mark=",", scientific=FALSE)
m.b = 0.25
par(mar=c(m.b,4,2,m.b))
# for (p in 1:ncol(gen.xi)) {
p=1
ifelse(p==ncol(gen.xi),pt<-"B",pt<-paste("A",p,sep=""))
plot(gen.xi[,p],rep(c(YR),sim.list[[d]]$J/2),type="n",
     main=expression("Slope Error Estimates : 2PNO"),ylim=YR,
     xlab=NA,ylab=NA,xaxt="n",yaxt="n")
axis(2,at=ticks,labels=labels,las=2)
abline(h=1,lwd=2,lty=2)
points(gen.xi[,p],stats1D$A1_ERR_CLT_Ann/stats1D$A1_RMSE,pch=0,col=1)
text(0.6,0.9*YR[2],txt.labels[1])

par(mar=c(m.b,4,m.b,m.b))
plot(gen.xi[,p],rep(c(YR),sim.list[[d]]$J/2),type="n",
     main=NA,ylim=YR,xlab=NA,ylab=NA,xaxt="n",yaxt="n")
axis(2,at=ticks,labels=labels,las=2)
abline(h=1,lwd=2,lty=2)
points(gen.xi[,p],stats1D$A1_ERR_CLT_Burn/stats1D$A1_RMSE,pch=1,col=1)
text(0.6,0.9*YR[2],txt.labels[2])

plot(gen.xi[,p],rep(c(YR),sim.list[[d]]$J/2),type="n",
     main=NA,ylim=YR,xlab=NA,ylab=NA,xaxt="n",yaxt="n")
axis(2,at=ticks,labels=labels,las=2)
abline(h=1,lwd=2,lty=2)
points(gen.xi[,p],stats1D$A1_ERR_CLT_Emp/stats1D$A1_RMSE,pch=2,col=1)
text(0.6,0.9*YR[2],txt.labels[3])

plot(gen.xi[,p],rep(c(YR),sim.list[[d]]$J/2),type="n",
     main=NA,ylim=YR,xlab=NA,ylab=NA,xaxt="n",yaxt="n")
axis(2,at=ticks,labels=labels,las=2)
abline(h=1,lwd=2,lty=2)
points(gen.xi[,p],stats1D$A1_ERR_MCMC_Ann/stats1D$A1_RMSE,pch=8,col=1)
text(0.6,0.9*YR[2],txt.labels[4])

plot(gen.xi[,p],rep(c(YR),sim.list[[d]]$J/2),type="n",
     main=NA,ylim=YR,xlab=NA,ylab=NA,xaxt="n",yaxt="n")
axis(2,at=ticks,labels=labels,las=2)
abline(h=1,lwd=2,lty=2)
points(gen.xi[,p],stats1D$A1_ERR_MCMC_Burn/stats1D$A1_RMSE,pch=9,col=1)
text(0.6,0.9*YR[2],txt.labels[5])

plot(gen.xi[,p],rep(c(YR),sim.list[[d]]$J/2),type="n",
     main=NA,ylim=YR,xlab=NA,ylab=NA,yaxt="n",xaxt="n")
axis(2,at=ticks,labels=labels,las=2)
abline(h=1,lwd=2,lty=2)
points(gen.xi[,p],stats1D$A1_ERR_MCMC_Emp/stats1D$A1_RMSE,pch=10,col=1)
text(0.6,0.9*YR[2],txt.labels[6])
mtext(text = "Ratio of Error Approximation to RMSE",side = 2,line = 2.5,cex=0.7,adj=0)

plot(gen.xi[,p],rep(c(YR),sim.list[[d]]$J/2),type="n",
     main=NA,ylim=YR,xlab=NA,ylab=NA,xaxt="n",yaxt="n")
axis(2,at=ticks,labels=labels,las=2)
abline(h=1,lwd=2,lty=2)
points(gen.xi[,p],stats1D$A1_ERR_LMI_Iter/stats1D$A1_RMSE,pch=3,col=1)
text(0.6,0.9*YR[2],txt.labels[7])

plot(gen.xi[,p],rep(c(YR),sim.list[[d]]$J/2),type="n",
     main=NA,ylim=YR,xlab=NA,ylab=NA,xaxt="n",yaxt="n")
axis(2,at=ticks,labels=labels,las=2)
abline(h=1,lwd=2,lty=2)
points(gen.xi[,p],stats1D$A1_ERR_LMI_Iter_Emp_2PL/stats1D$A1_RMSE,pch=4,col=1)
text(0.6,0.9*YR[2],txt.labels[8])

plot(gen.xi[,p],rep(c(YR),sim.list[[d]]$J/2),type="n",
     main=NA,ylim=YR,xlab=NA,ylab=NA,xaxt="n",yaxt="n")
axis(2,at=ticks,labels=labels,las=2)
abline(h=1,lwd=2,lty=2)
points(gen.xi[,p],stats1D$A1_ERR_LMI_Iter_Emp_2PNO/stats1D$A1_RMSE,pch=5,col=1)
text(0.6,0.9*YR[2],txt.labels[9])

plot(gen.xi[,p],rep(c(YR),sim.list[[d]]$J/2),type="n",
     main=NA,ylim=YR,xlab=NA,ylab=NA,xaxt="n",yaxt="n")
axis(2,at=ticks,labels=labels,las=2)
abline(h=1,lwd=2,lty=2)
points(gen.xi[,p],stats1D$A1_ERR_Simple_2PL/stats1D$A1_RMSE,pch=6,col=1)
text(0.6,0.9*YR[2],txt.labels[10])

par(mar=c(4,4,m.b,m.b))
plot(gen.xi[,p],rep(c(YR),sim.list[[d]]$J/2),type="n",
     main=NA,ylim=YR,xlab=NA,ylab=NA,yaxt="n")
axis(2,at=ticks,labels=labels,las=2)
abline(h=1,lwd=2,lty=2)
points(gen.xi[,p],stats1D$A1_ERR_Simple_2PNO/stats1D$A1_RMSE,pch=7,col=1)
text(0.6,0.9*YR[2],txt.labels[11])


fitA1.CLT_Ann<-loess(A1_ERR_CLT_Ann/A1_RMSE~A1,data=stats1D,span=0.5)
fitA1.CLT_Ann
sd(fitA1.CLT_Ann$residuals)
#fitA1o<-lm(A1_sd_LMIo/A1_sd~A1,data=stats1D)
fitA1.CLT_Burn<-loess(A1_ERR_CLT_Burn/A1_RMSE~A1,data=stats1D,span=0.5)
fitA1.CLT_Burn
sd(fitA1.CLT_Burn$residuals)
#fitA1MC<-lm(A1_sd_MCE/A1_sd~A1,data=stats1D)
fitA1.CLT_Emp<-loess(A1_ERR_CLT_Emp/A1_RMSE~A1,data=stats1D,span=0.5)
fitA1.CLT_Emp
sd(fitA1.CLT_Emp$residuals)
#fitA1LM<-lm(A1_sd_LMI/A1_sd~A1,data=stats1D)
fitA1.LMI_Iter<-loess(A1_ERR_LMI_Iter/A1_RMSE~A1,data=stats1D,span=0.5)
fitA1.LMI_Iter
sd(fitA1.LMI_Iter$residuals)
fitA1.LMI_Iter_Emp_2PL<-loess(A1_ERR_LMI_Iter_Emp_2PL/A1_RMSE~A1,data=stats1D,span=0.5)
fitA1.LMI_Iter_Emp_2PL
sd(fitA1.LMI_Iter_Emp_2PL$residuals)
#fitA1o<-lm(A1_sd_LMIo/A1_sd~A1,data=stats1D)
fitA1.LMI_Iter_Emp_2PNO<-loess(A1_ERR_LMI_Iter_Emp_2PNO/A1_RMSE~A1,data=stats1D,span=0.5)
fitA1.LMI_Iter_Emp_2PNO
sd(fitA1.LMI_Iter_Emp_2PNO$residuals)
#fitA1MC<-lm(A1_sd_MCE/A1_sd~A1,data=stats1D)
fitA1.Simple_2PL<-loess(A1_ERR_Simple_2PL/A1_RMSE~A1,data=stats1D,span=0.5)
fitA1.Simple_2PL
sd(fitA1.Simple_2PL$residuals)
#fitA1LM<-lm(A1_sd_LMI/A1_sd~A1,data=stats1D)
fitA1.Simple_2PNO<-loess(A1_ERR_Simple_2PNO/A1_RMSE~A1,data=stats1D,span=0.5)
fitA1.Simple_2PNO
sd(fitA1.Simple_2PNO$residuals)
fitA1.MCMC_Ann<-loess(A1_ERR_MCMC_Ann/A1_RMSE~A1,data=stats1D,span=0.5)
fitA1.MCMC_Ann
sd(fitA1.MCMC_Ann$residuals)
#fitA1o<-lm(A1_sd_LMIo/A1_sd~A1,data=stats1D)
fitA1.MCMC_Burn<-loess(A1_ERR_MCMC_Burn/A1_RMSE~A1,data=stats1D,span=0.5)
fitA1.MCMC_Burn
sd(fitA1.MCMC_Burn$residuals)
#fitA1MC<-lm(A1_sd_MCE/A1_sd~A1,data=stats1D)
fitA1.MCMC_Emp<-loess(A1_ERR_MCMC_Emp/A1_RMSE~A1,data=stats1D,span=0.5)
fitA1.MCMC_Emp
sd(fitA1.MCMC_Emp$residuals)
######plot2
par(mar=c(4,6,4,m.b))
seqA1<-seq(range(stats1D$A1)[1],range(stats1D$A1)[2],length.out = 50)
plot(stats1D$A1,stats1D$A1_ERR_MCMC_Emp/stats1D$A1_RMSE,main=expression("Loess fit to"~hat(sigma)/sigma["RMSE"]),
     ylab="Ratio of Error Approximation to RMSE",xlab="Generated Slope (1D, 2PNO)",ylim=YR,type="n")
lines(seqA1,predict(fitA1.CLT_Ann,newdata = data.frame(A1=seqA1)),lty=1)
lines(seqA1,predict(fitA1.CLT_Burn,newdata = data.frame(A1=seqA1)),lty=2,col=2)
lines(seqA1,predict(fitA1.CLT_Emp,newdata = data.frame(A1=seqA1)),lty=3,col=3)
lines(seqA1,predict(fitA1.MCMC_Ann,newdata = data.frame(A1=seqA1)),lty=4,col=4,lwd=2)
lines(seqA1,predict(fitA1.MCMC_Burn,newdata = data.frame(A1=seqA1)),lty=5,col=6,lwd=2)
lines(seqA1,predict(fitA1.MCMC_Emp,newdata = data.frame(A1=seqA1)),lty=6,col=1)
legend("topright",c(expression("CLT"["Ann"]),
                    expression("CLT"["Burn"]),
                    expression("CLT"["Emp"]),
                    expression("MCMC"["Ann"]),
                    expression("MCMC"["Burn"]),
                    expression("MCMC"["Emp"])),inset = 0.05,lty=1:6,col=c(1:4,6,1),lwd=c(1,1,1,2,2,1))
abline(h=1)
par(mar=c(4,6,7,m.b))
plot(stats1D$A1,stats1D$A1_ERR_MCMC_Emp/stats1D$A1_RMSE,main=expression("Loess fit to Hessian"~hat(sigma)/sigma["RMSE"]),
     ylab="Ratio of Error Approximation to RMSE",xlab="Generated Slope (1D, 2PNO)",ylim=YR,type="n")
lines(seqA1,predict(fitA1.LMI_Iter,newdata = data.frame(A1=seqA1)),lty=1,col=1)
lines(seqA1,predict(fitA1.LMI_Iter_Emp_2PL,newdata = data.frame(A1=seqA1)),lty=2,col=2)
lines(seqA1,predict(fitA1.LMI_Iter_Emp_2PNO,newdata = data.frame(A1=seqA1)),lty=3,col=3,lwd=2)
lines(seqA1,predict(fitA1.Simple_2PL,newdata = data.frame(A1=seqA1)),lty=4,col=4)
lines(seqA1,predict(fitA1.Simple_2PNO,newdata = data.frame(A1=seqA1)),lty=5,col=6)
legend("topright",c(expression("ICE"),
                    expression("IPCE"["2PL"]),
                    expression("IPCE"["2PNO"]),
                    expression("SPCE"["2PL"]),
                    expression("SPCE"["2PNO"])),inset = 0.05,lty=1:5,col=c(1:4,6),lwd=c(1,1,2,1,1))
abline(h=1)

#########  Error Plots for B
p=2
YR<-range(c(stats1D$B_ERR_CLT_Ann,
            stats1D$B_ERR_CLT_Burn,
            stats1D$B_ERR_CLT_Emp,
            stats1D$B_ERR_LMI_Iter,
            stats1D$B_ERR_LMI_Iter_Emp_2PL,
            stats1D$B_ERR_LMI_Iter_Emp_2PNO,
            stats1D$B_ERR_Simple_2PL,
            stats1D$B_ERR_Simple_2PNO,
            stats1D$B_ERR_MCMC_Ann,
            stats1D$B_ERR_MCMC_Burn,
            stats1D$B_ERR_MCMC_Emp)/stats1D$B_RMSE)

ticks <- pretty(YR)[-c(length(pretty(YR))-1:0)]
labels <- format(ticks, big.mark=",", scientific=FALSE)
m.b = 0.25
par(mar=c(m.b,4,2,m.b))
# for (p in 1:ncol(gen.xi)) {
p=2
ifelse(p==ncol(gen.xi),pt<-"B",pt<-paste("A",p,sep=""))
plot(gen.xi[,p],rep(c(YR),sim.list[[d]]$J/2),type="n",
     main=expression("Intercept Error Estimates : 2PNO"),ylim=YR,
     xlab=NA,ylab=NA,xaxt="n",yaxt="n")
axis(2,at=ticks,labels=labels,las=2)
abline(h=1,lwd=2,lty=2)
points(gen.xi[,p],stats1D$B_ERR_CLT_Ann/stats1D$B_RMSE,pch=0,col=1)
text(-1,0.9*YR[2],txt.labels[1])

par(mar=c(m.b,4,m.b,m.b))
plot(gen.xi[,p],rep(c(YR),sim.list[[d]]$J/2),type="n",
     main=NA,ylim=YR,xlab=NA,ylab=NA,xaxt="n",yaxt="n")
axis(2,at=ticks,labels=labels,las=2)
abline(h=1,lwd=2,lty=2)
points(gen.xi[,p],stats1D$B_ERR_CLT_Burn/stats1D$B_RMSE,pch=1,col=1)
text(-1,0.9*YR[2],txt.labels[2])

plot(gen.xi[,p],rep(c(YR),sim.list[[d]]$J/2),type="n",
     main=NA,ylim=YR,xlab=NA,ylab=NA,xaxt="n",yaxt="n")
axis(2,at=ticks,labels=labels,las=2)
abline(h=1,lwd=2,lty=2)
points(gen.xi[,p],stats1D$B_ERR_CLT_Emp/stats1D$B_RMSE,pch=2,col=1)
text(-1,0.9*YR[2],txt.labels[3])

plot(gen.xi[,p],rep(c(YR),sim.list[[d]]$J/2),type="n",
     main=NA,ylim=YR,xlab=NA,ylab=NA,xaxt="n",yaxt="n")
axis(2,at=ticks,labels=labels,las=2)
abline(h=1,lwd=2,lty=2)
points(gen.xi[,p],stats1D$B_ERR_MCMC_Ann/stats1D$B_RMSE,pch=8,col=1)
text(-1,0.9*YR[2],txt.labels[4])

plot(gen.xi[,p],rep(c(YR),sim.list[[d]]$J/2),type="n",
     main=NA,ylim=YR,xlab=NA,ylab=NA,xaxt="n",yaxt="n")
axis(2,at=ticks,labels=labels,las=2)
abline(h=1,lwd=2,lty=2)
points(gen.xi[,p],stats1D$B_ERR_MCMC_Burn/stats1D$B_RMSE,pch=9,col=1)
text(-1,0.9*YR[2],txt.labels[5])

plot(gen.xi[,p],rep(c(YR),sim.list[[d]]$J/2),type="n",
     main=NA,ylim=YR,xlab=NA,ylab=NA,yaxt="n",xaxt="n")
axis(2,at=ticks,labels=labels,las=2)
abline(h=1,lwd=2,lty=2)
points(gen.xi[,p],stats1D$B_ERR_MCMC_Emp/stats1D$B_RMSE,pch=10,col=1)
text(-1,0.9*YR[2],txt.labels[6])
mtext(text = "Ratio of Error Approximation to RMSE",side = 2,line = 2.5,cex=0.7,adj=0)

plot(gen.xi[,p],rep(c(YR),sim.list[[d]]$J/2),type="n",
     main=NA,ylim=YR,xlab=NA,ylab=NA,xaxt="n",yaxt="n")
axis(2,at=ticks,labels=labels,las=2)
abline(h=1,lwd=2,lty=2)
points(gen.xi[,p],stats1D$B_ERR_LMI_Iter/stats1D$B_RMSE,pch=3,col=1)
text(-1,0.9*YR[2],txt.labels[7])

plot(gen.xi[,p],rep(c(YR),sim.list[[d]]$J/2),type="n",
     main=NA,ylim=YR,xlab=NA,ylab=NA,xaxt="n",yaxt="n")
axis(2,at=ticks,labels=labels,las=2)
abline(h=1,lwd=2,lty=2)
points(gen.xi[,p],stats1D$B_ERR_LMI_Iter_Emp_2PL/stats1D$B_RMSE,pch=4,col=1)
text(-1,0.9*YR[2],txt.labels[8])

plot(gen.xi[,p],rep(c(YR),sim.list[[d]]$J/2),type="n",
     main=NA,ylim=YR,xlab=NA,ylab=NA,xaxt="n",yaxt="n")
axis(2,at=ticks,labels=labels,las=2)
abline(h=1,lwd=2,lty=2)
points(gen.xi[,p],stats1D$B_ERR_LMI_Iter_Emp_2PNO/stats1D$B_RMSE,pch=5,col=1)
text(-1,0.9*YR[2],txt.labels[9])

plot(gen.xi[,p],rep(c(YR),sim.list[[d]]$J/2),type="n",
     main=NA,ylim=YR,xlab=NA,ylab=NA,xaxt="n",yaxt="n")
axis(2,at=ticks,labels=labels,las=2)
abline(h=1,lwd=2,lty=2)
points(gen.xi[,p],stats1D$B_ERR_Simple_2PL/stats1D$B_RMSE,pch=6,col=1)
text(-1,0.9*YR[2],txt.labels[10])

par(mar=c(4,4,m.b,m.b))
plot(gen.xi[,p],rep(c(YR),sim.list[[d]]$J/2),type="n",
     main=NA,ylim=YR,xlab=NA,ylab=NA,yaxt="n")
axis(2,at=ticks,labels=labels,las=2)
abline(h=1,lwd=2,lty=2)
points(gen.xi[,p],stats1D$B_ERR_Simple_2PNO/stats1D$B_RMSE,pch=7,col=1)
text(-1,0.9*YR[2],txt.labels[11])


fitB.CLT_Ann<-loess(B_ERR_CLT_Ann/B_RMSE~B,data=stats1D,span=0.5)
fitB.CLT_Ann
sd(fitB.CLT_Ann$residuals)
#fitBo<-lm(B_sd_LMIo/B_sd~B,data=stats1D)
fitB.CLT_Burn<-loess(B_ERR_CLT_Burn/B_RMSE~B,data=stats1D,span=0.5)
fitB.CLT_Burn
sd(fitB.CLT_Burn$residuals)
#fitBMC<-lm(B_sd_MCE/B_sd~B,data=stats1D)
fitB.CLT_Emp<-loess(B_ERR_CLT_Emp/B_RMSE~B,data=stats1D,span=0.5)
fitB.CLT_Emp
sd(fitB.CLT_Emp$residuals)
#fitBLM<-lm(B_sd_LMI/B_sd~B,data=stats1D)
fitB.LMI_Iter<-loess(B_ERR_LMI_Iter/B_RMSE~B,data=stats1D,span=0.5)
fitB.LMI_Iter
sd(fitB.LMI_Iter$residuals)
fitB.LMI_Iter_Emp_2PL<-loess(B_ERR_LMI_Iter_Emp_2PL/B_RMSE~B,data=stats1D,span=0.5)
fitB.LMI_Iter_Emp_2PL
sd(fitB.LMI_Iter_Emp_2PL$residuals)
#fitBo<-lm(B_sd_LMIo/B_sd~B,data=stats1D)
fitB.LMI_Iter_Emp_2PNO<-loess(B_ERR_LMI_Iter_Emp_2PNO/B_RMSE~B,data=stats1D,span=0.5)
fitB.LMI_Iter_Emp_2PNO
sd(fitB.LMI_Iter_Emp_2PNO$residuals)
#fitBMC<-lm(B_sd_MCE/B_sd~B,data=stats1D)
fitB.Simple_2PL<-loess(B_ERR_Simple_2PL/B_RMSE~B,data=stats1D,span=0.5)
fitB.Simple_2PL
sd(fitB.Simple_2PL$residuals)
#fitBLM<-lm(B_sd_LMI/B_sd~B,data=stats1D)
fitB.Simple_2PNO<-loess(B_ERR_Simple_2PNO/B_RMSE~B,data=stats1D,span=0.5)
fitB.Simple_2PNO
sd(fitB.Simple_2PNO$residuals)
fitB.MCMC_Ann<-loess(B_ERR_MCMC_Ann/B_RMSE~B,data=stats1D,span=0.5)
fitB.MCMC_Ann
sd(fitB.MCMC_Ann$residuals)
#fitBo<-lm(B_sd_LMIo/B_sd~B,data=stats1D)
fitB.MCMC_Burn<-loess(B_ERR_MCMC_Burn/B_RMSE~B,data=stats1D,span=0.5)
fitB.MCMC_Burn
sd(fitB.MCMC_Burn$residuals)
#fitBMC<-lm(B_sd_MCE/B_sd~B,data=stats1D)
fitB.MCMC_Emp<-loess(B_ERR_MCMC_Emp/B_RMSE~B,data=stats1D,span=0.5)
fitB.MCMC_Emp
sd(fitB.MCMC_Emp$residuals)
######plot2
par(mar=c(4,6,4,m.b))
seqB<-seq(range(stats1D$B)[1],range(stats1D$B)[2],length.out = 50)
plot(stats1D$B,stats1D$B_ERR_MCMC_Emp/stats1D$B_RMSE,main=expression("Loess fit to"~hat(sigma)/sigma["RMSE"]),
     ylab="Ratio of Error Approximation to RMSE",xlab="Generated Slope (1D, 2PNO)",ylim=YR,type="n")
lines(seqB,predict(fitB.CLT_Ann,newdata = data.frame(B=seqB)),lty=1)
lines(seqB,predict(fitB.CLT_Burn,newdata = data.frame(B=seqB)),lty=2,col=2)
lines(seqB,predict(fitB.CLT_Emp,newdata = data.frame(B=seqB)),lty=3,col=3)
lines(seqB,predict(fitB.MCMC_Ann,newdata = data.frame(B=seqB)),lty=4,col=4,lwd=2)
lines(seqB,predict(fitB.MCMC_Burn,newdata = data.frame(B=seqB)),lty=5,col=6,lwd=2)
lines(seqB,predict(fitB.MCMC_Emp,newdata = data.frame(B=seqB)),lty=6,col=1)
legend("topright",c(expression("CLT"["Ann"]),
                    expression("CLT"["Burn"]),
                    expression("CLT"["Emp"]),
                    expression("MCMC"["Ann"]),
                    expression("MCMC"["Burn"]),
                    expression("MCMC"["Emp"])),inset = 0.05,lty=1:6,col=c(1:4,6,1),lwd=c(1,1,1,2,2,1))
abline(h=1)
par(mar=c(4,6,7,m.b))
plot(stats1D$B,stats1D$B_ERR_MCMC_Emp/stats1D$B_RMSE,main=expression("Loess fit to Hessian"~hat(sigma)/sigma["RMSE"]),
     ylab="Ratio of Error Approximation to RMSE",xlab="Generated Slope (1D, 2PNO)",ylim=YR,type="n")
lines(seqB,predict(fitB.LMI_Iter,newdata = data.frame(B=seqB)),lty=1,col=1)
lines(seqB,predict(fitB.LMI_Iter_Emp_2PL,newdata = data.frame(B=seqB)),lty=2,col=2)
lines(seqB,predict(fitB.LMI_Iter_Emp_2PNO,newdata = data.frame(B=seqB)),lty=3,col=3,lwd=2)
lines(seqB,predict(fitB.Simple_2PL,newdata = data.frame(B=seqB)),lty=4,col=4)
lines(seqB,predict(fitB.Simple_2PNO,newdata = data.frame(B=seqB)),lty=5,col=6)
legend("topright",c(expression("ICE"),
                    expression("IPCE"["2PL"]),
                    expression("IPCE"["2PNO"]),
                    expression("SPCE"["2PL"]),
                    expression("SPCE"["2PNO"])),inset = 0.05,lty=1:5,col=c(1:4,6),lwd=c(1,1,2,1,1))
abline(h=1)

###########################
# DEFINITIONS
#
# iError        = IRT  approx Louis Missing Info
# oError        = 2PNO approx Louis Missing Info
# EmpSE$VARLMIi = IRT  approx Louis Missing Info
# EmpSE$VARLMIo = 2PNO approx Louis Missing Info
# MCERR         = MCMC variance 

MCErrFunc1Run<-function(x,burnin=burn,start=400) {
  for (iii in 1:dim(x$Aiter)[2]) {
    Atrim = start+10*1:(floor((burnin-start)/10))
    if (iii==1) {
      AE<-apply(x$Aiter[,1,Atrim],1,sd)
    } else {
      AE<-cbind(AE,apply(x$Aiter[,1,Atrim],1,sd))
    }
  }
  Btrim = start+7*1:(floor((burnin-start)/7))
  BE<-apply(x$Biter[,Btrim],1,sd)
  1.96*cbind(AE,BE)
}

###  
# Wanna test statistics about the negative variances?
###
# ng = c(0,0)
# for (iii in c(20,50,80)) {
#   for (jjj in c(2000,5000,8000)) {
#     for (kkk in 1:25) {
#       estfile=paste("Error500/FitError1D_",iii,"_",jjj,"_s1_r",kkk,".rda",sep="")
#       estfile=paste("Error500/FitError1D_",iii,"_",jjj,"_s1_r",kkk,".rda",sep="")
#       load(estfile)
#       ng = c(colSums(FitDATA$xiError<0)+ng)
#     }
#   }
# }
# print(ng/(25*20))
par(mfrow=c(1,1))
R=c(); C=c()
ITEMS = 20; STUDS = 2000; SIM = 1
for (jj in 1:50) {
estfile=paste("Error500/FitError1D-",ITEMS,"-",STUDS,"-s",SIM,"-r",jj,".rda",sep="")
load(estfile)
attach(FitDATA)
YY = cbind(FitDATA$XI,FitDATA$xi,FitDATA$EmpSE$VARLMIo)
YY = YY[order(YY[,2]),]
if (jj==1) {
  plot(YY[,2],YY[,5],ylim = c(-0.05,0.05))
} else {
    points(YY[,2],YY[,5])
    lines(YY[,2],YY[,5],col=rainbow(20,start = 0,end=0.7)[jj])
}
R = c(R,cor(abs(YY[,2]),YY[,6]))
C = c(C,sum(log(1+YY[,6])))
detach(FitDATA)
}
plot(C,R)


for (i in 1:50) {
estfile=paste("Error500/FitError1D-20-2000-s2-r",i,".rda",sep="")
load(estfile)
attach(FitDATA)
if(sum(EmpSE$VARLMIo[,1]>0)==20) print(i)
detach(FitDATA)
}


Plot1DRunErrors = function(J,N,S,R) {
  estfile=paste("Error500/FitError1D-",J,"-",N,"-s",S,"-r",R,".rda",sep="")
  load(estfile)
  attach(FitDATA)
  # sqrt(EmpSE$VARLMIo[,1])
  # detach(FitDATA)
  RMSE = stats1D[(stats1D$J==J&stats1D$N==N),"A1_sd"]
  MCERR = MCErrFunc1Run(MCMCDATA)
  pdf(paste0("Stats1DFit_",J,"_",N,".pdf"),width=10,height=10)
  par(mfrow=c(2,2),mar=c(4,5,1,1))
  plot(XI[,1],XI[,1]-xi[,1],pch=19,cex=1,xlim=range(XI[,1]),ylim=range(c((XI[,1]-xi[,1])-2*RMSE,(XI[,1]-xi[,1])+2*RMSE)),
     main=NA,yaxt="n",xaxt="n",
     xlab="Generated Slope",ylab=expression(hat(A)-A %+-% 2*sigma[RMSE]))
  axis(2,cex.axis=0.8,las=2)
  axis(1,cex.axis=0.8)
  arrows(XI[,1],(XI[,1]-xi[,1])-2*RMSE,XI[,1],(XI[,1]-xi[,1])+2*RMSE,code=3,angle=90,length=0.07)
  abline(h=0,lty=2)

  par(mar=c(4,2,1,1))
  plot(density(sqrt(EmpSE$VARLMIo[,1])/RMSE),yaxt="n",xaxt="n",ylab=NA,lty=3,
     xlim=c(0.25,2),main=NA,xlab=expression(hat(sigma)[A]/sigma[RMSE]))
  mtext("Density Approximation",2,line=1,cex=0.7)
  axis(1,cex.axis=0.8)
  lines(density(sqrt(oError[,1])/RMSE),lty=5)
  lines(density(sqrt(EmpSE$VARLMIi[,1])/RMSE),lty=6)
  lines(density(sqrt(iError[,1])/RMSE),lty=4)
  lines(density(sqrt(xiError[,1])/RMSE),lty=2)
  lines(density(sqrt(EmpSE$VARLMIo[,1])/RMSE),lty=3)
  lines(density(MCERR[,1]/RMSE),lty=1)
  legend("topright",legend = c("MCMC","SAEM","2PLC","2PNC","2PLI","2PNI"),
       lty=c(1,2,4,5,6,3),title="Error Estimation")

  RMSE = stats1D[(stats1D$J==J&stats1D$N==N),"B_sd"]
  par(mar=c(4,5,1,1))
  plot(XI[,2],XI[,2]-xi[,2],pch=19,cex=1,xlim=range(XI[,2]),
     ylim=range(c((XI[,2]-xi[,2])-2*RMSE,(XI[,2]-xi[,2])+2*RMSE)),
     main=NA,yaxt="n",xaxt="n",
     xlab="Generated Intercept",ylab=expression(hat(b)-b %+-% 2*sigma[RMSE]))
  axis(2,cex.axis=0.8,las=2)
  axis(1,cex.axis=0.8)
  arrows(XI[,2],(XI[,2]-xi[,2])-2*RMSE,XI[,2],(XI[,2]-xi[,2])+2*RMSE,code=3,angle=90,length=0.07)
  abline(h=0,lty=2)

  par(mar=c(4,2,1,1))
  plot(density(sqrt(iError[,2])/RMSE),yaxt="n",xaxt="n",ylab=NA,lty=4,
     xlim=c(0.2,1.9),main=NA,
     xlab=expression(hat(sigma)[b]/sigma[RMSE]))
  axis(1,cex.axis=0.8)
  lines(density(MCERR[,2]/RMSE),lty=1)
  mtext("Density Approximation",2,line=1,cex=0.7)
  lines(density(sqrt(EmpSE$VARLMIo[,2])/RMSE),lty=3)
  lines(density(sqrt(xiError[,2])/RMSE),lty=2)
  lines(density(sqrt(oError[,2])/RMSE),lty=5)
  lines(density(sqrt(EmpSE$VARLMIi[,2])/RMSE),lty=6)
  legend("topright",legend = c("MCMC","SAEM","2PLC","2PNC","2PLI","2PNI"),
       lty=c(1,2,4,5,6,3),title="Error Estimation")
  dev.off()

  detach(FitDATA)
}

Plot1DRunErrors(J=50,N=5000,S=2,R=5)
# 
# 
# estfile=paste("Error500/FitError1D-20-2000-s1-r5.rda",sep="")
# load(estfile)
# 
# 
# 
# 
# nf <- layout(matrix(c(1,2,3,4,5,6),3,2), c(7,7), c(1.5,3,4), TRUE)
# layout.show(nf)
# par(cex=0.8)
# par(mar=c(m.b,5,2,m.b))
# plot(FitDATA$XI[,1],FitDATA$A,main="A",xaxt="n",ylab="Estimated Slope",las=2)  #,xlab="Generated Slope"
# abline(a=0,b=1)
# par(mar=c(m.b,5,m.b,m.b))
# plot(FitDATA$XI[,1],sqrt(FitDATA$xiError[,1]),xaxt="n",ylim=c(0.0,0.2),las=2,ylab="Calculated Standard Error")#,xlab="Generated Slope")
# points(FitDATA$XI[,1],sqrt(FitDATA$iError[,1]),pch=2)
# points(FitDATA$XI[,1],sqrt(FitDATA$oError[,1]),pch=3)
# points(FitDATA$XI[,1],FitDATA$EmpSE$SEA,pch=4)
# points(FitDATA$XI[,1],FitDATA$EmpSE$MCSA,pch=5)
# points(FitDATA$XI[,1],sqrt(FitDATA$EmpSE$VARLMIi[,1]),pch=6)
# points(FitDATA$XI[,1],sqrt(FitDATA$EmpSE$VARLMIo[,1]),pch=22)
# points(stats1D[(stats1D$J==20&stats1D$N==2000),c("A1","A1_sd")],pch=18,cex=1.5)
# legend("topleft",c("Iterative (RM)","2PL (At Convergence)","2PNO (At Convergence)","MCMC Chain (Post Convergence)","MCMC CLT","2PL (Post Convergence)","2PNO (Post Convergence)","RMSE"),pch=c(1,2,3,4,5,6,22,18),cex=0.85)
# par(mar=c(5,5,m.b,m.b))
# plot(FitDATA$XI[,1],sqrt(FitDATA$xiError[,1])/stats1D[(stats1D$J==20&stats1D$N==2000),c("A1_sd")],las=2,xlab="Generated Slope",ylab=expression(paste("Ratio ",hat(sigma)/sigma[RMSE])),ylim=c(0,2.5))
# points(FitDATA$XI[,1],sqrt(FitDATA$iError[,1])/stats1D[(stats1D$J==20&stats1D$N==2000),c("A1_sd")],pch=2)
# points(FitDATA$XI[,1],sqrt(FitDATA$oError[,1])/stats1D[(stats1D$J==20&stats1D$N==2000),c("A1_sd")],pch=3)
# points(FitDATA$XI[,1],FitDATA$EmpSE$SEA/stats1D[(stats1D$J==20&stats1D$N==2000),c("A1_sd")],pch=4)
# points(FitDATA$XI[,1],FitDATA$EmpSE$MCSA/stats1D[(stats1D$J==20&stats1D$N==2000),c("A1_sd")],pch=5)
# points(FitDATA$XI[,1],sqrt(FitDATA$EmpSE$VARLMIi[,1])/stats1D[(stats1D$J==20&stats1D$N==2000),c("A1_sd")],pch=6)
# points(FitDATA$XI[,1],sqrt(FitDATA$EmpSE$VARLMIo[,1])/stats1D[(stats1D$J==20&stats1D$N==2000),c("A1_sd")],pch=22)
# abline(h=1,lty=2)
# par(mar=c(m.b,m.b,2,5))
# plot(FitDATA$XI[,2],FitDATA$B,main="b",yaxt="n",xaxt="n")#,xlab="Generated Intercept"
# abline(a=0,b=1)
# axis(4,las=2)
# mtext(text = "Estimated Intercept",side = 4,line = 3,cex=0.8)
# par(mar=c(m.b,m.b,m.b,5))
# plot(FitDATA$XI[,2],sqrt(FitDATA$xiError[,2]),yaxt="n",xaxt="n",ylim=c(0.0,0.25),xlim=c(-2,2))
# axis(4,las=2)
# mtext(text ="Calculated Standard Error",side = 4,line = 3,cex=0.8)
# points(FitDATA$XI[,2],sqrt(FitDATA$iError[,2]),pch=2)
# points(FitDATA$XI[,2],sqrt(FitDATA$oError[,2]),pch=3)
# points(FitDATA$XI[,2],FitDATA$EmpSE$SEB,pch=4)
# points(FitDATA$XI[,2],FitDATA$EmpSE$MCSB,pch=5)
# points(FitDATA$XI[,2],sqrt(FitDATA$EmpSE$VARLMIi[,2]),pch=6)
# points(FitDATA$XI[,2],sqrt(FitDATA$EmpSE$VARLMIo[,2]),pch=22)
# points(stats1D[(stats1D$J==20&stats1D$N==2000),c("B","B_sd")],pch=18,cex=1.5)
# 
# par(mar=c(5,m.b,m.b,5))
# plot(FitDATA$XI[,2],sqrt(FitDATA$xiError[,2])/stats1D[(stats1D$J==20&stats1D$N==2000),c("B_sd")],xlab="Generated Intercept",ylim=c(0,2.5),yaxt="n")
# axis(4,las=2)
# mtext(text = expression(paste("Ratio ",hat(sigma)/sigma[RMSE])),side = 4,line = 3,cex=0.8)
# points(FitDATA$XI[,2],sqrt(FitDATA$iError[,2])/stats1D[(stats1D$J==20&stats1D$N==2000),c("B_sd")],pch=2)
# points(FitDATA$XI[,2],sqrt(FitDATA$oError[,2])/stats1D[(stats1D$J==20&stats1D$N==2000),c("B_sd")],pch=3)
# points(FitDATA$XI[,2],FitDATA$EmpSE$SEB/stats1D[(stats1D$J==20&stats1D$N==2000),c("B_sd")],pch=4)
# points(FitDATA$XI[,2],FitDATA$EmpSE$MCSB/stats1D[(stats1D$J==20&stats1D$N==2000),c("B_sd")],pch=5)
# points(FitDATA$XI[,2],sqrt(FitDATA$EmpSE$VARLMIi[,2])/stats1D[(stats1D$J==20&stats1D$N==2000),c("B_sd")],pch=6)
# points(FitDATA$XI[,2],sqrt(FitDATA$EmpSE$VARLMIo[,2])/stats1D[(stats1D$J==20&stats1D$N==2000),c("B_sd")],pch=22)
# abline(h=1,lty=2)
# #points(,pch=16)
# estfile=paste("Error500/FitError1D_20_2000_s1_r13.rda",sep="")
# load(estfile)
# RMSE = stats1D[(stats1D$J==20&stats1D$N==2000),"A1_sd"]
# MCERR = MCErrFunc1Run(MCMCDATA)
# 
# A20<-cbind(FitDATA$XI[,1],FitDATA$A,RMSE,
#       sqrt(FitDATA$xiError[,1])/RMSE,
#       sqrt(FitDATA$iError[,1])/RMSE,
#       sqrt(FitDATA$oError[,1])/RMSE,
#       MCERR[,1]/RMSE,
#       sqrt(FitDATA$EmpSE$VARLMIi[,1])/RMSE,
#       sqrt(FitDATA$EmpSE$VARLMIo[,1])/RMSE)
# colnames(A20)<-c("A_Gen","A_Fit","RMSE","IterativeRM","TwoPLAtConv","TwoPNOAtConv","MCMCChain","TwoPLPostConv","TwoPNOPostConv")
# write.csv(A20,"A20_Errors.csv")
# RMSE = stats1D[(stats1D$J==20&stats1D$N==2000),"B_sd"]
# 
# B20<-cbind(FitDATA$XI[,2],FitDATA$B,RMSE,
#       sqrt(FitDATA$xiError[,2])/RMSE,
#       sqrt(FitDATA$iError[,2])/RMSE,
#       sqrt(FitDATA$oError[,2])/RMSE,
#       MCERR[,2]/RMSE,
#       sqrt(FitDATA$EmpSE$VARLMIi[,2])/RMSE,
#       sqrt(FitDATA$EmpSE$VARLMIo[,2])/RMSE)
# colnames(B20)<-c("B_Gen","B_Fit","RMSE","IterativeRM","TwoPLAtConv","TwoPNOAtConv","MCMCChain","TwoPLPostConv","TwoPNOPostConv")
# write.csv(B20,"B20_Errors.csv")
# 
#nrow(stats1D)
# stats[361:720,]
# rbind(c(50,2000,colMeans(stats[stats$N==2000&stats$Burn==50,c(5,6,8,9)])),
#       c(50,3000,colMeans(stats[stats$N==3000&stats$Burn==50,c(5,6,8,9)])),
#       c(50,4000,colMeans(stats[stats$N==4000&stats$Burn==50,c(5,6,8,9)])),
#       c(100,2000,colMeans(stats[stats$N==2000&stats$Burn==100,c(5,6,8,9)])),
#       c(100,3000,colMeans(stats[stats$N==3000&stats$Burn==100,c(5,6,8,9)])),
#       c(100,4000,colMeans(stats[stats$N==4000&stats$Burn==100,c(5,6,8,9)])))

#detach(FitDATA)
###############
# RUN 2D Test #
###############
sims=50
# items=c(30,40,50)
# examinees=c(2000,4000) # still have to run 50-3000 ...3000,
# burn = 100
#b = 1; n = 1; j = 1; i = 1 # debug

items=c(20,50,80)
examinees=c(2000,5000,8000) #,10000)
burn = 500
SIM = 1

#for (b in 1:length(burn)) {
for (j in 1:length(items)) {
  for (n in 1:length(examinees)) {
    sims=50
    simfile=paste("Error500/SimError2D-",items,"-",examinees,"-s",SIM,".rda",sep="")    # Save to an ".rda" file with this name
    load(simfile)
    estfile=paste("Error500/FitError2D-",items,"-",examinees,"-s",SIM,"-r1.rda",sep="") 
    load(estfile)
    XIiter<-array(FitDATA$xi, dim=c(items,1+FitDATA$settings$Adim,sims))        
    ERRiter<-array(FitDATA$xiError, dim=c(items,1+FitDATA$settings$Adim,sims))        
    iERR<-array(FitDATA$iError, dim=c(items,1+FitDATA$settings$Adim,sims))        
    oERR<-array(FitDATA$oError, dim=c(items,1+FitDATA$settings$Adim,sims))        
    MCERR<-array(MCErrorFunction(MCMCDATA), dim=c(items,1+FitDATA$settings$Adim,sims))
      #array(cbind(FitDATA$EmpSE$MCSA,FitDATA$EmpSE$MCSB), dim=c(items,1+FitDATA$settings$Adim,sims))
    EmpLMIi<-array(FitDATA$EmpSE$VARLMIi, dim=c(items,1+FitDATA$settings$Adim,sims))
    EmpLMIo<-array(FitDATA$EmpSE$VARLMIo, dim=c(items,1+FitDATA$settings$Adim,sims))
    k=0
    for (i in 1:sims) {
      estfile=paste("Error500/FitError2D-",items,"-",examinees,"-s",SIM,"-r",i,".rda",sep="") 
      if (file.exists(estfile)) {
        k=k+1
        load(estfile)
        XIiter[,,k]<-as.matrix(FitDATA$xi) 
        ERRiter[,,k]<-as.matrix(FitDATA$xiError) 
        iERR[,,k]<-as.matrix(FitDATA$iError) 
        oERR[,,k]<-as.matrix(FitDATA$oError) 
        MCERR[,,k]<-as.matrix(MCErrorFunction(MCMCDATA))
          #as.matrix(cbind(FitDATA$EmpSE$MCSA,FitDATA$EmpSE$MCSB)) 
        EmpLMIi[,,k]<-as.matrix(FitDATA$EmpSE$VARLMIi) 
        EmpLMIo[,,k]<-as.matrix(FitDATA$EmpSE$VARLMIo) 
      }
      # if (i==25 & k<25) {
      #   XIiter<-XIiter[,,1:k]
      #   ERRiter<-ERRiter[,,1:k]
      #   iERR<-iERR[,,1:k]
      #   oERR<-oERR[,,1:k]
      #   MCERR<-MCERR[,,1:k]
      #   EmpLMIi<-EmpLMIi[,,1:k]
      #   EmpLMIo<-EmpLMIo[,,1:k]
      # }
    }
    sims<-k
    mXI<-apply(XIiter,c(1,2),mean)
    mMCE<-apply(MCERR,c(1,2),mean)
    miERR<-apply(iERR,c(1,2),mean)
    moERR<-apply(oERR,c(1,2),mean)
    mLMI<-apply(ERRiter,c(1,2),mean)
    mLMIi<-apply(EmpLMIi,c(1,2),mean)
    mLMIo<-apply(EmpLMIo,c(1,2),mean)
    
    eiERR<-matrix(sapply(miERR,negsqrt),nrow=items,ncol=1+FitDATA$settings$Adim)
    eoERR<-matrix(sapply(moERR,negsqrt),nrow=items,ncol=1+FitDATA$settings$Adim)
    eLMI<-matrix(sapply(mLMI,negsqrt),nrow=items,ncol=1+FitDATA$settings$Adim)
    eLMIi<-matrix(sapply(mLMIi,negsqrt),nrow=items,ncol=1+FitDATA$settings$Adim)
    eLMIo<-matrix(sapply(mLMIo,negsqrt),nrow=items,ncol=1+FitDATA$settings$Adim)
    
    bias<-mXI-gen.xi
    print(dim(XIiter))
    sXI<-apply(XIiter-array(rep(mXI,sims),dim=c(items,1+FitDATA$settings$Adim,sims)),c(1,2),sd)
    sMCE<-apply(MCERR-array(rep(mMCE,sims),dim=c(items,1+FitDATA$settings$Adim,sims)),c(1,2),sd)
    siERR<-apply(iERR-array(rep(miERR,sims),dim=c(items,1+FitDATA$settings$Adim,sims)),c(1,2),sd)
    soERR<-apply(oERR-array(rep(moERR,sims),dim=c(items,1+FitDATA$settings$Adim,sims)),c(1,2),sd)
    sLMI<-apply(ERRiter-array(rep(mLMI,sims),dim=c(items,1+FitDATA$settings$Adim,sims)),c(1,2),sd)
    sLMIi<-apply(EmpLMIi-array(rep(mLMIi,sims),dim=c(items,1+FitDATA$settings$Adim,sims)),c(1,2),sd)
    sLMIo<-apply(EmpLMIo-array(rep(mLMIo,sims),dim=c(items,1+FitDATA$settings$Adim,sims)),c(1,2),sd)
    par(mfrow=c(ncol(gen.xi),2))
    for (p in 1:ncol(gen.xi)) {
      ifelse(p==ncol(gen.xi),pt<-"B",pt<-paste("A",p,sep=""))
      plot(gen.xi[,p],bias[,p],pch=19,cex=1.5,xlim=range(gen.xi[,p]),ylim=range(c(bias[,p]-2*sXI[,p],bias[,p]+2*sXI[,p])),
           main=paste("Parameter",pt,"Bias & Error (",sims," Runs)\n",items,"items; ",examinees,"examinees; ",burn[b],"burn-in"),
           xlab=pt,ylab=paste(pt,"hat"))
      arrows(gen.xi[,p],bias[,p]-2*sXI[,p],gen.xi[,p],bias[,p]+2*sXI[,p],code=3,angle=90,length=0.07)
      plot(gen.xi[,p],eLMI[,p],type="n",main=paste("Parameter",pt,"Approx Errors"),ylim=range(c(0,eLMI[,p],mMCE[,p],eLMIi[,p],eLMIo[,p])),xlab=pt,ylab="Approximate Error")
      points(gen.xi[,p],eiERR[,p],pch=16,col=3)
      points(gen.xi[,p],eoERR[,p],pch=16,col=1)
      points(gen.xi[,p],eLMI[,p],pch=16,col=2)
      points(gen.xi[,p],mMCE[,p],pch=16,col=4)
      points(gen.xi[,p],eLMIi[,p],pch=16,col=6)
      points(gen.xi[,p],eLMIo[,p],pch=16,col=5)
      points(gen.xi[,p],sXI[,p])
      legend("top",c("ConvIRT","MCMC","PostIRT","Post2PNO","OneIRT","One2PNO","RMSE"),pch=c(16,16,16,16,16,16,1),col=c(2,4,6,5,3,1,1))
    }
    if (j==1&n==1) {
      stats2D<-data.frame(N=rep(examinees,items),J=rep(items,items),Burn=rep(burn[b],items))
      for(p in 1:ncol(gen.xi)) {
        ifelse(p==ncol(gen.xi),pt<-"B",pt<-paste("A",p,sep=""))
        stats2D<-cbind(stats2D,gen.xi[,p])
        colnames(stats2D)[ncol(stats2D)]<-pt
        stats2D<-cbind(stats2D,bias[,p])
        colnames(stats2D)[ncol(stats2D)]<-paste(pt,"bias",sep="_")
        stats2D<-cbind(stats2D,sXI[,p])
        colnames(stats2D)[ncol(stats2D)]<-paste(pt,"sd",sep="_")
        stats2D<-cbind(stats2D,mLMI[,p])
        colnames(stats2D)[ncol(stats2D)]<-paste(pt,"VAR_LMI",sep="_")
        stats2D<-cbind(stats2D,eLMI[,p])
        colnames(stats2D)[ncol(stats2D)]<-paste(pt,"sd_LMI",sep="_")
        stats2D<-cbind(stats2D,sLMI[,p])
        colnames(stats2D)[ncol(stats2D)]<-paste(pt,"sd_VAR_LMI",sep="_")
        stats2D<-cbind(stats2D,miERR[,p])
        colnames(stats2D)[ncol(stats2D)]<-paste(pt,"VAR_iERR",sep="_")
        stats2D<-cbind(stats2D,eiERR[,p])
        colnames(stats2D)[ncol(stats2D)]<-paste(pt,"sd_iERR",sep="_")
        stats2D<-cbind(stats2D,siERR[,p])
        colnames(stats2D)[ncol(stats2D)]<-paste(pt,"sd_VAR_iERR",sep="_")
        stats2D<-cbind(stats2D,moERR[,p])
        colnames(stats2D)[ncol(stats2D)]<-paste(pt,"VAR_oERR",sep="_")
        stats2D<-cbind(stats2D,eoERR[,p])
        colnames(stats2D)[ncol(stats2D)]<-paste(pt,"sd_oERR",sep="_")
        stats2D<-cbind(stats2D,soERR[,p])
        colnames(stats2D)[ncol(stats2D)]<-paste(pt,"sd_VAR_oERR",sep="_")
        stats2D<-cbind(stats2D,mMCE[,p])
        colnames(stats2D)[ncol(stats2D)]<-paste(pt,"sd_MCE",sep="_")
        stats2D<-cbind(stats2D,sMCE[,p])
        colnames(stats2D)[ncol(stats2D)]<-paste(pt,"sd_VAR_MCE",sep="_")
        stats2D<-cbind(stats2D,mLMIi[,p])
        colnames(stats2D)[ncol(stats2D)]<-paste(pt,"VAR_LMIi",sep="_")
        stats2D<-cbind(stats2D,eLMIi[,p])
        colnames(stats2D)[ncol(stats2D)]<-paste(pt,"sd_LMIi",sep="_")
        stats2D<-cbind(stats2D,sLMIi[,p])
        colnames(stats2D)[ncol(stats2D)]<-paste(pt,"sd_VAR_LMIi",sep="_")
        stats2D<-cbind(stats2D,mLMIo[,p])
        colnames(stats2D)[ncol(stats2D)]<-paste(pt,"VAR_LMIo",sep="_")
        stats2D<-cbind(stats2D,eLMIo[,p])
        colnames(stats2D)[ncol(stats2D)]<-paste(pt,"sd_LMIo",sep="_") 
        stats2D<-cbind(stats2D,sLMIo[,p])
        colnames(stats2D)[ncol(stats2D)]<-paste(pt,"sd_VAR_LMIo",sep="_")
      }
    } else {
      st<-cbind(rep(examinees,items),rep(items,items),rep(burn[b],items))
      for(p in 1:ncol(gen.xi)) {
        st<-cbind(st,gen.xi[,p],bias[,p],sXI[,p],mLMI[,p],sqrt(mLMI[,p]),sLMI[,p],
                  miERR[,p],sqrt(miERR[,p]),siERR[,p],moERR[,p],sqrt(moERR[,p]),
                  soERR[,p],mMCE[,p],sMCE[,p],
                  mLMIi[,p],eLMIi[,p],sLMIi[,p],mLMIo[,p],eLMIo[,p],sLMIo[,p])
      }
      colnames(st)<-colnames(stats2D)
      stats2D<-rbind(stats2D,st)
    }
  }
}

write.csv(stats2D,"Stats2D.csv")

for ()
FitDATA$EmpSE$VARLMIi
# estfile=paste("Error500/FitError2D-50-5000-s1-r2.rda",sep="")
# load(estfile)
# attach(FitDATA)
# cbind(XI[,1]-xi[,1],XI[,2]-xi[,2],XI[,3]-xi[,3])
# XI


Plot2DRunErrors = function(J,N,S,R,sMC = 350) {
  estfile=paste("Error500/FitError2D-",J,"-",N,"-s",S,"-r",R,".rda",sep="")
  load(estfile)
  attach(FitDATA)
  pdf(paste0("Stats2DFit_",J,"_",N,".pdf"),width=10,height=15)
  #detach(FitDATA)
  RMSE = stats2D[(stats2D$J==J&stats2D$N==N),"A1_sd"]
  MCERR = MCErrFunc1Run(MCMCDATA,start=sMC)
  par(mfrow=c(3,2),mar=c(4,5,1,1))
  plot(XI[,1],XI[,1]-xi[,1],pch=19,cex=1,xlim=range(XI[,1]),ylim=range(c((XI[,1]-xi[,1])-2*RMSE,(XI[,1]-xi[,1])+2*RMSE)),
     main=NA,yaxt="n",xaxt="n",
     xlab="Generated Slope A1",ylab=expression(hat(A1)-A1 %+-% 2*sigma[RMSE]))
  axis(2,cex.axis=0.8,las=2)
  axis(1,cex.axis=0.8)
  arrows(XI[,1],(XI[,1]-xi[,1])-2*RMSE,XI[,1],(XI[,1]-xi[,1])+2*RMSE,code=3,angle=90,length=0.07)
  abline(h=0,lty=2)

  par(mar=c(4,2,1,1))
  plot(density(sqrt(oError[,1])/RMSE),yaxt="n",xaxt="n",ylab=NA,lty=5,
       xlim=c(0.25,2),main=NA,xlab=expression(hat(sigma)[A1]/sigma[RMSE]))
  mtext("Density Approximation",2,line=1,cex=0.7)
  axis(1,cex.axis=0.8)
#  lines(density(sqrt(EmpSE$VARLMIo[,1])/RMSE),lty=3)
#  lines(density(sqrt(EmpSE$VARLMIi[,1])/RMSE),lty=6)
  lines(density(sqrt(iError[,1])/RMSE),lty=4)
  lines(density(sqrt(xiError[,1])/RMSE),lty=2)
  lines(density(sqrt(EmpSE$VARLMIo[,1])/RMSE),lty=3)
  lines(density(MCERR[,1]/RMSE),lty=1)
  legend("topright",legend = c("MCMC","SAEM","2PLC","2PNC","2PNI"), #"2PLI",
         lty=c(1,2,4,5,3),title="Error Estimation") #6,
  
  par(mar=c(4,5,1,1))
  RMSE = stats2D[(stats1D$J==J&stats1D$N==N),"A2_sd"]
  plot(XI[,2],XI[,2]-xi[,2],pch=19,cex=1,xlim=range(XI[,2]),ylim=range(c((XI[,2]-xi[,2])-2*RMSE,(XI[,2]-xi[,2])+2*RMSE)),
     main=NA,yaxt="n",xaxt="n",
     xlab="Generated Slope A2",ylab=expression(hat(A2)-A2 %+-% 2*sigma[RMSE]))
  axis(2,cex.axis=0.8,las=2)
  axis(1,cex.axis=0.8)
  arrows(XI[,2],(XI[,2]-xi[,2])-2*RMSE,XI[,2],(XI[,2]-xi[,2])+2*RMSE,code=3,angle=90,length=0.07)
  abline(h=0,lty=2)

  par(mar=c(4,2,1,1))
  plot(density(sqrt(EmpSE$VARLMIo[,2])/RMSE),yaxt="n",xaxt="n",ylab=NA,lty=3,
       xlim=c(0.25,2),main=NA,xlab=expression(hat(sigma)[A2]/sigma[RMSE]))
  mtext("Density Approximation",2,line=1,cex=0.7)
  axis(1,cex.axis=0.8)
  lines(density(sqrt(oError[,2])/RMSE),lty=5)
#  lines(density(sqrt(EmpSE$VARLMIi[,2])/RMSE),lty=6)
  lines(density(sqrt(iError[,2])/RMSE),lty=4)
  lines(density(sqrt(xiError[,2])/RMSE),lty=2)
  lines(density(sqrt(EmpSE$VARLMIo[,2])/RMSE),lty=3)
  lines(density(MCERR[,2]/RMSE),lty=1)
  legend("topright",legend = c("MCMC","SAEM","2PLC","2PNC","2PNI"), #"2PLI",
         lty=c(1,2,4,5,3),title="Error Estimation") #6,
  
  par(mar=c(4,5,1,1))
  RMSE = stats1D[(stats1D$J==J&stats1D$N==N),"B_sd"]
  plot(XI[,3],XI[,3]-xi[,3],pch=19,cex=1,xlim=range(XI[,3]),
     ylim=range(c((XI[,3]-xi[,3])-2*RMSE,(XI[,3]-xi[,3])+2*RMSE)),
     main=NA,yaxt="n",xaxt="n",
     xlab="Generated Intercept",ylab=expression(hat(b)-b %+-% 2*sigma[RMSE]))
  axis(2,cex.axis=0.8,las=2)
  axis(1,cex.axis=0.8)
  arrows(XI[,3],(XI[,3]-xi[,3])-2*RMSE,XI[,3],(XI[,3]-xi[,3])+2*RMSE,code=3,angle=90,length=0.07)
  abline(h=0,lty=2)

  par(mar=c(4,2,1,1))
  plot(density(sqrt(iError[,3])/RMSE),yaxt="n",xaxt="n",ylab=NA,lty=4,
       xlim=c(0.25,2),main=NA,xlab=expression(hat(sigma)[b]/sigma[RMSE]))
  mtext("Density Approximation",2,line=1,cex=0.7)
  axis(1,cex.axis=0.8)
  lines(density(sqrt(oError[,3])/RMSE),lty=5)
  #lines(density(sqrt(EmpSE$VARLMIi[,3])/RMSE),lty=6)
#  lines(density(sqrt(EmpSE$VARLMIo[,3])/RMSE),lty=3)
  lines(density(sqrt(xiError[,3])/RMSE),lty=2)
  lines(density(sqrt(EmpSE$VARLMIo[,3])/RMSE),lty=3)
  lines(density(MCERR[,3]/RMSE),lty=1)
  legend("topright",legend = c("MCMC","SAEM","2PLC","2PNC","2PNI"), #"2PLI",
         lty=c(1,2,4,5,3),title="Error Estimation") #6,
  
  dev.off()

  detach(FitDATA)
}
Plot2DRunErrors(J=50,N=5000,S=1,R=7,sMC = 300)
# colnames(stats2D)
# Tree2D<-stats2D
# Tree2D$LMI_NAN<-factor((is.nan(stats2D$B_sd_LMI)|is.nan(stats2D$A1_sd_LMI)|is.nan(stats2D$A2_sd_LMI))+0)
# print(rpart(LMI_NAN~.,data=Tree2D))
# stats2D[is.nan(stats2D$B_sd_LMI),c(4,24,44)]
# library(bnlearn)
# ?mmpc

#stats2D<-stats2D[!is.nan(stats2D$B_sd_LMI)&!is.nan(stats2D$A2_sd_LMI)&!is.nan(stats2D$A1_sd_LMI),]

pdf("ErrorFitStats2D.pdf",width=8,height=12)
CUTS = !is.nan(stats2D$B_sd_LMI)&!is.nan(stats2D$A2_sd_LMI)&!is.nan(stats2D$A1_sd_LMI)
CUTS1 = stats2D$B_sd_LMIo>0&stats2D$A2_sd_LMIo>0&stats2D$A1_sd_LMIo>0
nf <- layout(matrix(c(1,2,3,4,6,7,8,9,11,12,13,14,5,5,5,5,10,10,10,10,15,15,15,15),12,2), c(16,16), c(4,3,3,5,4,3,3,5,4,3,3,5), TRUE)
#layout.show(nf)
SIZE=10000
par(cex=0.6)
YR<-range(c(stats2D$A1_sd_oERR,stats2D$A1_sd_LMIo,stats2D$A1_sd_MCE)/stats2D$A1_sd)#,stats2D$A1_sd_LMI[CUTS]
fit2A1oE<-loess(A1_sd_oERR/A1_sd~A1,data=stats2D,span=0.5)
fit2A1oE
sd(fit2A1oE$residuals)
#fit2A1o<-lm(A1_sd_LMIo/A1_sd~A1,data=stats2D)
fit2A1o<-loess(A1_sd_LMIo/A1_sd~A1,data=stats2D[CUTS1,],span=0.5)
fit2A1o
sd(fit2A1o$residuals)
#fit2A1MC<-lm(A1_sd_MCE/A1_sd~A1,data=stats2D)
fit2A1MC<-loess(A1_sd_MCE/A1_sd~A1,data=stats2D,span=0.5)
fit2A1MC
sd(fit2A1MC$residuals)
#fit2A1LM<-lm(A1_sd_LMI/A1_sd~A1,data=stats2D)
fit2A1LM<-loess(A1_sd_LMI/A1_sd~A1,data=stats2D[CUTS,],span=0.5)
fit2A1LM
sd(fit2A1LM$residuals)

ticks <- pretty(YR)[-length(pretty(YR))]
labels <- format(ticks, big.mark=",", scientific=FALSE)

par(mar=c(m.b,4,2,m.b))
plot(stats2D$A1,stats2D$A1_sd_MCE/stats2D$A1_sd,pch=16,col=2,#cex=(SIZE-5000+stats2D$N)/SIZE,
     ylim=YR,xaxt="n",ylab=NA,yaxt="n")
axis(2,at=ticks,labels=labels,las=2)
text(range(stats2D$A1)[1]+0.7*(range(stats2D$A1)[2]-range(stats2D$A1)[1]),0.8*YR[2],"MCMC")
abline(h=1,lwd=2,lty=2)
lines(seq(range(stats2D$A1)[1],range(stats2D$A1)[2],length.out = 50),predict(fit2A1MC,newdata = data.frame(A1=seq(range(stats2D$A1)[1],range(stats2D$A1)[2],length.out = 50))),lty=4)
par(mar=c(m.b,4,m.b,m.b))
plot(stats2D$A1[CUTS],stats2D$A1_sd_LMI[CUTS]/stats2D$A1_sd[CUTS],
     yaxt="n",pch=16,col=4,#cex=(SIZE-5000+stats2D$N)/SIZE,
     ylim=YR,ylab=NA,xaxt="n")
axis(2,at=ticks,labels=labels,las=2)
text(range(stats2D$A1)[1]+0.7*(range(stats2D$A1)[2]-range(stats2D$A1)[1]),0.8*YR[2],"SAEM")
abline(h=1,lwd=2,lty=2)
mtext(text = "Ratio of Error Approximation to RMSE",side = 2,line = 2.5,cex = 0.6,adj = 1)
lines(seq(range(stats2D$A1)[1],range(stats2D$A1)[2],length.out = 50),predict(fit2A1LM,newdata = data.frame(A1=seq(range(stats2D$A1)[1],range(stats2D$A1)[2],length.out = 50))),lty=5)
plot(stats2D$A1,stats2D$A1_sd_oERR/stats2D$A1_sd,xaxt="n",pch=16,#cex=(SIZE-5000+stats2D$N)/SIZE,
     ylim=YR,ylab=NA,yaxt="n")#ylab="Ratio of Error Approximation to RMSE",xlab="Generated Slope (2D, 2PNO)")
axis(2,at=ticks,labels=labels,las=2)
text(range(stats2D$A1)[1]+0.7*(range(stats2D$A1)[2]-range(stats2D$A1)[1]),0.8*YR[2],"2PNC")
abline(h=1,lwd=2,lty=2)
lines(seq(range(stats2D$A1)[1],range(stats2D$A1)[2],length.out = 50),predict(fit2A1oE,newdata = data.frame(A1=seq(range(stats2D$A1)[1],range(stats2D$A1)[2],length.out = 50))),lty=1)
par(mar=c(4,4,m.b,m.b))
plot(stats2D$A1[CUTS1],stats2D$A1_sd_LMIo[CUTS1]/stats2D$A1_sd[CUTS1],pch=16,#cex=(SIZE-5000+stats2D$N)/SIZE,
     col=3,ylim=YR,ylab=NA,yaxt="n",xlab="Generated Slope (2D, 2PNO)")
axis(2,at=ticks,labels=labels,las=2)
text(range(stats2D$A1)[1]+0.7*(range(stats2D$A1)[2]-range(stats2D$A1)[1]),0.8*YR[2],"2PNI")
abline(h=1,lwd=2,lty=2)
lines(seq(range(stats2D$A1)[1],range(stats2D$A1)[2],length.out = 50),predict(fit2A1o,newdata = data.frame(A1=seq(range(stats2D$A1)[1],range(stats2D$A1)[2],length.out = 50))),lty=3)
#fit2A1oE<-lm(A1_sd_oERR/A1_sd~A1,data=stats2D)
# legend("topleft",c("2PNO (At Convergence)","2PNO (Post Convergence)","Iterative","MCMC (Post Convergence)","N = 2000","N = 5000","N = 8000"),
#        pch=c(16,16,3,2,1,1,1),col=c(1,3,4,2,1,1,1),pt.cex=c(1,1,1,1,(SIZE-5000+c(2000,5000,8000))/SIZE))
par(mar=c(4,5,2,m.b))
plot(stats2D$A1,stats2D$A1_sd_oERR/stats2D$A1_sd,ylab="Ratio of Error Approximation to RMSE",xlab="Generated Slope (2D, 2PNO)",ylim=c(0,2),type="n")
lines(seq(range(stats2D$A1)[1],range(stats2D$A1)[2],length.out = 50),predict(fit2A1oE,newdata = data.frame(A1=seq(range(stats2D$A1)[1],range(stats2D$A1)[2],length.out = 50))),lty=1)
lines(seq(range(stats2D$A1)[1],range(stats2D$A1)[2],length.out = 50),predict(fit2A1o,newdata = data.frame(A1=seq(range(stats2D$A1)[1],range(stats2D$A1)[2],length.out = 50))),lty=3)
lines(seq(range(stats2D$A1)[1],range(stats2D$A1)[2],length.out = 50),predict(fit2A1MC,newdata = data.frame(A1=seq(range(stats2D$A1)[1],range(stats2D$A1)[2],length.out = 50))),lty=4)
lines(seq(range(stats2D$A1)[1],range(stats2D$A1)[2],length.out = 50),predict(fit2A1LM,newdata = data.frame(A1=seq(range(stats2D$A1)[1],range(stats2D$A1)[2],length.out = 50))),lty=5)
legend("topright",c("MCMC","SAEM","2PNC","2PNI"),inset = 0.05,
       title=expression(paste("LOESS FIT TO ",hat(sigma)/sigma[RMSE])),lty=c(4,5,1,3))

CUTS2 = stats2D$A2>0
YR<-range(c(stats2D$A2_sd_oERR,stats2D$A2_sd_LMIo,stats2D$A2_sd_MCE)/stats2D$A2_sd)
fit2A2oE<-loess(A2_sd_oERR/A2_sd~A2,data=stats2D[CUTS2,],span=0.51)
fit2A2oE
sd(fit2A2oE$residuals)
#fit2A2o<-lm(A2_sd_LMIo/A2_sd~A2,data=stats2D)
fit2A2o<-loess(A2_sd_LMIo/A2_sd~A2,data=stats2D[CUTS2&CUTS1,],span=0.51)
fit2A2o
sd(fit2A2o$residuals)
#fit2A2MC<-lm(A2_sd_MCE/A2_sd~A2,data=stats2D)
fit2A2MC<-loess(A2_sd_MCE/A2_sd~A2,data=stats2D[CUTS2,],span=0.51)
fit2A2MC
sd(fit2A2MC$residuals)
#fit2A2LM<-lm(A2_sd_LMI/A2_sd~A2,data=stats2D)
fit2A2LM<-loess(A2_sd_LMI/A2_sd~A2,data=stats2D[CUTS&CUTS2,],span=0.51)
fit2A2LM
sd(fit2A2LM$residuals)
ticks <- pretty(YR)[-c(length(pretty(YR))-1:0)]
labels <- format(ticks, big.mark=",", scientific=FALSE)

par(mar=c(m.b,4,2,m.b))
plot(stats2D$A2,stats2D$A2_sd_MCE/stats2D$A2_sd,pch=16,col=2,#cex=(SIZE-5000+stats2D$N)/SIZE,
     ylim=YR,xaxt="n",ylab=NA,yaxt="n")
text(range(stats2D$A2)[1]+0.7*(range(stats2D$A2)[2]-range(stats2D$A2)[1]),0.8*YR[2],"MCMC")
axis(2,at=ticks,labels=labels,las=2)
abline(h=1,lwd=2,lty=2)
lines(seq(range(stats2D$A2)[1],range(stats2D$A2)[2],length.out = 50),predict(fit2A2MC,newdata = data.frame(A2=seq(range(stats2D$A2)[1],range(stats2D$A2)[2],length.out = 50))),lty=4)
par(mar=c(m.b,4,m.b,m.b))
plot(stats2D$A2[CUTS],stats2D$A2_sd_LMI[CUTS]/stats2D$A2_sd[CUTS],
     pch=16,col=4,#cex=(SIZE-5000+stats2D$N)/SIZE,
     ylim=YR,ylab=NA,yaxt="n",xaxt="n")
text(range(stats2D$A2)[1]+0.7*(range(stats2D$A2)[2]-range(stats2D$A2)[1]),0.8*YR[2],"SAEM")
axis(2,at=ticks,labels=labels,las=2)
abline(h=1,lwd=2,lty=2)
lines(seq(range(stats2D$A2)[1],range(stats2D$A2)[2],length.out = 50),predict(fit2A2LM,newdata = data.frame(A2=seq(range(stats2D$A2)[1],range(stats2D$A2)[2],length.out = 50))),lty=5)
mtext(text = "Ratio of Error Approximation to RMSE",side = 2,line = 2.5,cex = 0.6,adj = 1)
plot(stats2D$A2,stats2D$A2_sd_oERR/stats2D$A2_sd,xaxt="n",pch=16,#cex=(SIZE-5000+stats2D$N)/SIZE,
     ylim=YR,ylab=NA,yaxt="n")#ylab="Ratio of Error Approximation to RMSE",xlab="Generated Slope (2D, 2PNO)")
axis(2,at=ticks,labels=labels,las=2)
text(range(stats2D$A2)[1]+0.7*(range(stats2D$A2)[2]-range(stats2D$A2)[1]),0.8*YR[2],"2PNC")
abline(h=1,lwd=2,lty=2)
lines(seq(range(stats2D$A2)[1],range(stats2D$A2)[2],length.out = 50),predict(fit2A2oE,newdata = data.frame(A2=seq(range(stats2D$A2)[1],range(stats2D$A2)[2],length.out = 50))),lty=1)
par(mar=c(4,4,m.b,m.b))
plot(stats2D$A2[CUTS1],stats2D$A2_sd_LMIo[CUTS1]/stats2D$A2_sd[CUTS1],pch=16,#cex=(SIZE-5000+stats2D$N)/SIZE,
     col=3,ylim=YR,ylab=NA,yaxt="n",xlab="Generated Slope (2D, 2PNO)")
text(range(stats2D$A2)[1]+0.7*(range(stats2D$A2)[2]-range(stats2D$A2)[1]),0.8*YR[2],"2PNI")
axis(2,at=ticks,labels=labels,las=2)
abline(h=1,lwd=2,lty=2)
lines(seq(range(stats2D$A2)[1],range(stats2D$A2)[2],length.out = 50),predict(fit2A2o,newdata = data.frame(A2=seq(range(stats2D$A2)[1],range(stats2D$A2)[2],length.out = 50))),lty=3)
#fit2A2oE<-lm(A2_sd_oERR/A2_sd~A2,data=stats2D)
# legend("topleft",c("2PNO (At Convergence)","2PNO (Post Convergence)","Iterative","MCMC (Post Convergence)","N = 2000","N = 5000","N = 8000"),
#        pch=c(16,16,3,2,1,1,1),col=c(1,3,4,2,1,1,1),pt.cex=c(1,1,1,1,(SIZE-5000+c(2000,5000,8000))/SIZE))
par(mar=c(4,5,2,m.b))
plot(stats2D$A2,stats2D$A2_sd_oERR/stats2D$A2_sd,ylab="Ratio of Error Approximation to RMSE",xlab="Generated Slope (2D, 2PNO)",ylim=YR,type="n")
lines(seq(range(stats2D$A2)[1],range(stats2D$A2)[2],length.out = 50),predict(fit2A2oE,newdata = data.frame(A2=seq(range(stats2D$A2)[1],range(stats2D$A2)[2],length.out = 50))),lty=1)
lines(seq(range(stats2D$A2)[1],range(stats2D$A2)[2],length.out = 50),predict(fit2A2o,newdata = data.frame(A2=seq(range(stats2D$A2)[1],range(stats2D$A2)[2],length.out = 50))),lty=3)
lines(seq(range(stats2D$A2)[1],range(stats2D$A2)[2],length.out = 50),predict(fit2A2MC,newdata = data.frame(A2=seq(range(stats2D$A2)[1],range(stats2D$A2)[2],length.out = 50))),lty=4)
lines(seq(range(stats2D$A2)[1],range(stats2D$A2)[2],length.out = 50),predict(fit2A2LM,newdata = data.frame(A2=seq(range(stats2D$A2)[1],range(stats2D$A2)[2],length.out = 50))),lty=5)
legend("topright",c("MCMC","SAEM","2PNC","2PNI"),inset = 0.05,
       title=expression(paste("LOESS FIT TO ",hat(sigma)/sigma[RMSE])),lty=c(4,5,1,3))

#CUTS = !is.nan(stats2D$B_sd_LMI)&!is.nan(stats2D$A2_sd_LMI)&!is.nan(stats2D$A1_sd_LMI)
CUTS3 = stats2D$B_sd_LMIo/stats2D$B_sd<5
YR<-c(0,max(c(stats2D$B_sd_oERR,stats2D$B_sd_LMIo[CUTS3],stats2D$B_sd_MCE)/c(stats2D$B_sd,stats2D$B_sd[CUTS3],stats2D$B_sd)))
#stats2D$B2<-stats2D$B^2

#fit2BoE<-lm(B_sd_oERR/B_sd~B+B2,data=stats2D)
fit2BoE<-loess(B_sd_oERR/B_sd~B,data=stats2D,span=0.5)
fit2BoE
sd(fit2BoE$residuals)
#fit2Bo<-lm(B_sd_LMIo/B_sd~B+B2,data=stats2D)

fit2Bo<-loess(B_sd_LMIo/B_sd~B,data=stats2D[CUTS3&CUTS1,],span=0.5)
fit2Bo
sd(fit2Bo$residuals)
#fit2BMC<-lm(B_sd_MCE/B_sd~B+B2,data=stats2D)
fit2BMC<-loess(B_sd_MCE/B_sd~B,data=stats2D,span=0.5)
fit2BMC
sd(fit2BMC$residuals)
#fit2BLM<-lm(B_sd_LMI/B_sd~B+B2,data=stats2D)
fit2BLM<-loess(B_sd_LMI/B_sd~B,data=stats2D[CUTS,],span=0.5)
fit2BLM
sd(fit2BLM$residuals)
ticks <- pretty(YR)[-c(length(pretty(YR))-1:0)]
labels <- format(ticks, big.mark=",", scientific=FALSE)
par(mar=c(m.b,4,2,m.b))
plot(stats2D$B,stats2D$B_sd_MCE/stats2D$B_sd,pch=16,col=2,#cex=(SIZE-5000+stats2D$N)/SIZE,
     ylim=YR,xaxt="n",ylab=NA,yaxt="n")
text(range(stats2D$B)[1]+0.7*(range(stats2D$B)[2]-range(stats2D$B)[1]),0.8*YR[2],"MCMC")
axis(2,at=ticks,labels=labels,las=2)
abline(h=1,lwd=2,lty=2)
lines(seq(range(stats2D$B)[1],range(stats2D$B)[2],length.out=50),predict(fit2BMC,newdata = data.frame(B=seq(range(stats2D$B)[1],range(stats2D$B)[2],length.out=50),B2=seq(range(stats2D$B)[1],range(stats2D$B)[2],length.out=50)^2)),lty=4)
par(mar=c(m.b,4,m.b,m.b))
plot(stats2D$B[CUTS],stats2D$B_sd_LMI[CUTS]/stats2D$B_sd[CUTS],
     xaxt="n", pch=16,col=4,#cex=(SIZE-5000+stats2D$N)/SIZE,
     ylim=YR,ylab=NA,yaxt="n")
text(range(stats2D$B)[1]+0.7*(range(stats2D$B)[2]-range(stats2D$B)[1]),0.8*YR[2],"SAEM")
axis(2,at=ticks,labels=labels,las=2)
abline(h=1,lwd=2,lty=2)
lines(seq(range(stats2D$B)[1],range(stats2D$B)[2],length.out=50),predict(fit2BLM,newdata = data.frame(B=seq(range(stats2D$B)[1],range(stats2D$B)[2],length.out=50),B2=seq(range(stats2D$B)[1],range(stats2D$B)[2],length.out=50)^2)),lty=5)
mtext(text = "Ratio of Error Approximation to RMSE",side = 2,line = 2.5,cex = 0.6,adj = 1)
plot(stats2D$B,stats2D$B_sd_oERR/stats2D$B_sd,pch=16,#cex=(SIZE-5000+stats2D$N)/SIZE,
     ylim=YR,ylab=NA,xaxt="n",yaxt="n")#,ylab="Ratio of Error Approximation to RMSE",xlab="Generated Intercepts (2D, 2PNO)")
axis(2,at=ticks,labels=labels,las=2)
text(range(stats2D$B)[1]+0.7*(range(stats2D$B)[2]-range(stats2D$B)[1]),0.8*YR[2],"2PNC")
abline(h=1,lwd=2,lty=2)
lines(seq(range(stats2D$B)[1],range(stats2D$B)[2],length.out=50),predict(fit2BoE,newdata = data.frame(B=seq(range(stats2D$B)[1],range(stats2D$B)[2],length.out=50),B2=seq(range(stats2D$B)[1],range(stats2D$B)[2],length.out=50)^2)),lty=1)
par(mar=c(4,4,m.b,m.b))
plot(stats2D$B[CUTS1],stats2D$B_sd_LMIo[CUTS1]/stats2D$B_sd[CUTS1],pch=16,#cex=(SIZE-5000+stats2D$N)/SIZE,
     col=3,ylim=YR,ylab=NA,yaxt="n",xlab="Generated Intercepts (2D, 2PNO)")
text(range(stats2D$B)[1]+0.7*(range(stats2D$B)[2]-range(stats2D$B)[1]),0.8*YR[2],"2PNI")
axis(2,at=ticks,labels=labels,las=2)
abline(h=1,lwd=2,lty=2)
lines(seq(range(stats2D$B)[1],range(stats2D$B)[2],length.out=50),predict(fit2Bo,newdata = data.frame(B=seq(range(stats2D$B)[1],range(stats2D$B)[2],length.out=50),B2=seq(range(stats2D$B)[1],range(stats2D$B)[2],length.out=50)^2)),lty=3)
# legend("topleft",c("2PNO (At Convergence)","2PNO (Post Convergence)","Iterative","MCMC (Post Convergence)","N = 2000","N = 5000","N = 8000"),
#        pch=c(16,16,3,2,1,1,1),col=c(1,3,4,2,1,1,1),pt.cex=c(1,1,1,1,(SIZE-5000+c(2000,5000,8000))/SIZE))
par(mar=c(4,5,2,m.b))
plot(stats2D$B,stats2D$B_sd_oERR/stats2D$B_sd,ylab="Ratio of Error Approximation to RMSE",xlab="Generated Intercepts (2D, 2PNO)",ylim=c(0,2),type="n")
lines(seq(range(stats2D$B)[1],range(stats2D$B)[2],length.out=50),predict(fit2BoE,newdata = data.frame(B=seq(range(stats2D$B)[1],range(stats2D$B)[2],length.out=50),B2=seq(range(stats2D$B)[1],range(stats2D$B)[2],length.out=50)^2)),lty=1)
lines(seq(range(stats2D$B)[1],range(stats2D$B)[2],length.out=50),predict(fit2Bo,newdata = data.frame(B=seq(range(stats2D$B)[1],range(stats2D$B)[2],length.out=50),B2=seq(range(stats2D$B)[1],range(stats2D$B)[2],length.out=50)^2)),lty=3)
lines(seq(range(stats2D$B)[1],range(stats2D$B)[2],length.out=50),predict(fit2BMC,newdata = data.frame(B=seq(range(stats2D$B)[1],range(stats2D$B)[2],length.out=50),B2=seq(range(stats2D$B)[1],range(stats2D$B)[2],length.out=50)^2)),lty=4)
lines(seq(range(stats2D$B)[1],range(stats2D$B)[2],length.out=50),predict(fit2BLM,newdata = data.frame(B=seq(range(stats2D$B)[1],range(stats2D$B)[2],length.out=50),B2=seq(range(stats2D$B)[1],range(stats2D$B)[2],length.out=50)^2)),lty=5)
legend("topright",c("MCMC","SAEM","2PNC","2PNI"),inset = 0.05,
       title=expression(paste("LOESS FIT TO ",hat(sigma)/sigma[RMSE])),lty=c(4,5,1,3))



#par(mfrow=c(3,2),mar=c(5,5,1,2))
# SIZE=10000
# YR<-range(c(stats2D$A1_sd_oERR,stats2D$A1_sd_LMIo,stats2D$A1_sd_MCE,stats2D$A1_sd_LMI)/stats2D$A1_sd)
# plot(stats2D$A1,stats2D$A1_sd_oERR/stats2D$A1_sd,pch=16,cex=(SIZE-5000+stats2D$N)/SIZE,ylim=YR,ylab="Ratio of Error Approximation to RMSE",xlab="Generated Slope 1 (2D, 2PNO)")
# points(stats2D$A1,stats2D$A1_sd_LMIo/stats2D$A1_sd,pch=16,cex=(SIZE-5000+stats2D$N)/SIZE,col=3)
# points(stats2D$A1,stats2D$A1_sd_MCE/stats2D$A1_sd,pch=2,col=2,cex=(SIZE-5000+stats2D$N)/SIZE)
# points(stats2D$A1,stats2D$A1_sd_LMI/stats2D$A1_sd,pch=3,col=4,cex=(SIZE-5000+stats2D$N)/SIZE)
# abline(h=1,lwd=2,lty=2)
# fit2A1oE<-lm(A1_sd_oERR/A1_sd~A1,data=stats2D)
# fit2A1oE
# sd(fit2A1oE$residuals)
# fit2A1o<-lm(A1_sd_LMIo/A1_sd~A1,data=stats2D)
# fit2A1o
# sd(fit2A1o$residuals)
# fit2A1MC<-lm(A1_sd_MCE/A1_sd~A1,data=stats2D)
# fit2A1MC
# sd(fit2A1MC$residuals)
# fit2A1LM<-lm(A1_sd_LMI/A1_sd~A1,data=stats2D)
# fit2A1LM
# sd(fit2A1LM$residuals)
# lines(range(stats2D$A1),predict.lm(fit2A1oE,newdata = data.frame(A1=range(stats2D$A1))),lty=1)
# lines(range(stats2D$A1),predict.lm(fit2A1o,newdata = data.frame(A1=range(stats2D$A1))),lty=3)
# lines(range(stats2D$A1),predict.lm(fit2A1MC,newdata = data.frame(A1=range(stats2D$A1))),lty=4)
# lines(range(stats2D$A1),predict.lm(fit2A1LM,newdata = data.frame(A1=range(stats2D$A1))),lty=5)
# legend("topleft",c("2PNO (At Convergence)","2PNO (Post Convergence)","Iterative","MCMC (Post Convergence)","N = 2000","N = 5000","N = 8000"),
#        pch=c(16,16,3,2,1,1,1),col=c(1,3,4,2,1,1,1),pt.cex=c(1,1,1,1,(SIZE-5000+c(2000,5000,8000))/SIZE))
# plot(stats2D$A1,stats2D$A1_sd_oERR/stats2D$A1_sd,ylab="Ratio of Error Approximation to RMSE",xlab="Generated Slope 1 (2D, 2PNO)",ylim=YR,type="n")
# lines(range(stats2D$A1),predict.lm(fit2A1oE,newdata = data.frame(A1=range(stats2D$A1))),lty=1)
# lines(range(stats2D$A1),predict.lm(fit2A1o,newdata = data.frame(A1=range(stats2D$A1))),lty=3)
# lines(range(stats2D$A1),predict.lm(fit2A1MC,newdata = data.frame(A1=range(stats2D$A1))),lty=4)
# lines(range(stats2D$A1),predict.lm(fit2A1LM,newdata = data.frame(A1=range(stats2D$A1))),lty=5)
# legend("topright",c(expression(paste("2PNO (At Convergence): ",beta[0] == 1.17,"  ", beta[1] == -0.44)),
#                     expression(paste("2PNO (Post Convergence): ",beta[0] == 1.26,"  ", beta[1] == -0.23)),
#                     expression(paste("Iterative: ",beta[0] == 1.11,"  ", beta[1] == -0.34)),
#                     expression(paste("MCMC: ",beta[0] == 0.44,"  ", beta[1] == 0.65))),lty=c(1,3,5,4))
# 
# 
# YR<-range(c(stats2D$A2_sd_oERR,stats2D$A2_sd_LMIo,stats2D$A2_sd_MCE,stats2D$A2_sd_LMI)/stats2D$A2_sd)
# plot(stats2D$A2,stats2D$A2_sd_oERR/stats2D$A2_sd,pch=16,cex=(SIZE-5000+stats2D$N)/SIZE,ylim=YR,ylab="Ratio of Error Approximation to RMSE",xlab="Generated Slope 2 (2D, 2PNO)")
# points(stats2D$A2,stats2D$A2_sd_LMIo/stats2D$A2_sd,pch=16,cex=(SIZE-5000+stats2D$N)/SIZE,col=3)
# points(stats2D$A2,stats2D$A2_sd_MCE/stats2D$A2_sd,pch=2,col=2,cex=(SIZE-5000+stats2D$N)/SIZE)
# points(stats2D$A2,stats2D$A2_sd_LMI/stats2D$A2_sd,pch=3,col=4,cex=(SIZE-5000+stats2D$N)/SIZE)
# abline(h=1,lwd=2,lty=2)
# fit2A2oE<-lm(A2_sd_oERR/A2_sd~A2,data=stats2D)
# fit2A2oE
# sd(fit2A2oE$residuals)
# fit2A2o<-lm(A2_sd_LMIo/A2_sd~A2,data=stats2D)
# fit2A2o
# sd(fit2A2o$residuals)
# fit2A2MC<-lm(A2_sd_MCE/A2_sd~A2,data=stats2D)
# fit2A2MC
# sd(fit2A2MC$residuals)
# fit2A2LM<-lm(A2_sd_LMI/A2_sd~A2,data=stats2D)
# fit2A2LM
# sd(fit2A2LM$residuals)
# lines(range(stats2D$A2),predict.lm(fit2A2oE,newdata = data.frame(A2=range(stats2D$A2))),lty=1)
# lines(range(stats2D$A2),predict.lm(fit2A2o,newdata = data.frame(A2=range(stats2D$A2))),lty=3)
# lines(range(stats2D$A2),predict.lm(fit2A2MC,newdata = data.frame(A2=range(stats2D$A2))),lty=4)
# lines(range(stats2D$A2),predict.lm(fit2A2LM,newdata = data.frame(A2=range(stats2D$A2))),lty=5)
# legend("topleft",c("2PNO (At Convergence)","2PNO (Post Convergence)","Iterative","MCMC (Post Convergence)","N = 2000","N = 5000","N = 8000"),
#        pch=c(16,16,3,2,1,1,1),col=c(1,3,4,2,1,1,1),pt.cex=c(1,1,1,1,(SIZE-5000+c(2000,5000,8000))/SIZE))
# plot(stats2D$A2,stats2D$A2_sd_oERR/stats2D$A2_sd,ylab="Ratio of Error Approximation to RMSE",xlab="Generated Slope 2 (2D, 2PNO)",ylim=YR,type="n")
# lines(range(stats2D$A2),predict.lm(fit2A2oE,newdata = data.frame(A2=range(stats2D$A2))),lty=1)
# lines(range(stats2D$A2),predict.lm(fit2A2o,newdata = data.frame(A2=range(stats2D$A2))),lty=3)
# lines(range(stats2D$A2),predict.lm(fit2A2MC,newdata = data.frame(A2=range(stats2D$A2))),lty=4)
# lines(range(stats2D$A2),predict.lm(fit2A2LM,newdata = data.frame(A2=range(stats2D$A2))),lty=5)
# legend("topright",c(expression(paste("2PNO (At Convergence): ",beta[0] == 1.17,"  ", beta[1] == -0.44)),
#                     expression(paste("2PNO (Post Convergence): ",beta[0] == 1.26,"  ", beta[1] == -0.23)),
#                     expression(paste("Iterative: ",beta[0] == 1.11,"  ", beta[1] == -0.34)),
#                     expression(paste("MCMC: ",beta[0] == 0.44,"  ", beta[1] == 0.65))),lty=c(1,3,5,4))
# 
# YR<-range(c(stats2D$B_sd_oERR,stats2D$B_sd_LMIo,stats2D$B_sd_MCE,stats2D$B_sd_LMI)/stats2D$B_sd)
# plot(stats2D$B,stats2D$B_sd_oERR/stats2D$B_sd,pch=16,cex=(SIZE-5000+stats2D$N)/SIZE,ylim=YR,ylab="Ratio of Error Approximation to RMSE",xlab="Generated Intercepts (2D, 2PNO)")
# points(stats2D$B,stats2D$B_sd_LMIo/stats2D$B_sd,pch=16,cex=(SIZE-5000+stats2D$N)/SIZE,col=3)
# points(stats2D$B,stats2D$B_sd_MCE/stats2D$B_sd,pch=2,col=2,cex=(SIZE-5000+stats2D$N)/SIZE)
# points(stats2D$B,stats2D$B_sd_LMI/stats2D$B_sd,pch=3,col=4,cex=(SIZE-5000+stats2D$N)/SIZE)
# abline(h=1,lwd=2,lty=2)
# stats2D$B2<-stats2D$B^2
# fit2BoE<-lm(B_sd_oERR/B_sd~B+B2,data=stats2D)
# fit2BoE
# sd(fit2BoE$residuals)
# fit2Bo<-lm(B_sd_LMIo/B_sd~B+B2,data=stats2D)
# fit2Bo
# sd(fit2Bo$residuals)
# fit2BMC<-lm(B_sd_MCE/B_sd~B+B2,data=stats2D)
# fit2BMC
# sd(fit2BMC$residuals)
# fit2BLM<-lm(B_sd_LMI/B_sd~B+B2,data=stats2D)
# fit2BLM
# sd(fit2BLM$residuals)
# lines(seq(range(stats2D$B)[1],range(stats2D$B)[2],length.out=50),predict.lm(fit2BoE,newdata = data.frame(B=seq(range(stats2D$B)[1],range(stats2D$B)[2],length.out=50),B2=seq(range(stats2D$B)[1],range(stats2D$B)[2],length.out=50)^2)),lty=1)
# lines(seq(range(stats2D$B)[1],range(stats2D$B)[2],length.out=50),predict.lm(fit2Bo,newdata = data.frame(B=seq(range(stats2D$B)[1],range(stats2D$B)[2],length.out=50),B2=seq(range(stats2D$B)[1],range(stats2D$B)[2],length.out=50)^2)),lty=3)
# lines(seq(range(stats2D$B)[1],range(stats2D$B)[2],length.out=50),predict.lm(fit2BMC,newdata = data.frame(B=seq(range(stats2D$B)[1],range(stats2D$B)[2],length.out=50),B2=seq(range(stats2D$B)[1],range(stats2D$B)[2],length.out=50)^2)),lty=4)
# lines(seq(range(stats2D$B)[1],range(stats2D$B)[2],length.out=50),predict.lm(fit2BLM,newdata = data.frame(B=seq(range(stats2D$B)[1],range(stats2D$B)[2],length.out=50),B2=seq(range(stats2D$B)[1],range(stats2D$B)[2],length.out=50)^2)),lty=5)
# legend("topleft",c("2PNO (At Convergence)","2PNO (Post Convergence)","Iterative","MCMC (Post Convergence)","N = 2000","N = 5000","N = 8000"),
#        pch=c(16,16,3,2,1,1,1),col=c(1,3,4,2,1,1,1),pt.cex=c(1,1,1,1,(SIZE-5000+c(2000,5000,8000))/SIZE))
# plot(stats2D$B,stats2D$B_sd_oERR/stats2D$B_sd,ylab="Ratio of Error Approximation to RMSE",xlab="Generated Intercepts (2D, 2PNO)",ylim=YR,type="n")
# lines(seq(range(stats2D$B)[1],range(stats2D$B)[2],length.out=50),predict.lm(fit2BoE,newdata = data.frame(B=seq(range(stats2D$B)[1],range(stats2D$B)[2],length.out=50),B2=seq(range(stats2D$B)[1],range(stats2D$B)[2],length.out=50)^2)),lty=1)
# lines(seq(range(stats2D$B)[1],range(stats2D$B)[2],length.out=50),predict.lm(fit2Bo,newdata = data.frame(B=seq(range(stats2D$B)[1],range(stats2D$B)[2],length.out=50),B2=seq(range(stats2D$B)[1],range(stats2D$B)[2],length.out=50)^2)),lty=3)
# lines(seq(range(stats2D$B)[1],range(stats2D$B)[2],length.out=50),predict.lm(fit2BMC,newdata = data.frame(B=seq(range(stats2D$B)[1],range(stats2D$B)[2],length.out=50),B2=seq(range(stats2D$B)[1],range(stats2D$B)[2],length.out=50)^2)),lty=4)
# lines(seq(range(stats2D$B)[1],range(stats2D$B)[2],length.out=50),predict.lm(fit2BLM,newdata = data.frame(B=seq(range(stats2D$B)[1],range(stats2D$B)[2],length.out=50),B2=seq(range(stats2D$B)[1],range(stats2D$B)[2],length.out=50)^2)),lty=5)
# legend("topright",c(expression(paste("2PNO (At Convergence): ",beta[0] == 1.14,"  ", beta[1] == 0.012,"  ", beta[2] == -0.119)),
#                     expression(paste("2PNO (Post Convergence): ",beta[0] == 1.10,"  ", beta[1] == 0.023,"  ", beta[2] == -0.035)),
#                     expression(paste("Iterative: ",beta[0] == 1.28,"  ", beta[1] == 0.001,"  ", beta[2] == -0.123)),
#                     expression(paste("MCMC: ",beta[0] == 0.79,"  ", beta[1] == 0.044,"  ", beta[2] == 0.41))),lty=c(1,3,5,4))

# 
# #for (b in 1:length(burn)) {
#   for (j in 1:length(items)) {
#     for (n in 1:length(examinees)) {
#       simfile=paste("Error-",100,"/SimError2D-",items,"-",examinees,".rda",sep="")    # Save to an ".rda" file with this name
#       load(simfile)
#       estfile=paste("Error-",100,"/FitError2D-",items,"-",examinees,"-",i,".rda",sep="") 
#       load(estfile)
#       XIiter<-array(FitDATA$xi, dim=c(items,1+settings$Adim,sims))        
#       ERRiter<-array(FitDATA$xiError, dim=c(items,1+settings$Adim,sims))        
#       for (i in 1:sims) {
#         estfile=paste("Error-",100,"/FitError2D-",items,"-",examinees,"-",i,".rda",sep="") 
#         load(estfile)
#         XIiter[,,i]<-as.matrix(FitDATA$xi) 
#         ERRiter[,,i]<-as.matrix(FitDATA$xiError) 
#       }
#       mXI<-apply(XIiter,c(1,2),mean)
#       mLMI<-apply(ERRiter,c(1,2),mean)
#       bias<-mXI-gen.xi
#       print(dim(XIiter))
#       sXI<-apply(XIiter-array(rep(mXI,sims),dim=c(items,1+settings$Adim,sims)),c(1,2),sd)
#       sLMI<-apply(ERRiter-array(rep(mLMI,sims),dim=c(items,1+settings$Adim,sims)),c(1,2),sd)
#       par(mfrow=c(ncol(gen.xi),1))
#       for (p in 1:ncol(gen.xi)) {
#         ifelse(p==ncol(gen.xi),pt<-"B",pt<-paste("A",p,sep=""))
#         plot(gen.xi[,p],bias[,p],pch=19,cex=1.5,xlim=c(0.9,1.1)*range(gen.xi[,p]),ylim=range(c(bias[,p]-2*sXI[,p],bias[,p]+2*sXI[,p])),
#              main=paste("Parameter",pt,"Bias & Error (",sims," Runs)\n",items,"items; ",examinees,"examinees; ",burn[b],"burn-in"),
#              xlab=pt,ylab=paste(pt,"hat"))
#         arrows(gen.xi[,p],bias[,p]-2*sXI[,p],gen.xi[,p],bias[,p]+2*sXI[,p],code=3,angle=90,length=0.07)
#       }
#       if (b==1&j==1&n==1) {
#         stats2D<-data.frame(N=rep(examinees,items),J=rep(items,items),Burn=rep(100,items))
#         for(p in 1:ncol(gen.xi)) {
#           ifelse(p==ncol(gen.xi),pt<-"B",pt<-paste("A",p,sep=""))
#           stats2D<-cbind(stats2D,gen.xi[,p])
#           colnames(stats2D)[ncol(stats2D)]<-pt
#           stats2D<-cbind(stats2D,bias[,p])
#           colnames(stats2D)[ncol(stats2D)]<-paste(pt,"bias",sep="_")
#           stats2D<-cbind(stats2D,sXI[,p])
#           colnames(stats2D)[ncol(stats2D)]<-paste(pt,"sd",sep="_")
#           stats2D<-cbind(stats2D,mLMI[,p])
#           colnames(stats2D)[ncol(stats2D)]<-paste(pt,"VAR_LMI",sep="_")
#           stats2D<-cbind(stats2D,sqrt(mLMI[,p]))
#           colnames(stats2D)[ncol(stats2D)]<-paste(pt,"sd_LMI",sep="_")
#           stats2D<-cbind(stats2D,sLMI[,p])
#           colnames(stats2D)[ncol(stats2D)]<-paste(pt,"sd_VAR_LMI",sep="_")
#         }
#       } else {
#         stdat<-cbind(rep(examinees,items),rep(items,items),rep(100,items))
#         for(p in 1:ncol(gen.xi)) {
#           stdat<-cbind(stdat,gen.xi[,p],bias[,p],sXI[,p],mLMI[,p],sqrt(mLMI[,p]),sLMI[,p])
#         }
#         colnames(stdat)<-colnames(stats2D)
#         stats2D<-rbind(stats2D,stdat)
#       }
#     }
#   }
# #}
# 
# par(mfrow=c(1,2))
# plot(stats[stats$N==4000&stats$J==50,c("A1","A1_sd")],ylab="Empirical Standard Error",xlab="Slope Parameter",
#      pch=19,ylim=c(0,max(stats[stats$N==4000&stats$J==50,"A1_sd"])),xlim=c(0,1.8),
#      main="Comparison of Standard Errors\nof Slope Parameters")
# abline(lm(stats[stats$N==4000&stats$J==50,c("A1_sd","A1")])$coefficients,lwd=1.5)
# points(STATS[STATS$N==4000&STATS$J==50,c("A1","A1_sd")],pch=19,col=2)
# abline(lm(STATS[STATS$N==4000&STATS$J==50,c("A1_sd","A1")])$coefficients,lwd=1.5,col=2)
# points(STATS[STATS$N==4000&STATS$J==50,c("A2","A2_sd")],pch=19,col=4)
# abline(lm(STATS[STATS$N==4000&STATS$J==50,c("A2_sd","A2")])$coefficients,lwd=1.5,col=4)
# legend("topleft",c("Slope-1D","Slope1-2D","Slope2-2D"),pch=19,lty=1,col=c(1,2,4))
# plot(stats[stats$N==4000&stats$J==50,c("B","B_sd")],ylab="Empirical Standard Error",xlab="Slope Parameter",
#      pch=19,ylim=c(0,max(stats[stats$N==4000&stats$J==50,"B_sd"])),xlim=1.1*range(c(stats[stats$N==4000&stats$J==50,"B"],STATS[STATS$N==4000&STATS$J==50,"B"])),
#      main="Comparison of Standard Errors\nof Intercepts")
# lines(stats[stats$N==4000&stats$J==50,"B"][order(stats[stats$N==4000&stats$J==50,"B"])],predict(lm(stats[stats$N==4000&stats$J==50,"B_sd"]~stats[stats$N==4000&stats$J==50,"B"]+I(stats[stats$N==4000&stats$J==50,"B"]^2)))[order(stats[stats$N==4000&stats$J==50,"B"])],lwd=1.5)
# points(STATS[STATS$N==4000&STATS$J==50,c("B","B_sd")],pch=19,col=2)
# lines(STATS[STATS$N==4000&STATS$J==50,"B"][order(STATS[STATS$N==4000&STATS$J==50,"B"])],predict(lm(STATS[STATS$N==4000&STATS$J==50,"B_sd"]~STATS[STATS$N==4000&STATS$J==50,"B"]+I(STATS[STATS$N==4000&STATS$J==50,"B"]^2)))[order(STATS[STATS$N==4000&STATS$J==50,"B"])],lwd=1.5,col=2)
# legend("top",c("Intercept-1D","Intercept-2D"),pch=19,lty=1,col=c(1,2))
# 
