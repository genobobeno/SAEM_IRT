GetRMSE<-function(condition) {
  source("CreateSimulationStructure.R") #condition<-"S4"
  d<-condition;    simdir<-paste0(gen.dir,"/",d)
  SimList<-readRDS(paste0(simdir,"/",SFileString(sim.list[[d]],gen=TRUE),"_1.rds"))
  gen.xi<-SimList$gen.xi;    fitdir<-paste0(fit.dir,"/",d)
  sims=sim.list[[d]]$Reps;    items=sim.list[[d]]$J;    examinees=sim.list[[d]]$N 
  XIiter<-array(0, dim=c(items,1+SimList$gen.structure$Adim+(SimList$gen.structure$guess+0),sims))
  if (sim.list[[d]]$K>2) {
    TAUiter<-array(0, dim=c(items,sim.list[[d]]$K-1,sims))
  }
  for (i in 1:sims) { #i<-1
    cat(":",i)
    FitDATA<-readRDS(paste0(fitdir,"/",SFileString(sim.list[[d]],gen=FALSE,r = i),".rds"))
    if (sum(colSums(FitDATA$xi[,1:FitDATA$settings$Adim])<0)>0) {
      for (q in which(colSums(FitDATA$xi[,1:FitDATA$settings$Adim])<0)) {
        FitDATA$xi[,q]<-(-1)*FitDATA$xi[,q]
      }
    }
    XIiter[,,i]<-as.matrix(FitDATA$xi[,1:(1+FitDATA$settings$Adim+(FitDATA$settings$guess+0))])   # converged
    if ("tau" %in% names(FitDATA)) {
      TAUiter[,,i]<-as.matrix(FitDATA$tau)
    }
  }
  RMSE<-sXI<-apply(XIiter,c(1,2),sd)  # RMSE
  if (FitDATA$settings$guess) {
    colnames(RMSE)<-c(paste0("A",1:FitDATA$settings$Adim,"_RMSE"),"B_RMSE","C_RMSE")
  } else {
    colnames(RMSE)<-c(paste0("A",1:FitDATA$settings$Adim,"_RMSE"),"B_RMSE")      
  }
  if (sim.list[[d]]$K>2) {
    gen.tau<-SimList$gen.tau
    TAU_RMSE<-sTAU<-apply(TAUiter,c(1,2),sd)  # RMSE
    colnames(TAU_RMSE)<-paste0("Tau",1:(FitDATA$settings$ncat-1),"_RMSE")
    list(RMSE=cbind(RMSE,TAU_RMSE),gen.xi=gen.xi,gen.tau=gen.tau,settings=FitDATA$settings)
  } else {
    list(RMSE=RMSE,gen.xi=gen.xi,settings=FitDATA$settings)
  }
}
