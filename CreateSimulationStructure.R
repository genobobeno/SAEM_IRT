
# Set programs options so this can run on autopilot
gen.dir<-"GeneratedFiles"
fit.dir<-"ConvergedModelFits"
filestructure<-FALSE
options(scipen=999)

sim.list<-list("S1"=list(Guessing=FALSE,J=100,K=2,N=5000,Q=1,mType="beta",Reps=50),
               "S2"=list(Guessing=TRUE, J=100,K=2,N=5000,Q=1,mType="beta",Reps=50),
               "S3"=list(Guessing=FALSE,J=100,K=4,N=5000,Q=1,mType="beta",Reps=50),
               "S4"=list(Guessing=FALSE,J=30, K=4,N=5000,Q=3,mType="bifactor",Reps=50), #bifactor
               "S5"=list(Guessing=FALSE,J=30, K=4,N=5000,Q=3,mType="subscale",Reps=50), #subscale
               "S6"=list(Guessing=FALSE,J=100,K=4,N=10000,Q=5,mType="bifactor",Reps=5), #bifactor
               "S7"=list(Guessing=FALSE,J=100,K=4,N=10000,Q=5,mType="subscale",Reps=5), #subscale
               "S8"=list(Guessing=FALSE,J=100,K=4,N=100000,Q=10,mType="bifactor",Reps=5), #bifactor
               "S9"=list(Guessing=FALSE,J=100,K=4,N=100000,Q=10,mType="subscale",Reps=5)) #subscale

if (!"Simulations.rds" %in% dir(gen.dir)) {
# Set Simulation Conditions' Parameters
  if (!gen.dir %in% dir()) dir.create(gen.dir)
  if (!fit.dir %in% dir()) dir.create(fit.dir)
  # for (d in names(sim.list)) {
  #   if (length(dir(gen.dir))==0 | !d %in% dir(gen.dir)) dir.create(paste0(gen.dir,"/",d))
  # }
  saveRDS(sim.list,paste0(gen.dir,"/Simulations.rds"))
} else {
  print("Generated Files directory and Simulations.rds already created")
  if (exists("sim.list")) {
    NS<-do.call(rbind,sim.list)
    OS<-do.call(rbind,readRDS(paste0(gen.dir,"/Simulations.rds")))
    if (sum(unlist(NS)==unlist(OS),na.rm = TRUE)!=length(unlist(OS))) {
      print("But Simulation list has changed")
      saveRDS(sim.list,paste0(gen.dir,"/Simulations.rds"))
    }
  }
  sim.list<-readRDS(paste0(gen.dir,"/Simulations.rds"))
  print("Simulation.rds read into sim.list")
}
#######################################################

SFileString<-function(l,gen,r=NA) {
  f<-paste0(sapply(names(l),function(x) (paste0(substr(x,1,1),l[[x]]))),collapse="_")
  ifelse(gen,paste0("Gen_",f),paste0("Sim_",f,"_",r))
}

# Check directories
for (d in names(sim.list)) {
  ## Check generated directories
  if (length(dir(gen.dir))==0 | !d %in% dir(gen.dir)) { 
    dir.create(paste0(gen.dir,"/",d))
  } else {
    print(paste0(gen.dir,"/",d," already exists."))
  }
  ## Check fit directories
  if (length(dir(fit.dir))==0 | !d %in% dir(fit.dir)) { 
    dir.create(paste0(fit.dir,"/",d))
  } else {
    print(paste0(fit.dir,"/",d," already exists."))
  }
}

