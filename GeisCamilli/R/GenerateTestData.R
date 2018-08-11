GenerateTestData <-
function(j=30,n=1000,structure=structure,gen.list=NA,verbose = FALSE) {
  #j=J;n=N
  #cat("start")
  stopifnot(tolower(structure$icc) %in% c("logistic","ogive"),j>4,
            length(structure$Aparams)>=1,n>50,structure$Adim==length(structure$tmu))
  if (verbose)  print("***********************  Generating Parameters  **************************")
  #print(unlist(structure))
  gen.structure<-structure
  if (is.na(gen.list)[1]) {
    A<-GenerateA(j,structure$Adim,tolower(structure$Adist),structure$Aparams)
    #cat("generated A")
    if (!is.na(structure$ncat) && structure$ncat!=2) {
      tau = GenerateTau(j,ncat=structure$ncat,taudist=tolower(structure$taudist),tauparams=structure$tauparams)
      b<-as.vector(rowMeans(tau))
    } else {
      tau<-NA
      b<-GenerateB(j,tolower(structure$bdist),structure$bparams)
    }
    #cat("generated Tau and/or B")
    if (length(A)==j & verbose) cat(" ...IRT pseudo-difficulties of Ogive: ",b/A)
    xi<-cbind(A,b)
    if (structure$guess) {
      c<-GenerateC(j,tolower(structure$cdist),structure$cparams)
      xi<-cbind(xi,c)
    }
    gen.xi<-xi
    gen.tau<-tau
  } else {
    gen.xi<-xi<-gen.list$gen.xi
    gen.tau<-tau<-gen.list$gen.tau
  } 
  t<-GenerateTheta(n,structure$tmu,structure$tsigma)
  #cat("generated Theta")
  rp<-GenerateRP(xi=gen.xi,theta=t,structure=structure,tau=gen.tau)
  #cat("generated RP")
  gen.rp<-rp
  gen.theta<-t
  gen.list<-list(gen.xi=gen.xi,gen.tau=gen.tau)
  if (!is.na(structure$simfile)) {
    if (grepl("\\.[Rr][Dd][Aa]",structure$simfile)) {
      save(gen.rp,gen.xi,gen.theta,gen.tau,gen.structure,gen.list,file=structure$simfile)
    } else {save(gen.rp,gen.xi,gen.theta,gen.tau,gen.structure,gen.list,file=paste(structure$simfile,".rda",sep=""))}
    if (grepl("\\.[Rr][Dd][Ss]",structure$simfile)) {
      saveRDS(list(gen.rp=gen.rp,gen.xi=gen.xi,gen.theta=gen.theta,gen.tau=gen.tau,
                   gen.structure=gen.structure,gen.list=gen.list),structure$simfile)
    } else {saveRDS(list(gen.rp=gen.rp,gen.xi=gen.xi,gen.theta=gen.theta,gen.tau=gen.tau,
                         gen.structure=gen.structure,gen.list=gen.list),paste(structure$simfile,".rds",sep=""))}
  }
  #cat("Saved")
  return(list(RP=rp,XI=xi,THETA=t,TAU=tau,gen.list=gen.list))
}
