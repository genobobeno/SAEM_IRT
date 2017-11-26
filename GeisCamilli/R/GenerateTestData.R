GenerateTestData <-
function(j=30,n=1000,structure=structure,verbose = FALSE) {
  #j=J;n=N
  stopifnot(tolower(structure$icc) %in% c("logistic","ogive"),j>4,
            length(structure$Aparams)>=1,n>50,structure$Adim==length(structure$tmu))
  if (verbose)  print("***********************  Generating Parameters  **************************")
  #print(unlist(structure))
  A<-GenerateA(j,structure$Adim,tolower(structure$Adist),structure$Aparams)
  if (!is.na(structure$ncat) && structure$ncat!=2) {
    tau = GenerateTau(j,ncat=structure$ncat,taudist=tolower(structure$taudist),tauparams=structure$tauparams)
    b<-as.vector(rowMeans(tau))
  } else {
    tau<-NA
    b<-GenerateB(j,tolower(structure$bdist),structure$bparams)
  }
  if (length(A)==j & verbose) cat(" ...IRT pseudo-difficulties of Ogive: ",b/A)
  xi<-cbind(A,b)
  if (structure$guess) {
    c<-GenerateC(j,tolower(structure$cdist),structure$cparams)
    xi<-cbind(xi,c)
  }  
  t<-GenerateTheta(n,structure$tmu,structure$tsigma)
  rp<-GenerateRP(xi=xi,theta=t,structure=structure,tau=tau)
  gen.rp<-rp
  gen.xi<-xi
  gen.theta<-t
  gen.structure<-structure
  gen.tau<-tau
  if (!is.na(structure$simfile)) {
    if (grepl("\\.[Rr][Dd][Aa]",structure$simfile)) {
      save(gen.rp,gen.xi,gen.theta,gen.tau,gen.structure,file=structure$simfile)
    } else {save(gen.rp,gen.xi,gen.theta,gen.tau,gen.structure,file=paste(structure$simfile,".rda",sep=""))}
  }
  return(list(RP=rp,XI=xi,THETA=t,TAU=tau))
}
