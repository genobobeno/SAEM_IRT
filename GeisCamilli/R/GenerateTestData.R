GenerateTestData <-
function(j=30,n=1000,structure=structure) {
  j=30;n=500
  structure=list(icc="ogive",          # Item Char Curve; "ogive" or "logistic" 
                 Adist="beta",         # prior distribution of A's/loadings
                 Aparams=c(0.2,1.7),   # parameters of A's/loadings' prior distribution
                 Adim=1,               # 1 (univariate) or 2, 3, etc. multiple dimensions for multivariate
                 bdist="norm",         # distribution of B/intercept
                 bparams=c(0,1),     # parameters of B
                 guess=FALSE,          # guessing ? TRUE/FALSE
                 ncat=3,
                 taudist="unif",
                 tauparams=c(-1,1),
                 cdist="unif",         # guessing parameter distribution for 3PNO or 3PL
                 cparams=c(0.05,0.3),  # bounds
                 tmu= 0,  #c(0,0),                # Theta Prior... e.g. 0, or multivariate c(0,0) # can be multidimensional
                 tsigma= 1) #matrix(c(1,0,0,1),nrow=2,ncol=2))
  
  stopifnot(tolower(structure$icc) %in% c("logistic","ogive"),j>4,
            length(structure$Aparams)>=1,n>50,structure$Adim==length(structure$tmu))
  print("***********************  Generating Parameters  **************************")
  #print(unlist(structure))
  A<-GenerateA(j,structure$Adim,tolower(structure$Adist),structure$Aparams)
  b<-GenerateB(j,tolower(structure$bdist),structure$bparams)
  if (length(A)==j) cat(" ...IRT pseudo-difficulties of Ogive: ",b/A)
  xi<-cbind(A,b)
  if (structure$guess) {
    c<-GenerateC(j,tolower(structure$cdist),structure$cparams)
    xi<-cbind(xi,c)
  }  
  t<-GenerateTheta(n,structure$tmu,structure$tsigma)
  if (!is.na(structure$ncat) && structure$ncat!=2) {
    tau = GenerateTau(j,ncat=structure$ncat,taudist=tolower(structure$taudist),tauparams=structure$tauparams)
  } else {
    tau<-NA
  }
  rp<-GenerateRP(xi=xi,theta=t,structure=structure,tau=tau)
  gen.rp<-rp
  gen.xi<-xi
  gen.theta<-t
  gen.structure<-structure
  if (!is.na(structure$simfile)) {
    if (grepl("\\.[Rr][Dd][Aa]",structure$simfile)) {
      save(gen.rp,gen.xi,gen.theta,gen.structure,file=structure$simfile)
    } else {save(gen.rp,gen.xi,gen.theta,gen.structure,file=paste(structure$simfile,".rda",sep=""))}
  }
  return(list(RP=rp,XI=xi,THETA=t))  
}
