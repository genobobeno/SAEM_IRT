if (!"knitr" %in% installed.packages()[,1]) { install.packages("knitr",dependencies = TRUE)} 
knitr::opts_chunk$set(echo=FALSE,warning = FALSE)
source("Hello.R", chdir = T)
Hello()
PACK = c(ifelse(get_os()=="windows","snow","parallel"),"foreach",
         ifelse(get_os()=="windows","doSNOW","doParallel"),
         "rlecuyer","GPArotation","mvnfast","psych","modeest","MASS","stringr",
         "Rcpp","RcppArmadillo","devtools","compiler","doRNG","abind",
         "truncdist","truncnorm","modeest","combinat","fastGHQuad","bdsmatrix","mcmc",
         "RMTstat","maptools","quantreg","plyr","ggplot2","Hmisc","fields","gplots",
         "mvtnorm","psych","binom")

### Install packages not already installed ###
for (i in 1:length(PACK)){
  if(PACK[i] %in% rownames(installed.packages()) == FALSE) {
    install.packages(PACK[i],dependencies = TRUE)}
}
required<-lapply(PACK, require, character.only = TRUE)                       # Call the libraries (will deprecate)
Rcpp::sourceCpp("MVNormRand.cpp")                                  # C++ MV Norm Random (fast)
file.sources = list.files(path = "./GeisCamilli/R/",pattern="*.R") # Grab the R functions!
sourced<-sapply(paste0("./GeisCamilli/R/",file.sources),source,.GlobalEnv)  # Source them!
