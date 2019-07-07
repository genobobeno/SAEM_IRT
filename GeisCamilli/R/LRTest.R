LRTest<-function(fit.data1,fit.data2) {
  F1<-GetLikelihood(fit.data1,TRUE)
  F2<-GetLikelihood(fit.data2,TRUE)
  if (F1["DOF"]>F2["DOF"]) {
    dev<-2*(F1["LogL"]-F2["LogL"])
    dof<-F1["DOF"]-F2["DOF"]
    p<-pchisq(dev,df = dof)
  } else {
    dev<-2*(F2["LogL"]-F1["LogL"])
    dof<-F2["DOF"]-F1["DOF"]
    p<-pchisq(dev,df = dof)
  }
  print(paste("Ratio =",dev,": DoF =",dof,": P =",1-p))
  setNames(1-p,"p-Value")
}
