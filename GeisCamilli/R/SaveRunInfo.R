SaveRunInfo <-
function(gen.xi,gen.theta,gen.RP,xi,xiErr,theta,thetaErr,settings=settings,structure=structure,tag) {
  SimData=list(gen.xi=gen.xi,gen.theta=gen.theta,gen.RP=gen.RP,
               xi=xi,xiErr=xiErr,theta=theta,thetaErr=thetaErr,settings=settings,structure=structure)
  save(SimData,file=paste("D",settings$dimF,"_J",nrow(xi),"_N",nrow(as.matrix(theta)),"_",tag,".Rda",sep=""))
  print(paste("Saved file: D",settings$dimF,"_J",nrow(xi),"_N",nrow(as.matrix(theta)),"_",tag,".Rda",sep=""))
}
