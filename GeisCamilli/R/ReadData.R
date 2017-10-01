ReadData <-
function(GParam=FALSE) {
  TTrue<-as.matrix(read.table("Thetas.dat",header=F))
  A_Gen<-as.matrix(read.table("Aparams.dat",header=F))  
  B_Gen<-as.matrix(read.table("Bparams.dat",header=F))
  if (GParam==TRUE) {
    C_Gen<-as.matrix(read.table("Cparams.dat",header=F))  
  }
  Response<-read.table("Response.dat",header=F)
  if (GParam==FALSE)  {
    RList<-list(A_Gen=A_Gen,B_Gen=B_Gen,Response=Response,TTrue=TTrue)
  } else {
    RList<-list(A_Gen=A_Gen,B_Gen=B_Gen,C_Gen=C_Gen,Response=Response,TTrue=TTrue)
  }
  return(RList)
}
