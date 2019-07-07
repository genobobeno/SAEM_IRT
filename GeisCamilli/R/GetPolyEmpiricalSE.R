GetPolyEmpiricalSE <-
function(FitDATA,rp,IT=settings$EmpIT,estgain=settings$estgain,thinA=settings$thinA,thinB=settings$thinB) {
  # FitDATA<-list(RP=rp,xi=cbind(ARot,b),A=A,AR=AR,B=b,C=C,
  #               tau=d+matrix(rep(b,ncat-1),length(b),ncat-1),
  #               EZ=Z,EZZ=covZ,settings=settings,That=THAT$THETA,Trot=TROT)
  cat("\n")
  print(paste("Starting",IT,"iterations of Empirical SEs, Thinning A:",thinA,"; Thinning B:",thinB))
  settings=FitDATA$settings
  ncat = settings$ncat
  Q<-settings$Adim
  if (Q>1) {
    if (!is.na(FitDATA$AR)) {
      A=as.matrix(FitDATA$AR$loadings)
      THat=as.matrix(FitDATA$Trot[,1:settings$Adim])
    } else {
      A=as.matrix(FitDATA$A)
      THat=as.matrix(FitDATA$That[,1:settings$Adim])
    }
    B=FitDATA$B
    D=FitDATA$tau - matrix(rep(B,ncat-1),length(B),ncat-1)
    if (settings$fm %in% c("new","pca")) settings$fm <- "camilli"
  } else {
    A=FitDATA$A
    B=FitDATA$B
    D=FitDATA$tau - matrix(rep(B,ncat-1),length(B),ncat-1)
    THat=FitDATA$That[,1]
  }
  # Alines<-A
  # Blines<-B
  prior<-list(tmu=settings$tmu,tsigma=settings$tsigma)
  J<-ncol(rp)
  N<-nrow(rp)
  XIMean0<-rep(1,J*(1+settings$Adim))  

  ASEiter<-array(A, dim=c(J,settings$Adim,1))    
  BSEiter<-matrix(B, nrow=J, ncol=1)
  DSEiter<-array(D, dim=c(J,ncat-1,1))
  #TSEiter<-array(THat, dim=c(N,settings$Adim,1))
  LLSEiter<-vector()
  #print(paste("Running Standard Errors using",IT,"iterations"))
  
  IQ <- diag(Q) 
  # Set up n x m x (ncat-1) data array of binary values for
  y.a = data.matrix(rp)
  # cumulative option indicators
  y.b = array(NA,c(J,ncat-1,N))
  for (i in 1:N) { 
    for (j in 1:J) {
      # for missing responses, item propensities = 9 
      if (y.a[i,j]==9)	{
        y.b[j,,i] <- rep(9,(ncat-1))
      } else {
        for (k in 0:(ncat-2))
        {		
          y.b[j,(k+1),i] <- ifelse(y.a[i,j]<=k,0,1) 
        }
      }
    }
  }
  Y <- y.b
  n1cat = (ncat-1)
  MY = lapply(1:(dim(Y)[1]),function(x) (which(Y[x,,]==9,arr.ind = TRUE)))
  MN = lapply(1:(dim(Y)[1]),function(x) (which(Y[x,,]!=9,arr.ind = TRUE)))
  R  = lapply(1:(dim(Y)[1]),function(x) Y[x,,])
  missList = lapply(MY,function(x) (list(miss = nrow(x),nmis = N*n1cat-nrow(x),mcol = (N*n1cat-nrow(x))/n1cat)))
  # if (!is.na(settings$simfile)) {
  #   ifelse(grepl("\\.[Rr][Dd][Aa]",settings$simfile),
  #          simfile<-settings$simfile,
  #          simfile<-paste(settings$simfile,".rda",sep=""))
  #   if (file.exists(simfile)) {
  #     load(file=simfile)
  #   } else {
  #     print("settings$simfile does not exist.")
  #     settings$simfile<-NA
  #   }
  # }
  # Values for the double truncated WrapZ() function  
  indL<-list();indU<-list()
  for (j in 1:J) {
    yA<-y.a[,j]
    ### This missing handling needs to be revisited... 
    # Best solution is a sample from a different kind of augmented Z... thresholds off?
    if (sum(yA==9)>0) {
      yA[yA==9]<-floor(mean(yA[yA!=9]))
    }
    yL <- matrix(yA+1,N,1)
    indL[[j]] <- cbind(1:N,yL); 		indU[[j]] <- cbind(1:N,yL + 1)
  }
  if (settings$parallel) {
    clusterExport(cl,c("Y","Q","n1cat","N","J","MY","MN","R","missList","indL","indU"),envir=environment())
    clusterEvalQ(cl,c("WrapX","WrapT","WrapZ"))
  }
  ATA 		<- t(A)%*%A #*4
  BTB_INV	<- solve(IQ + ATA)
  d<-D; b<-B

  # Now iterate
  # if (settings$record) {
  #   Aiter = array(0,dim=c(J,Q,niter+1))
  #   Aiter[,,1]<-A
  #   Viter = array(0,dim=c(Q,niter+1)) 
  #   Viter[,1]<-eigen(covZ)$values[1:Q]
  #   Biter = array(0,dim=c(J,niter+1)) 
  #   Biter[,1]<-B
  #   Diter = array(0,dim=c(J,n1cat,niter+1))
  #   Diter[,,1]<-D
  # }
  #Tstart<-Sys.time()
  atemp<-list()
  atemp$Avec<-Avec0<-rep(0,Q)
  atemp$Atemp<-A
  prevA=mat.or.vec(J,Q)
  pRotT<-NA
  pRotV<-NA
  for (i in 1:IT) {
    #if (i%%100==0) cat(nproc,i,"|")
    if (i%%10==1) cat(".")
    if (i%%100==1) cat(":")
    if (i%%500==1) cat("\n",i,"\t : ")
    
    if (settings$parallel) {
      X2		<- simplify2array(parSapply(cl,1:J,WrapX,simplify=FALSE,A=A,b=b,d=d,theta=THat), higher=TRUE)	
    } else {
      X2	<- lapply(1:J,function (x) (WrapX1(x,A=A,b=b,d=d,theta=THat,
                                            settings=settings,R=R,MN=MN,
                                            missList=missList,MY=MY)))
      X2  <- simplify2array(X2,higher=TRUE)	
    }
    X3		<- t(apply(X2,c(1,3),mean))
    b <- -(t(t(rowMeans(X3))))
    d <- -(t(apply(X3, 1, scale, scale=FALSE)))
    #b         <- -meanB;    d        <- -meanD   
    #Z		<- apply(X2,c(2,3),mean)
    if (!settings$dbltrunc) {
      Z		<- colMeans(X2)
    } else {
      if (settings$parallel) {
        X1 		<- parSapply(cl,1:J,WrapZ,simplify=FALSE,A=A,b=b,d=d,theta=THat)
        Z		<- simplify2array(X1, higher=TRUE)	
      } else {
        Z 		<- simplify2array(sapply(1:J,function (x) (WrapZ1(x,A=A,b=b,d=d,theta=THat,
                                                              settings=settings,indL=indL,
                                                              indU=indU))),higher=TRUE)
      } 
    }
    covZ	<- cov(Z)
    # print(i)
    # print("covZ")
    # print(covZ)
    # out<-eigen( covZ-diag(J), symmetric=TRUE)
    # print(out)
    # print(out$vectors[,1:Q]%*%sqrt(diag(out$values[1:Q])))
    #if (i>20) {
    prevA <- A 
    Avec0 <- atemp$Avec[1:Q]
    #}
    if (settings$drawA=="lowertriangular") {
      #compue covariances
      covTMC	<- cov(THat)
      covTZMC	<- cov(THat,Z)
      covT		<- covTMC
      covTZ		<- covTZMC
      atemp	<- DrawALowerDiag(covT=covT,covTZ=covTZ,Q,J,N)
    } else if (settings$drawA=="eigen") {
      atemp	<- tryCatch({
        DrawAEigen(covZ = covZ-diag(J),Q)
      }, error = function(e) {
        atemp
      })
    } else {
      atemp	<- DrawA(covZ-diag(J)/n1cat,Q,a=prevA)
    }
    A <- Anew <- as.matrix(atemp$Atemp)
    # print("A")
    # print(A)
    ATA 		<- t(A)%*%A
    BTB_INV	<- solve(IQ + ATA)
    # print("BTB_INV")
    # print(BTB_INV)
    
    if (settings$parallel) {
      if (Q==1) {
        THat		<- thetanew	<- as.matrix(parSapply(cl,1:N,WrapT,A=A,Z = Z,BTB_INV=BTB_INV,
                                                 b=b,dbltrunc=settings$dbltrunc))
      } else {
        THat		<- thetanew	<- t(parSapply(cl,1:N,WrapTmv,A=A,Z = Z,BTB_INV=BTB_INV,
                                         b=b,dbltrunc=settings$dbltrunc))
      }
    } else {
      if (Q==1) {
        THat	<- thetanew	<- as.matrix(sapply(1:N,function (x) (WrapT(x,A=A,Z = Z,BTB_INV=BTB_INV,
                                                                      b=b,dbltrunc=settings$dbltrunc))))
      } else {
        THat	<- thetanew	<- t(sapply(1:N,function (x) (WrapTmv(x,A=A,Z = Z,BTB_INV=BTB_INV,
                                                                b=b,dbltrunc=settings$dbltrunc))))
      }
    }
    ASEiter<-abind(ASEiter,as.matrix(A),along=3)    
    BSEiter<-cbind(BSEiter,B)
    DSEiter<-abind(DSEiter,as.matrix(d),along=3)
    #TSEiter<-abind(TSEiter,as.matrix(THat),along=3)
  }
  MCthin<-MCthinA<-1:floor(IT/thinA)*thinA
  MCthinB<-1:floor(IT/thinB)*thinB
  if (Q>1) {
    SEA<-as.matrix(apply(ASEiter[,,MCthinA],c(1,2),sd))
    MEA<-as.matrix(apply(ASEiter[,,MCthinA],c(1,2),mean))
    MCSA<-sqrt(as.matrix(apply(ASEiter[,,MCthinA],c(1,2),function(x) (initseq(x[is.finite(x)])$var.pos))))
  } else {
    SEA<-as.vector(apply(as.matrix(ASEiter[,1,MCthinA]),1,sd))
    MEA<-as.vector(apply(as.matrix(ASEiter[,1,MCthinA]),1,mean))
    MCSA<-sqrt(as.vector(apply(ASEiter[,1,MCthinA],1,function(x) (initseq(x[is.finite(x)])$var.pos))))
  } 
  SEB<-as.vector(apply(BSEiter[,MCthinB],1,sd))
  MEB<-as.vector(apply(BSEiter[,MCthinB],1,mean))
  MCSB<-sqrt(as.vector(apply(BSEiter[,MCthinB],1,function(x) (initseq(x[is.finite(x)])$var.pos))))
  SED<-as.vector(apply(as.matrix(DSEiter[,1,MCthinB]),1,sd))
  MED<-as.vector(apply(as.matrix(DSEiter[,1,MCthinB]),1,mean))
  MCSD<-sqrt(as.vector(apply(DSEiter[,1,MCthinB],1,function(x) (initseq(x[is.finite(x)])$var.pos))))
  
  # if (Q>1) {
  #   SET<-as.matrix(apply(TSEiter[,,MCthin],c(1,2),sd))
  #   MET<-as.matrix(apply(TSEiter[,,MCthin],c(1,2),mean))
  #   MCST<-sqrt(as.matrix(apply(TSEiter[,,MCthin],c(1,2),function(x) (initseq(x[is.finite(x)])$var.pos))))
  # } else {
  #   SET<-as.vector(apply(as.matrix(TSEiter[,1,MCthin]),1,sd))
  #   MET<-as.vector(apply(as.matrix(TSEiter[,1,MCthin]),1,mean))
  #   MCST<-sqrt(as.vector(apply(TSEiter[,1,MCthin],1,function(x) (initseq(x[is.finite(x)])$var.pos))))
  # } 
  return(list(SEA=SEA,MEA=MEA,MCSA=MCSA,SEB=SEB,MEB=MEB,MCSB=MCSB,SED=SED,MED=MED,MCSD=MCSD))#,
              #SET=SET,MET=MET,MCST=MCST))
}
