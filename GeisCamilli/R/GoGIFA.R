GoGIFA <- function(rp,init=Init,settings=settings,TargetA=NA,timed=NA) {
  library(rlecuyer)		# for rand num generation
  if (get_os()=="windows") library(snow)			# for parallel processing
  library(GPArotation)	# for rotations
  library(mvnfast)		# for function mvrnorm
  library(psych)			# for ML factor analysis
  library(foreach)
  ttt<-Sys.time()
  if (!is.na(timed)[1] && timed$TF) clock<-Timing()
  stopifnot(length(settings$tmu)==settings$Adim,ncol(init$XI)==(1+settings$Adim+(settings$guess+0)))
  prior<-list(tmu=settings$tmu,tsigma=settings$tsigma)
  J<-ncol(rp)
  N<-nrow(rp)
  rp<-data.matrix(rp)
  CV<-mat.or.vec(J,J)+1
  pxi<-J*(1+settings$Adim+settings$guess)
  indL<-indU<-list()
  for (j in 1:J) {
    indL[[j]] <- cbind(1:N,rp[,j]+1)
    indU[[j]] <- cbind(1:N,rp[,j]+2)
  }
  if (settings$parallel) {
    clusterExport(cl,c("J","N","rp","indL","indU"),envir=environment())
    clusterEvalQ(cl,c("pMVNarma"))
  }
  if (settings$plots) {
    ItCC<-vector() 
    ItAC<-vector() 
  }
  if (tolower(settings$fm)=="licai") {
    ZCai<-array(rep(rnorm(N*J),settings$chains), dim=c(N,J,settings$chains))
    TCai<-array(rep(as.vector(init$THat),settings$chains), dim=c(N,settings$Adim,settings$chains))
  } else {
    Z<-matrix(rnorm(N*J),N,J)
    THat<-init$THat
  }
  if (settings$guess) {
    A<-init$XI[,1:(ncol(init$XI)-2)]
    B<-init$XI[,ncol(init$XI)-1]
    C<-init$XI[,ncol(init$XI)]
  } else {
    A<-init$XI[,1:(ncol(init$XI)-1)]
    B<-init$XI[,ncol(init$XI)]
    C<-NA #for placekeeping
  }
  if (settings$record) {
    Aiter<-array(A, dim=c(J,settings$Adim,1))    
    Biter<-matrix(B, nrow=J, ncol=1)
    ifelse(settings$guess,Citer<-matrix(C, nrow=J, ncol=1),Citer<-NA)
    Viter = array(0,dim=c(settings$Adim,1)) 
    Viter[,1]<-rep(1,settings$Adim)
    #Titer<-array(init$THat, dim=c(N,settings$Adim,1))
    LLiter<-vector()
  }
  #Gain constant array
  gain <- GainConstant(settings=settings)
  
  Dnew = rep(0,(settings$Adim+1+settings$guess)*J) # Jacobian Delta
  Gnew = mat.or.vec((settings$Adim+1+settings$guess)*J,
                    (settings$Adim+1+settings$guess)*J) # Hessian Covariance G
  Jnew = mat.or.vec((settings$Adim+1+settings$guess)*J,
                    (settings$Adim+1+settings$guess)*J) # Hessian D
  
  Deltak0<-rep(0,(settings$Adim+1)*J) # Jacobian Delta
  Hess<-array(mat.or.vec((settings$Adim+1+settings$guess)*J,
                         (settings$Adim+1+settings$guess)*J), 
              dim=c((settings$Adim+1+settings$guess)*J,
                    (settings$Adim+1+settings$guess)*J,1))
  Dk0<-mat.or.vec((settings$Adim+1+settings$guess)*J,
                  (settings$Adim+1+settings$guess)*J) # Hessian D
  Gk0<-mat.or.vec((settings$Adim+1+settings$guess)*J,
                  (settings$Adim+1+settings$guess)*J) # Hessian Covariance G
  XIMean0<-rep(1,J*(1+settings$Adim))
  test<-1
  It<-1  
  if (!is.na(settings$simfile)) {
    ifelse(grepl("\\.[Rr][Dd][Aa]",settings$simfile),
           simfile<-settings$simfile,
           simfile<-paste(settings$simfile,".rda",sep=""))
    if (file.exists(simfile)) {
      load(file=simfile)
    } else {
      print("settings$simfile does not exist.")
      settings$simfile<-NA
    }
  }
  Jacob=matrix(0,pxi,1)
  Hess=matrix(rep(0,(pxi)^2),pxi,pxi) 
  JacobT=matrix(0,pxi,1)
  HessT=matrix(rep(0,(pxi)^2),pxi,pxi) 
  Gamma=matrix(rep(0,(pxi)^2),pxi,pxi) 
  GammaT=matrix(rep(0,(pxi)^2),pxi,pxi) 
  W<-NA
  cat("Iterations")
  passing=0
  if (!is.na(timed)[1] && timed$TF) settings$timed.SAEM_Init<-Timing(clock)
  while (passing!=3) {    
    #print(paste(It,"Next Iteration"))
    if (It%%10==1) cat(".")
    if (It%%100==1) cat(":")
    if (It%%500==1) cat("\n",It,"\t : ")
    if (tolower(settings$fm)=="licai") {
      for (i in 1:settings$chains) {
        ZCai[,,i]<-SampZ(aa=A,bb=B,that=as.matrix(TCai[,1:settings$Adim,i]),rp=rp,w=NA)  
        TCai[,1:settings$Adim,i]<-SampT(aa=A,bb=B,zz=ZCai[,,i],rp=rp,prior=prior)
        ifelse(settings$Adim==1,th<-as.vector(TCai[,1:settings$Adim,i]),th<-as.matrix(TCai[,1:settings$Adim,i]))
        JHT<-GetJacobHessTheta(A,B,Z=ZCai[,,i],TH=th)
        JH<-GetJacobHess(A=A,B=B,Z=ZCai[,,i])
        JacobT=JacobT+JHT$Jacob/settings$chains
        HessT=HessT+JHT$Hess/settings$chains
        Jacob=Jacob+JH$Jacob/settings$chains
        Hess=Hess+JH$Hess/settings$chains
      }
      if (It==1) {
        Gamma = Hess*(-1)
        GammaT = HessT*(-1)
      } else {
        ifelse(It>settings$burnin,gain<-1.0/(It-settings$burnin)^settings$estgain,gain<-1) #
        Gamma = Gamma + (Hess*(-1)-Gamma)*gain
        GammaT = GammaT + (HessT*(-1)-GammaT)*gain
        psi = as.vector(as.matrix(t(cbind(B,A))))
        psi = psi + gain*ginv(GammaT)%*%JacobT
        if (It>settings$burnin) {
          J2    <- matrix(Jacob,pxi,1) %*% matrix(Jacob,1,pxi)
          Dnew = Dnew + JacobT
          Gnew = Gnew + Hess
          Jnew = Jnew + J2  
          Dold = Dnew*gain
          Gold = Gnew*gain
          Jold = Jnew*gain
          H1 = Gold
          H2 = Jold - matrix(Dold,pxi,1)%*%matrix(Dold,1,pxi)
          Hk0 = -(H1+H2)
        }        
        PSI<-matrix(psi,nrow=J,ncol=settings$Adim+1,byrow=TRUE)
        test<-c(as.vector(A-PSI[,1+1:settings$Adim]),B-PSI[,1])
        A<-PSI[,1+1:settings$Adim]
        B<-PSI[,1]
      }
      Aiter<-abind(Aiter,as.matrix(A),along=3)    
      Biter<-cbind(Biter,B)
    } else {
      if (settings$guess) {
        W<-DrawW(aa=A,bb=B,cc=C,tt=THat,rp=rp)
        for (j in 1:J) {
          indL[[j]] <- cbind(1:N,W[,j]+1)
          indU[[j]] <- cbind(1:N,W[,j]+2)
        }
        if (settings$parallel) {
          clusterExport(cl,c("indL","indU"),envir=environment())
        }
      }
      # Z<-SampZ(aa=A,bb=B,that=THat,rp=rp,w=W)    
      #print(THat)
      Z<-SampZFast(aa=A,bb=B,that=THat,indL=indL,indU=indU,srp=rp,w=W,prl=settings$parallel)    
      LL<-GIFAFullLL(A,B,Z,THat,prior=prior)
      if (It<=settings$burnin | tolower(settings$est)!="rm") {
        if (settings$guess) {
          PSI<-GIFAEstimate(aa=A,bb=B,zz=Z,tt=THat,settings=settings,gain=gain[It],w=W,rp=rp)    
        } else {
          PSI<-GIFAEstimate(aa=A,bb=B,zz=Z,tt=THat,settings=settings,gain=gain[It])
        }
      } else {
        if (settings$guess) {
          PSI<-GIFAEstimate(aa=A,bb=B,zz=Z,tt=THat,settings=settings,gain=gain[It],w=W,rp=rp,ez=PSI$EZ,ezz=PSI$EZZ)    
        } else {
          #print("Burned")
          PSI<-GIFAEstimate(aa=A,bb=B,zz=Z,tt=THat,settings=settings,gain=gain[It],ez=PSI$EZ,ezz=PSI$EZZ)
        }    
      }
      ##########################################
      A0<-A  # Fill old values for the "while(test)" Test
      B0<-B
      C0<-C
      A<-PSI$A
      B<-PSI$B
      C<-PSI$C
      ####################### inserting the occasional target rotation ######################
      if (settings$target.rotate.slopes) {
        targR<-TargetRotate(settings,TargetA,aa,that,it,mod.it=100)
        A<-targR$A; AR<-targR$AR; r.matrix<-targR$rotate.matrix
      }
      ##################################
      if (settings$record) {
        Aiter<-abind(Aiter,as.matrix(A),along=3)    
        Biter<-cbind(Biter,B)
        if (settings$guess) Citer<-cbind(Citer,C)
        if (settings$Adim>1) Viter<-cbind(Viter,PSI$Avec[1:Q])
        LLiter<-c(LLiter,LL)
        #if (It%%(10*settings$Adim)==0 && settings$plots && ncol(Biter)>40 && settings$Adim>1) PlotsCCF(Aiter,Biter,settings,ItCC,ItAC)
        if (It%%(10*settings$Adim)==0 && settings$plots && ncol(Biter)>settings$burnin+1 && !is.na(settings$simfile)) PlotsChain(Aiter,Biter,settings,A,init)
      }
      if (tolower(settings$est)=="rm") {
        if (tolower(settings$converge)=="a"|grepl("slop",tolower(settings$converge))) {
          test<-(A-A0)
        } else if (grepl("eig.+val",tolower(settings$converge))|grepl("e?(.+)val",tolower(settings$converge))) {
          if (settings$Adim==1) {
            test<-(sum(A) - sum(A0))
          } else if (prod(A==A0)==1) {
            test<-1
          } else {
            test<-(A-A0)
          }
        } else {
          test<-(A-A0)
        }
        passing<-ifelse(max(abs(test))<settings$eps,passing+1,0) 
        JH<-GetErrorLogitApp(A=A,B=B,C=C,TH=THat,RP=rp)
        Jacob <- JH$Jacob
        Hess  <- JH$Hess
        Dk0<-Dk0+gain[It]*(Hess-Dk0)
        Gk0<-Gk0+gain[It]*(as.matrix(Jacob)%*%t(as.matrix(Jacob))-Gk0)
        Deltak0<-Deltak0+gain[It]*(Jacob-Deltak0)
      } else if (tolower(settings$est)=="off") {
        if (settings$record) print("You need to set record=TRUE to use MCMC means to estimate convergence, i.e. for est='off'")
        if (It>settings$burnin) {
          llim<-floor((It-settings$burnin)/3)
          if (settings$Adim>1) {
            mA<-as.vector(apply(Aiter[,,(settings$burnin-10+llim):(It+1)],c(1,2),sum)/(It+12+llim-settings$burnin))
          } else {
            mA<-as.vector(rowMeans(as.matrix(Aiter[,1,(settings$burnin-10+llim):(It+1)])))
          } 
          mB<-rowMeans(Biter[,(settings$burnin-10):(It+1)])
          if (settings$guess) mC<-rowMeans(Citer[,(settings$burnin-10):(It+1)])
          XIMean<-c(mA,mB)
          test<-(abs(XIMean-XIMean0))
          passing<-ifelse(max(abs(test))<settings$eps,passing+1,0) 
          XIMean0<-XIMean
        }
      } else {
        print("You haven't chosen a valid estimation method... set settings$est='off' or 'rm'")
      }
      THat<-SampT(aa=A,bb=B,zz=Z,rp=rp,prior=prior,prl=settings$parallel,cores=settings$cores) 
    }
    It<-It+1
  }
  print(paste("Total Iterations:",It))
  if (!is.na(timed)[1] && timed$TF) {
    settings$timed.SAEM_Cycles<-Timing(clock)-settings$timed.SAEM_Init
    settings$timed.Iterations<-It
    print(paste("Time: ",settings$timed.SAEM_Cycles))
  }
  if (tolower(settings$est)=="off") {
    ifelse(settings$Adim==1,A<-mA,A<-matrix(mA,nrow=J,ncol=settings$Adim))
    B<-mB
    xiError<-NA
    if (settings$guess) C<-mC
  } else if (tolower(settings$est)=="rm") {
    Hk0<-Dk0+Gk0-as.matrix(Deltak0)%*%t(as.matrix(Deltak0))
    xiError<-diag(ginv(-1*Hk0))
    xiError<-matrix(xiError,nrow(init$XI),ncol(init$XI),byrow=T)
  }
  
  if (tolower(settings$fm)=="licai") {
    Z<-as.matrix(ZCai[,,1])
    THat<-as.matrix(TCai[,1:settings$Adim,1])
  }
  ####   Rotations
  # A is A_gen
  # B is estimated loading matrix, W is a weight matrix, the rotation target is bifactor 0s
  # pstT is partially specified target orthogonal rotation
  if (settings$Adim>1 & !is.na(settings$simfile)) {
    ifelse(grepl("\\.[Rr][Dd][Aa]",settings$simfile),
           load(file=settings$simfile),
           load(file=paste(settings$simfile,".rda",sep="")))
    #gen.rp, gen.xi, gen.theta, gen.structure
    RTS<-permn(1:settings$Adim)
    rtest<-vector()
    if (settings$rmethod=="targetT") {
      for (i in 1:length(RTS)) {
        rtest<-c(rtest,sum(abs(gen.xi[,RTS[[i]]]-targetT(A, Tmat=diag(ncol(A)), Target=gen.xi[,RTS[[i]]], normalize=FALSE, eps=1e-5, maxit=1000)$loadings)))
      }    
      Fctr<-RTS[[which.min(rtest)]]
      AR<-targetT(A, Tmat=diag(ncol(A)), Target=gen.xi[,Fctr], normalize=FALSE, eps=1e-5, maxit=1000)
      AR$APermute<-Fctr
      #Rotate Theta via %*%t(Th)
    } else if (settings$rmethod=="pstT") {
      for (i in 1:length(RTS)) {
        WR <- matrix(0,J,settings$Adim)
        WR[which(gen.xi[,1:settings$Adim]==0)] <- 1
        Tmat <- matrix(-1,settings$Adim,settings$Adim)
        Tmat[1,1] <-  1
        rtest<-c(rtest,sum(abs(gen.xi[,RTS[[i]]]-abs(pstT(A, Tmat=Tmat, W=WR, Target=as.matrix(gen.xi[,RTS[[i]]]), normalize=TRUE, eps=1e-8, maxit=1000)$loadings))))
      }    
      Fctr<-RTS[[which.min(rtest)]]
      AR <- pstT(A, Tmat=Tmat, W=WR, Target=as.matrix(gen.xi[,Fctr]), normalize=TRUE, eps=1e-8, maxit=1000)
      rits<-1
      while (min(apply((AR$loadings>0)+0,2,mean))<0.5) {
        rots<-rep(-1,settings$Adim^2)
        sr<-sample(1:settings$Adim^2)
        rots[sr[1:(rits%%(settings$Adim^2)+1)]]<-1
        Tmat<-matrix(rots,settings$Adim,settings$Adim)
        AR <- pstT(A, Tmat=Tmat, W=WR, Target=as.matrix(gen.xi[,Fctr]), normalize=TRUE, eps=1e-8, maxit=1000)
        rits<-rits+1
      }
      AR$APermute<-Fctr
    }
  } else if (settings$Adim>1) {#what is the proper rotation when there's no TARGET loading?
    ROT.Data<-RotateSlopes(slopes = A)
    if (settings$record) {
      if (grepl("\\.[Rr][Dd][Aa]",settings$estfile)) {
        rfilename=paste0(gsub("\\.[Rr][Dd][Aa]","",settings$estfile),"_ROT.rda")
      } else { 
        rfilename=paste0(settings$estfile,"_ROT.rda") 
      }
      #ROT.Data<-list(Bifactor=ARB,Varimax=ARV,Infomax=ARI)
      save(ROT.Data,settings,file=rfilename)
    }
    if (tolower(settings$rmethod) %in% tolower("Bifactor","Oblimin","Varimax","Infomax")) {
      AR<-ROT.Data[[grepl(tolower(settings$rmethod),tolower(names(ROT.Data)))]]
    } else {
      AR<-NA
    }
  } else {
    AR<-NA
  }
  if (settings$guess & settings$Adim>1) {
    #C<-rowMeans(Citer[,It-50:0])
    if (!is.na(AR)[1]) {
      xi=cbind(AR$loadings,B,C)
    } else {xi=cbind(A,B,C)}
  } else if (!settings$guess & settings$Adim>1 & tolower(settings$fm)!="licai") {
    if (!is.na(AR)[1]) {
      xi=cbind(AR$loadings,B)
    } else {
      xi=cbind(A,B)
    }
  } else if (settings$guess & settings$Adim==1) {
    #C<-rowMeans(Citer[,It-50:0])
    xi=cbind(A,B,C)
  } else {
    xi=cbind(A,B)
  }
  THAT<-GetThetaHat(aa=A,bb=B,cc=C,rp=rp,tHat=THat,zHat=Z,w=W,prior=prior,setting=settings,
                    indU=indU,indL=indL,RT=AR)
  if (settings$Adim>1 & tolower(settings$fm)!="licai" & !is.na(AR)[1] & 
      (!is.na(settings$nesttheta) | settings$thetamap)) {
    if (!is.na(settings$nesttheta)) {
      TROT<-THAT$THETA 
      #TROT<-cbind(THAT$THETA[,Fctr]%*%AR$Th,THAT$THETA[,settings$Adim+Fctr]%*%AR$Th,THAT$THETA[,2*settings$Adim+Fctr]%*%AR$Th)
      Z<-SampZFast(aa=AR$loadings,bb=B,that=TROT[,1:settings$Adim],indL=indL,indU=indU,srp=rp,w=W) 
    } else {
      TROT<-NA
    }
    if (settings$thetamap) {
      TMAPROT<-THAT$TMAP[,Fctr]%*%AR$Th
      Z<-SampZFast(aa=AR$loadings,bb=B,that=TMAPROT[,1:settings$Adim],indL=indL,indU=indU,srp=rp,w=W) 
    } else {
      TMAPROT<-NA
    }
    # Z<-SampZ(aa=AR$loadings,bb=B,that=TROT[,1:settings$Adim],rp=rp,w=W) 
    # Z<-SampZFast(aa=AR$loadings,bb=B,that=TROT[,1:settings$Adim],indL=indL,indU=indU,srp=rp,w=W) 
    oJH<-GetErrorOgive(A=AR$loadings,B=B,C=C,TH=TROT[,1:settings$Adim],Z=Z,RP=rp)
    oJacob <- oJH$Jacob
    oHess  <- oJH$Hess
    oHk0<-oHess-as.matrix(oJacob)%*%t(as.matrix(oJacob))
    iJH<-GetErrorLogitApp(A=AR$loadings,B=B,C=C,TH=TROT[,1:settings$Adim],RP=rp)
    iJacob <- iJH$Jacob
    iHess  <- iJH$Hess
    iHk0<-iHess-as.matrix(iJacob)%*%t(as.matrix(iJacob))
    iError<-diag(ginv(-1*iHk0))
    iError<-matrix(iError,nrow(init$XI),ncol(init$XI),byrow=T)
    oError<-diag(ginv(-1*oHk0))
    oError<-matrix(oError,nrow(init$XI),ncol(init$XI),byrow=T)
    if (is.na(settings$simfile)) {
      gen.xi=NA
      gen.theta=NA
    } else if (length(settings$thetaGen)==length(THAT$THETA[,1:(ncol(THAT$THETA)-2)])) {
      gen.theta=settings$thetaGen
    } else {
      cat("\n didn't pass second check to verify length of thetaGen=THETA\n")
      gen.theta=settings$thetaGen
      cat(length(settings$thetaGen), length(THAT$THETA[,1:(ncol(THAT$THETA)-2)]))
    }
    # FitDATA<-list(XI=gen.xi,RP=rp,THETA=gen.theta,A=A,AR=AR,B=B,C=C,xi=xi,
    #               xiError=xiError,iError=iError,oError=oError,gain=gain,EZ=PSI$EZ,EZZ=PSI$EZZ,
    #               That=THAT$THETA,Tmap=THAT$TMAP,Tmaprot=TMAPROT,TRmap=THAT$TRMAP,
    #               Theta=TROT[,1:settings$Adim],Trot=TROT,settings=settings)
  } else if (settings$Adim>1 & is.na(AR)[1] & 
             (!is.na(settings$nesttheta) | settings$thetamap)) {
    #Z<-SampZ(aa=A,bb=B,that=THAT$THETA[,1:settings$Adim],rp=rp,w=W) 
    if (!is.na(settings$nesttheta) & settings$target.rotate.slopes) {
      Z<-SampZFast(aa=AR$loadings,bb=B,that=THAT$THETA[,1:settings$Adim],indL=indL,indU=indU,srp=rp,w=W) 
      oJH<-GetErrorOgive(A=A,B=B,C=C,TH=THAT$THETA[,1:settings$Adim],Z=Z,RP=rp)
      iJH<-GetErrorLogitApp(A=A,B=B,C=C,TH=THAT$THETA[,1:settings$Adim],RP=rp)
    } else if (!is.na(settings$nesttheta)) {
      Z<-SampZFast(aa=A,bb=B,that=THAT$THETA[,1:settings$Adim],indL=indL,indU=indU,srp=rp,w=W) 
      oJH<-GetErrorOgive(A=A,B=B,C=C,TH=THAT$THETA[,1:settings$Adim],Z=Z,RP=rp)
      iJH<-GetErrorLogitApp(A=A,B=B,C=C,TH=THAT$THETA[,1:settings$Adim],RP=rp)
    } else {
      Z<-SampZFast(aa=AR$loadings,bb=B,that=THAT$TMAP[,1:settings$Adim],indL=indL,indU=indU,srp=rp,w=W) 
      oJH<-GetErrorOgive(A=A,B=B,C=C,TH=THAT$TMAP[,1:settings$Adim],Z=Z,RP=rp)
      iJH<-GetErrorLogitApp(A=A,B=B,C=C,TH=THAT$TMAP[,1:settings$Adim],RP=rp)
    }
    oJacob <- oJH$Jacob
    oHess  <- oJH$Hess
    oHk0<-oHess-as.matrix(oJacob)%*%t(as.matrix(oJacob))
    iJacob <- iJH$Jacob
    iHess  <- iJH$Hess
    iHk0<-iHess-as.matrix(iJacob)%*%t(as.matrix(iJacob))
    iError<-diag(ginv(-1*iHk0))
    iError<-matrix(iError,nrow(init$XI),ncol(init$XI),byrow=T)
    oError<-diag(ginv(-1*oHk0))
    oError<-matrix(oError,nrow(init$XI),ncol(init$XI),byrow=T)
    if (is.na(settings$simfile)) {
      gen.xi=NA
      gen.theta=NA
    } else if (length(settings$thetaGen)==length(THAT$THETA[,1:(ncol(THAT$THETA)-2)])) {
      gen.theta=settings$thetaGen
    }
    TROT<-NA; TMAPROT<-NA
    # FitDATA<-list(XI=gen.xi,RP=rp,THETA=gen.theta,A=A,AR=AR,B=B,C=C,xi=xi,
    #               xiError=xiError,iError=iError,oError=oError,gain=gain,EZ=PSI$EZ,EZZ=PSI$EZZ,
    #               That=THAT$THETA,Tmap=THAT$TMAP,Tmaprot=NA,TRmap=NA,
    #               Theta=NA,Trot=NA,settings=settings)
  } else {
    #Z<-SampZ(aa=A,bb=B,that=THAT$THETA[,1],rp=rp,w=W)
    Z<-SampZFast(aa=A,bb=B,that=THat,indL=indL,indU=indU,srp=rp,w=W)
    oJH<-GetErrorOgive(A=A,B=B,C=C,TH=THat,Z=Z,RP=rp)
    oJacob <- oJH$Jacob
    oHess  <- oJH$Hess
    oHk0<-oHess-as.matrix(oJacob)%*%t(as.matrix(oJacob))
    iJH<-GetErrorLogitApp(A=A,B=B,C=C,TH=THat,RP=rp)
    iJacob <- iJH$Jacob
    iHess  <- iJH$Hess
    iHk0<-iHess-as.matrix(iJacob)%*%t(as.matrix(iJacob))
    iError<-diag(ginv(-1*iHk0))
    iError<-matrix(iError,nrow(init$XI),ncol(init$XI),byrow=T)
    oError<-diag(ginv(-1*oHk0))
    oError<-matrix(oError,nrow(init$XI),ncol(init$XI),byrow=T)
    if (is.na(settings$simfile)) {
      gen.xi=NA
      gen.theta=NA
    } else if (length(settings$thetaGen)==length(THat)) {
      gen.theta=settings$thetaGen
    }
    TROT<-NA; TMAPROT<-NA
    # FitDATA<-list(XI=gen.xi,RP=rp,THETA=gen.theta,A=A,AR=NA,B=B,C=C,xi=xi,
    #               xiError=xiError,iError=iError,oError=oError,gain=gain,EZ=PSI$EZ,EZZ=PSI$EZZ,
    #               That=as.matrix(THat),Tmap=THAT$TMAP,Tmaprot=NA,TRmap=NA,
    #               Theta=NA,Trot=NA,settings=settings)    
  }
  FitDATA<-list(XI=gen.xi,RP=rp,THETA=gen.theta,A=A,AR=AR,B=B,C=C,xi=xi,fTheta=THat,
       gain=gain,EZ=PSI$EZ,EZZ=PSI$EZZ,That=THAT$THETA,Tmap=THAT$TMAP,
       Tmaprot=TMAPROT,TRmap=THAT$TRMAP,Trot=TROT,settings=settings)
  if (settings$esttheta & !is.na(settings$nesttheta)) {
    ThetaFix<-FixedParamTheta(FitDATA,rp=rp,indL=indL,indU=indU)
  } else {ThetaFix<-NA}
  if (settings$empiricalse) {
    EmpSE<-GetEmpiricalSE(FitDATA,rp=rp,indL=indL,indU=indU)
    #print(FitDATA$A)
    # if (settings$Adim>1 & !is.na(AR)[1]) {
    #   FitDATA<-list(XI=gen.xi,RP=rp,THETA=gen.theta,A=A,AR=AR,B=B,C=C,xi=xi,
    #                 xiError=xiError,iError=iError,oError=oError,gain=gain,EZ=PSI$EZ,EZZ=PSI$EZZ,
    #                 That=THAT$THETA,Tmap=THAT$TMAP,Tmaprot=TMAPROT,TRmap=THAT$TRMAP,
    #                 Theta=TROT[,1:settings$Adim],Trot=TROT,EmpSE=EmpSE,ThetaFix=ThetaFix,settings=settings)#      }
    # } else if (settings$Adim>1 & is.na(AR)[1]) {
    #   FitDATA<-list(XI=gen.xi,RP=rp,THETA=gen.theta,A=A,AR=NA,B=B,C=C,xi=xi,
    #                 xiError=xiError,iError=iError,oError=oError,gain=gain,EZ=PSI$EZ,EZZ=PSI$EZZ,
    #                 That=THAT$THETA,Tmap=THAT$TMAP,Tmaprot=NA,TRmap=NA,
    #                 Theta=NA,Trot=NA,EmpSE=EmpSE,ThetaFix=ThetaFix,settings=settings)#      }
    # } else {
    #   FitDATA<-list(XI=gen.xi,RP=rp,THETA=gen.theta,A=A,AR=NA,B=B,C=C,xi=xi,
    #                 xiError=xiError,iError=iError,oError=oError,gain=gain,EZ=PSI$EZ,EZZ=PSI$EZZ,
    #                 That=THAT$THETA,Tmap=THAT$TMAP,Tmaprot=NA,TRmap=NA,
    #                 Theta=NA,Trot=NA,EmpSE=EmpSE,ThetaFix=ThetaFix,settings=settings)    
    # }
  } else {EmpSE<-NA}
  if (settings$record) {
    Iterations<-list(Aiter=Aiter,Viter=Viter,Biter=Biter,Citer=Citer,LLiter=LLiter)
  } else {
    Iterations<-NA
  }
  
  FitDATA<-list(XI=gen.xi,RP=rp,THETA=gen.theta,A=A,AR=AR,B=B,C=C,xi=xi,fTheta=THat,
                xiError=xiError,iError=iError,oError=oError,gain=gain,EZ=PSI$EZ,EZZ=PSI$EZZ,
                That=THAT$THETA,Tmap=THAT$TMAP,Tmaprot=TMAPROT,TRmap=THAT$TRMAP,EmpSE=EmpSE,
                ThetaFix=ThetaFix,Trot=TROT,settings=settings,Iterations=Iterations)
  
  if (!is.na(TROT)[1]) {
    FitDATA$Theta=TROT[,1:settings$Adim]
  } else {FitDATA$Theta=NA}
  ########### Estimate Thetas with fixed parameters
  ## Write files
  if (grepl("\\.[Rr][Dd][Aa]",settings$estfile)) {
    filename=settings$estfile
  } else { filename=paste(settings$estfile,".rda",sep="") }
  save(FitDATA,file=filename)
  FitDATA  
}
