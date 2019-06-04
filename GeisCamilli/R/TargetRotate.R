TargetRotate<-function(settings,TargetA,aa,that,it,mod.it=100) {
  J<-nrow(gen.xi)
  if (settings$Adim>1 & !is.na(settings$simfile) & It%%mod.it==0) {
    ifelse(grepl("\\.[Rr][Dd][Aa]",settings$simfile),
           load(file=settings$simfile),
           load(file=paste(settings$simfile,".rda",sep="")))
    #gen.rp, gen.xi, gen.theta, gen.structure
    RTS<-permn(1:settings$Adim)
    rtest<-vector()
    if (tolower(settings$fm)=="pca") {
      AparPCA<-princomp(TargetA)
      Apar = AparPCA$scores
      print("PCA Scored A matrix")
      print(Apar)
      print("Theta Covariance Structure")
      print(cov(THat))
      print("A matrix as currently estimated")
      print(aa)
    } else {
      Apar<-TargetA
    }
    if (settings$rmethod=="targetT" & tolower(settings$fm)!="pca") {
      for (i in 1:length(RTS)) {
        ATest<-Apar[,RTS[[i]]]
        rtest<-c(rtest,sum(abs(ATest-targetT(aa, Tmat=diag(ncol(aa)), Target=ATest, normalize=FALSE, eps=1e-5, maxit=1000)$loadings)))
        print("Rotating A, permuting:")
        print(RTS[[i]])
        print("Generated:")
        print(ATest)
        print("Rotated A:")
        print(targetT(aa, Tmat=diag(ncol(aa)), Target=ATest, normalize=FALSE, eps=1e-5, maxit=1000)$loadings)
        print(rtest)
      }
      Fctr<-RTS[[which.min(rtest)]]
      AR<-targetT(aa, Tmat=diag(ncol(aa)), Target=Apar[,Fctr], normalize=FALSE, eps=1e-5, maxit=1000)
      AR$APermute<-Fctr
      aa<-AR$loadings
      output <- list(A=aa,AR=AR,rotate.matrix=AR$Th)
      #Rotate Theta via %*%t(Th)
    } else if (settings$rmethod=="pstT" & tolower(settings$fm)!="pca") {
      print("starting pstT rotation")
      for (i in 1:length(RTS)) {
        # A is A_gen
        # B is estimated loading matrix
        # W is a weight matrix. The rotation target is the bifactor 0’s
        # pstT is partially specified target orthogonal rotation
        ATest<-Apar[,RTS[[i]]]
        WR <- matrix(0,J,settings$Adim)
        WR[which(TargetA==0)] <- 1
        Tmat <- matrix(-1,settings$Adim,settings$Adim)
        Tmat[1,1] <-  1
        print("Created WR and Tmat")
        print(WR)
        print(Tmat)
        rtest<-c(rtest,sum(abs(ATest-abs(pstT(aa, Tmat=Tmat, W=WR, Target=as.matrix(ATest), normalize=TRUE, eps=1e-8, maxit=1000)$loadings))))
        # examine mean square residual of loadings. Not too shabby.
        print("Partially specified target rotation, MSE:")
      }
      Fctr<-RTS[[which.min(rtest)]]
      AR <- pstT(aa, Tmat=Tmat, W=WR, Target=as.matrix(Apar[,Fctr]), normalize=TRUE, eps=1e-8, maxit=1000)
      rits<-1
      while (min(apply((AR$loadings>0)+0,2,mean))<0.5) {
        rots<-rep(-1,settings$Adim^2)
        sr<-sample(1:settings$Adim^2)
        rots[sr[1:(rits%%(settings$Adim^2)+1)]]<-1
        Tmat<-matrix(rots,settings$Adim,settings$Adim)
        AR <- pstT(aa, Tmat=Tmat, W=WR, Target=as.matrix(Apar[,Fctr]), normalize=TRUE, eps=1e-8, maxit=1000)
        rits<-rits+1
      }
      print("Rotated")
      print(AR$loadings)
      print(Apar[,Fctr])
      AR$APermute<-Fctr
      print(sqrt(sum((Apar[,Fctr]-AR$loadings)^2/(settings$Adim*J))))
      aa<-AR$loadings
      output <- list(A=aa,AR=AR,rotate.matrix=AR$Th)
    } else if (tolower(settings$rmethod)=="explore") {#what is the proper rotation when there's no TARGET loading?
      AR<-bifactor(aa, Tmat=diag(ncol(aa)), normalize=FALSE, eps=1e-5, maxit=1000)
      aa<-AR$loadings
      output <- list(A=aa,AR=AR,rotate.matrix=AR$Th)
    } else {
      output <- list(A=aa,AR=NA,rotate.matrix=NA)
    }
  }
}