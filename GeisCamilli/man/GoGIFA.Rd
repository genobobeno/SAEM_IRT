\name{GoGIFA}
\alias{GoGIFA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
GoGIFA(rp, init = Init, settings = settings)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{rp}{
%%     ~~Describe \code{rp} here~~
}
  \item{init}{
%%     ~~Describe \code{init} here~~
}
  \item{settings}{
%%     ~~Describe \code{settings} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (rp, init = Init, settings = settings) 
{
    ttt <- Sys.time()
    stopifnot(length(settings$tmu) == settings$Adim, ncol(init$XI) == 
        (1 + settings$Adim + (settings$guess + 0)))
    prior <- list(tmu = settings$tmu, tsigma = settings$tsigma)
    J <- ncol(rp)
    N <- nrow(rp)
    CV <- mat.or.vec(J, J) + 1
    pxi <- J * (1 + settings$Adim)
    if (settings$plots) {
        ItCC <- vector()
        ItAC <- vector()
    }
    if (tolower(settings$fm) == "licai") {
        ZCai <- array(rep(rnorm(N * J), settings$chains), dim = c(N, 
            J, settings$chains))
        TCai <- array(rep(as.vector(init$THat), settings$chains), 
            dim = c(N, settings$Adim, settings$chains))
    }
    else {
        Z <- matrix(rnorm(N * J), N, J)
        THat <- init$THat
    }
    if (settings$guess) {
        A <- init$XI[, 1:(ncol(init$XI) - 2)]
        B <- init$XI[, ncol(init$XI) - 1]
        C <- init$XI[, ncol(init$XI)]
    }
    else {
        A <- init$XI[, 1:(ncol(init$XI) - 1)]
        B <- init$XI[, ncol(init$XI)]
        C <- NA
    }
    if (tolower(settings$record) == "on") {
        Aiter <- array(A, dim = c(J, settings$Adim, 1))
        Biter <- matrix(B, nrow = J, ncol = 1)
        ifelse(settings$guess, Citer <- matrix(C, nrow = J, ncol = 1), 
            Citer <- NA)
        LLiter <- vector()
    }
    RMwindow <- 100
    gain <- c(rep(1, settings$burnin + 1), runif(rep(1, RMwindow/2), 
        min = 1/(1:(RMwindow/2))^settings$estgain, max = 1), 
        1, runif(1, min = 1/((RMwindow/2 + 1):RMwindow)^settings$estgain, 
            max = 1/(1:(RMwindow/2))^settings$estgain), 1/RMwindow:ceiling(2/settings$eps))
    Dnew = rep(0, (settings$Adim + 1) * J)
    Gnew = mat.or.vec((settings$Adim + 1) * J, (settings$Adim + 
        1) * J)
    Jnew = mat.or.vec((settings$Adim + 1) * J, (settings$Adim + 
        1) * J)
    Deltak0 <- rep(0, (settings$Adim + 1) * J)
    Hess <- array(mat.or.vec((settings$Adim + 1) * J, (settings$Adim + 
        1) * J), dim = c((settings$Adim + 1) * J, (settings$Adim + 
        1) * J, 1))
    Dk0 <- mat.or.vec((settings$Adim + 1) * J, (settings$Adim + 
        1) * J)
    Gk0 <- mat.or.vec((settings$Adim + 1) * J, (settings$Adim + 
        1) * J)
    XIMean0 <- rep(1, J * (1 + settings$Adim))
    test <- rep(1, J * (1 + settings$Adim))
    It <- 1
    if (!is.na(settings$simfile)) {
        ifelse(grepl("\\.[Rr][Dd][Aa]", settings$simfile), simfile <- settings$simfile, 
            simfile <- paste(settings$simfile, ".rda", sep = ""))
        if (file.exists(simfile)) {
            load(file = simfile)
        }
        else {
            print("settings$simfile does not exist.")
            settings$simfile <- NA
        }
    }
    Jacob = matrix(0, pxi, 1)
    Hess = matrix(rep(0, (pxi)^2), pxi, pxi)
    JacobT = matrix(0, pxi, 1)
    HessT = matrix(rep(0, (pxi)^2), pxi, pxi)
    Gamma = matrix(rep(0, (pxi)^2), pxi, pxi)
    GammaT = matrix(rep(0, (pxi)^2), pxi, pxi)
    W <- NA
    while (max(abs(test)) > settings$eps) {
        if (It\%\%10 == 0) 
            cat(".")
        if (It\%\%100 == 0) 
            cat(":")
        if (It\%\%1000 == 0) 
            cat(" ", It, " \n")
        if (tolower(settings$fm) == "licai") {
            for (i in 1:settings$chains) {
                ZCai[, , i] <- SampZ(aa = A, bb = B, that = as.matrix(TCai[, 
                  1:settings$Adim, i]), rp = rp, w = NA)
                TCai[, 1:settings$Adim, i] <- SampT(aa = A, bb = B, 
                  zz = ZCai[, , i], rp = rp, prior = prior)
                ifelse(settings$Adim == 1, th <- as.vector(TCai[, 
                  1:settings$Adim, i]), th <- as.matrix(TCai[, 
                  1:settings$Adim, i]))
                JHT <- GetJacobHessTheta(A, B, Z = ZCai[, , i], 
                  TH = th)
                JH <- GetJacobHess(A = A, B = B, Z = ZCai[, , 
                  i])
                JacobT = JacobT + JHT$Jacob/settings$chains
                HessT = HessT + JHT$Hess/settings$chains
                Jacob = Jacob + JH$Jacob/settings$chains
                Hess = Hess + JH$Hess/settings$chains
            }
            if (It == 1) {
                Gamma = Hess * (-1)
                GammaT = HessT * (-1)
            }
            else {
                ifelse(It > settings$burnin, gain <- 1/(It - 
                  settings$burnin)^settings$estgain, gain <- 1)
                Gamma = Gamma + (Hess * (-1) - Gamma) * gain
                GammaT = GammaT + (HessT * (-1) - GammaT) * gain
                psi = as.vector(as.matrix(t(cbind(B, A))))
                psi = psi + gain * ginv(GammaT) \%*\% JacobT
                if (It > settings$burnin) {
                  J2 <- matrix(Jacob, pxi, 1) \%*\% matrix(Jacob, 
                    1, pxi)
                  Dnew = Dnew + JacobT
                  Gnew = Gnew + Hess
                  Jnew = Jnew + J2
                  Dold = Dnew * gain
                  Gold = Gnew * gain
                  Jold = Jnew * gain
                  H1 = Gold
                  H2 = Jold - matrix(Dold, pxi, 1) \%*\% matrix(Dold, 
                    1, pxi)
                  Hk0 = -(H1 + H2)
                }
                PSI <- matrix(psi, nrow = J, ncol = settings$Adim + 
                  1, byrow = TRUE)
                test <- c(as.vector(A - PSI[, 1 + 1:settings$Adim]), 
                  B - PSI[, 1])
                A <- PSI[, 1 + 1:settings$Adim]
                B <- PSI[, 1]
            }
            Aiter <- abind(Aiter, as.matrix(A), along = 3)
            Biter <- cbind(Biter, B)
        }
        else {
            if (settings$guess) 
                W <- DrawW(aa = A, bb = B, cc = C, tt = THat, 
                  rp = rp)
            Z <- SampZ(aa = A, bb = B, that = THat, rp = rp, 
                w = W)
            LL <- GIFAFullLL(A, B, Z, THat, prior = prior)
            if (It <= settings$burnin | tolower(settings$est) != 
                "rm") {
                if (settings$guess) {
                  PSI <- GIFAEstimate(aa = A, bb = B, zz = Z, 
                    tt = THat, settings = settings, gain = gain[It], 
                    w = W, rp = rp)
                }
                else {
                  PSI <- GIFAEstimate(aa = A, bb = B, zz = Z, 
                    tt = THat, settings = settings, gain = gain[It])
                }
            }
            else {
                if (settings$guess) {
                  PSI <- GIFAEstimate(aa = A, bb = B, zz = Z, 
                    tt = THat, settings = settings, gain = gain[It], 
                    w = W, rp = rp, ez = PSI$EZ, ezz = PSI$EZZ)
                }
                else {
                  PSI <- GIFAEstimate(aa = A, bb = B, zz = Z, 
                    tt = THat, settings = settings, gain = gain[It], 
                    ez = PSI$EZ, ezz = PSI$EZZ)
                }
            }
            A0 <- A
            B0 <- B
            C0 <- C
            A <- PSI$A
            B <- PSI$B
            C <- PSI$C
            if (tolower(settings$record) == "on") {
                Aiter <- abind(Aiter, as.matrix(A), along = 3)
                Biter <- cbind(Biter, B)
                if (settings$guess) 
                  Citer <- cbind(Citer, C)
                LLiter <- c(LLiter, LL)
                if (It\%\%(10 * settings$Adim) == 0 && settings$plots && 
                  ncol(Biter) > settings$burnin + 1 && !is.na(settings$simfile)) 
                  PlotsChain(Aiter, Biter, settings, A, init)
            }
            if (tolower(settings$est) == "rm") {
                test <- c(as.vector(A - A0), B - B0)
                JH <- GetErrorLogitApp(A = A, B = B, C = C, TH = THat, 
                  RP = rp)
                Jacob <- JH$Jacob
                Hess <- JH$Hess
                Dk0 <- Dk0 + gain[It] * (Hess - Dk0)
                Gk0 <- Gk0 + gain[It] * (as.matrix(Jacob) \%*\% 
                  t(as.matrix(Jacob)) - Gk0)
                Deltak0 <- Deltak0 + gain[It] * (Jacob - Deltak0)
            }
            else if (tolower(settings$est) == "off") {
                if (tolower(settings$record) != "on") 
                  print("You need to set record='on' to use MCMC means to estimate convergence, i.e. for est='off'")
                if (It > settings$burnin) {
                  llim <- floor((It - settings$burnin)/3)
                  if (settings$Adim > 1) {
                    mA <- as.vector(apply(Aiter[, , (settings$burnin - 
                      10 + llim):(It + 1)], c(1, 2), sum)/(It + 
                      12 + llim - settings$burnin))
                  }
                  else {
                    mA <- as.vector(rowMeans(as.matrix(Aiter[, 
                      1, (settings$burnin - 10 + llim):(It + 
                        1)])))
                  }
                  mB <- rowMeans(Biter[, (settings$burnin - 10):(It + 
                    1)])
                  if (settings$guess) 
                    mC <- rowMeans(Citer[, (settings$burnin - 
                      10):(It + 1)])
                  XIMean <- c(mA, mB)
                  test <- (abs(XIMean - XIMean0))
                  XIMean0 <- XIMean
                }
            }
            else {
                print("You haven't chosen a valid estimation method... set settings$est='off' or 'rm'")
            }
            THat <- SampT(aa = A, bb = B, zz = Z, rp = rp, prior = prior)
        }
        It <- It + 1
    }
    if (tolower(settings$est) == "off") {
        ifelse(settings$Adim == 1, A <- mA, A <- matrix(mA, nrow = J, 
            ncol = settings$Adim))
        B <- mB
        xiError <- NA
        if (settings$guess) 
            C <- mC
    }
    else if (tolower(settings$est) == "rm") {
        Hk0 <- Dk0 + Gk0 - as.matrix(Deltak0) \%*\% t(as.matrix(Deltak0))
        xiError <- diag(ginv(-1 * Hk0))
        xiError <- matrix(xiError, nrow(init$XI), ncol(init$XI), 
            byrow = T)
    }
    if (tolower(settings$fm) == "licai") {
        Z <- as.matrix(ZCai[, , 1])
        THat <- as.matrix(TCai[, 1:settings$Adim, 1])
    }
    if (settings$Adim > 1 & !is.na(settings$simfile)) {
        ifelse(grepl("\\.[Rr][Dd][Aa]", settings$simfile), load(file = settings$simfile), 
            load(file = paste(settings$simfile, ".rda", sep = "")))
        RTS <- permn(1:settings$Adim)
        rtest <- vector()
        if (settings$rmethod == "targetT") {
            for (i in 1:length(RTS)) {
                rtest <- c(rtest, sum(abs(gen.xi[, RTS[[i]]] - 
                  targetT(A, Tmat = diag(ncol(A)), Target = gen.xi[, 
                    RTS[[i]]], normalize = FALSE, eps = 1e-05, 
                    maxit = 1000)$loadings)))
            }
            Fctr <- RTS[[which.min(rtest)]]
            AR <- targetT(A, Tmat = diag(ncol(A)), Target = gen.xi[, 
                Fctr], normalize = FALSE, eps = 1e-05, maxit = 1000)
            AR$APermute <- Fctr
        }
        else if (settings$rmethod == "pstT") {
            for (i in 1:length(RTS)) {
                WR <- matrix(0, J, settings$Adim)
                WR[which(gen.xi[, 1:settings$Adim] == 0)] <- 1
                Tmat <- matrix(-1, settings$Adim, settings$Adim)
                Tmat[1, 1] <- 1
                rtest <- c(rtest, sum(abs(gen.xi[, RTS[[i]]] - 
                  abs(pstT(A, Tmat = Tmat, W = WR, Target = as.matrix(gen.xi[, 
                    RTS[[i]]]), normalize = TRUE, eps = 1e-08, 
                    maxit = 1000)$loadings))))
            }
            Fctr <- RTS[[which.min(rtest)]]
            AR <- pstT(A, Tmat = Tmat, W = WR, Target = as.matrix(gen.xi[, 
                Fctr]), normalize = TRUE, eps = 1e-08, maxit = 1000)
            rits <- 1
            while (min(apply((AR$loadings > 0) + 0, 2, mean)) < 
                0.5) {
                rots <- rep(-1, settings$Adim^2)
                sr <- sample(1:settings$Adim^2)
                rots[sr[1:(rits\%\%(settings$Adim^2) + 1)]] <- 1
                Tmat <- matrix(rots, settings$Adim, settings$Adim)
                AR <- pstT(A, Tmat = Tmat, W = WR, Target = as.matrix(gen.xi[, 
                  Fctr]), normalize = TRUE, eps = 1e-08, maxit = 1000)
                rits <- rits + 1
            }
            AR$APermute <- Fctr
        }
    }
    else if (settings$Adim > 1) {
        print("Doing a bifactor rotation")
        ARB <- bifactorT(A, Tmat = diag(ncol(A)), normalize = FALSE, 
            eps = 1e-05, maxit = 1000)
        print(ARB)
        print("Doing a Varimax rotation")
        ARV <- Varimax(A, Tmat = diag(ncol(A)), normalize = FALSE, 
            eps = 1e-05, maxit = 1000)
        print(ARV)
        print("Doing a Infomax rotation")
        ARI <- infomaxT(A, Tmat = diag(ncol(A)), normalize = FALSE, 
            eps = 1e-05, maxit = 1000)
        print(ARI)
        AR <- NA
    }
    else {
        AR <- NA
    }
    if (settings$guess & settings$Adim > 1) {
        if (!is.na(AR)) {
            xi = cbind(AR$loadings, B, C)
        }
        else {
            xi = cbind(A, B, C)
        }
    }
    else if (!settings$guess & settings$Adim > 1 & tolower(settings$fm) != 
        "licai") {
        if (!is.na(AR)) {
            xi = cbind(AR$loadings, B)
        }
        else {
            xi = cbind(A, B)
        }
    }
    else if (settings$guess & settings$Adim == 1) {
        xi = cbind(A, B, C)
    }
    else {
        xi = cbind(A, B)
    }
    THAT <- GetThetaHat(aa = A, bb = B, cc = C, rp = rp, tHat = THat, 
        zHat = Z, w = W, prior = prior, setting = settings, R = AR)
    if (settings$Adim > 1 & tolower(settings$fm) != "licai" & 
        !is.na(AR)) {
        TROT <- cbind(THAT$THETA[, Fctr] \%*\% AR$Th, THAT$THETA[, 
            settings$Adim + Fctr] \%*\% AR$Th, THAT$THETA[, 2 * 
            settings$Adim + Fctr] \%*\% AR$Th)
        if (settings$thetamap) {
            TMAPROT <- THAT$TMAP[, Fctr] \%*\% AR$Th
        }
        else {
            TMAPROT <- NA
        }
        Z <- SampZ(aa = AR$loadings, bb = B, that = TROT[, 1:settings$Adim], 
            rp = rp, w = W)
        oJH <- GetErrorOgive(A = AR$loadings, B = B, C = C, TH = TROT[, 
            1:settings$Adim], Z = Z, RP = rp)
        oJacob <- oJH$Jacob
        oHess <- oJH$Hess
        oHk0 <- oHess - as.matrix(oJacob) \%*\% t(as.matrix(oJacob))
        iJH <- GetErrorLogitApp(A = AR$loadings, B = B, C = C, 
            TH = TROT[, 1:settings$Adim], RP = rp)
        iJacob <- iJH$Jacob
        iHess <- iJH$Hess
        iHk0 <- iHess - as.matrix(iJacob) \%*\% t(as.matrix(iJacob))
        iError <- diag(ginv(-1 * iHk0))
        iError <- matrix(iError, nrow(init$XI), ncol(init$XI), 
            byrow = T)
        oError <- diag(ginv(-1 * oHk0))
        oError <- matrix(oError, nrow(init$XI), ncol(init$XI), 
            byrow = T)
        if (is.na(settings$simfile)) {
            gen.xi = NA
            gen.theta = NA
        }
        else if (length(settings$thetaGen) == length(THAT$THETA[, 
            1:(ncol(THAT$THETA) - 2)])) {
            gen.theta = settings$thetaGen
        }
        else {
            cat("\n didn't pass second check to verify length of thetaGen=THETA\n")
            gen.theta = settings$thetaGen
            cat(length(settings$thetaGen), length(THAT$THETA[, 
                1:(ncol(THAT$THETA) - 2)]))
        }
        FitDATA <- list(XI = gen.xi, RP = rp, THETA = gen.theta, 
            A = A, AR = AR, B = B, C = C, xi = xi, xiError = xiError, 
            iError = iError, oError = oError, gain = gain, EZ = PSI$EZ, 
            EZZ = PSI$EZZ, That = THAT$THETA, Tmap = THAT$TMAP, 
            Tmaprot = TMAPROT, TRmap = THAT$TRMAP, Theta = TROT[, 
                1:settings$Adim], Trot = TROT, settings = settings)
    }
    else if (settings$Adim > 1 & is.na(AR)) {
        Z <- SampZ(aa = A, bb = B, that = THAT$THETA[, 1:settings$Adim], 
            rp = rp, w = W)
        oJH <- GetErrorOgive(A = A, B = B, C = C, TH = THAT$THETA[, 
            1:settings$Adim], Z = Z, RP = rp)
        oJacob <- oJH$Jacob
        oHess <- oJH$Hess
        oHk0 <- oHess - as.matrix(oJacob) \%*\% t(as.matrix(oJacob))
        iJH <- GetErrorLogitApp(A = A, B = B, C = C, TH = THAT$THETA[, 
            1:settings$Adim], RP = rp)
        iJacob <- iJH$Jacob
        iHess <- iJH$Hess
        iHk0 <- iHess - as.matrix(iJacob) \%*\% t(as.matrix(iJacob))
        iError <- diag(ginv(-1 * iHk0))
        iError <- matrix(iError, nrow(init$XI), ncol(init$XI), 
            byrow = T)
        oError <- diag(ginv(-1 * oHk0))
        oError <- matrix(oError, nrow(init$XI), ncol(init$XI), 
            byrow = T)
        if (is.na(settings$simfile)) {
            gen.xi = NA
            gen.theta = NA
        }
        else if (length(settings$thetaGen) == length(THAT$THETA[, 
            1:(ncol(THAT$THETA) - 2)])) {
            gen.theta = settings$thetaGen
        }
        FitDATA <- list(XI = gen.xi, RP = rp, THETA = gen.theta, 
            A = A, AR = AR, B = B, C = C, xi = xi, xiError = xiError, 
            iError = iError, oError = oError, gain = gain, EZ = PSI$EZ, 
            EZZ = PSI$EZZ, That = THAT$THETA, Tmap = THAT$TMAP, 
            Tmaprot = NA, TRmap = NA, Theta = NA, Trot = NA, 
            settings = settings)
    }
    else {
        Z <- SampZ(aa = A, bb = B, that = THAT$THETA[, 1], rp = rp, 
            w = W)
        oJH <- GetErrorOgive(A = A, B = B, C = C, TH = THAT$THETA[, 
            1], Z = Z, RP = rp)
        oJacob <- oJH$Jacob
        oHess <- oJH$Hess
        oHk0 <- oHess - as.matrix(oJacob) \%*\% t(as.matrix(oJacob))
        iJH <- GetErrorLogitApp(A = A, B = B, C = C, TH = THAT$THETA[, 
            1], RP = rp)
        iJacob <- iJH$Jacob
        iHess <- iJH$Hess
        iHk0 <- iHess - as.matrix(iJacob) \%*\% t(as.matrix(iJacob))
        iError <- diag(ginv(-1 * iHk0))
        iError <- matrix(iError, nrow(init$XI), ncol(init$XI), 
            byrow = T)
        oError <- diag(ginv(-1 * oHk0))
        oError <- matrix(oError, nrow(init$XI), ncol(init$XI), 
            byrow = T)
        if (is.na(settings$simfile)) {
            gen.xi = NA
            gen.theta = NA
        }
        else if (length(settings$thetaGen) == length(THAT$THETA[, 
            1:(ncol(THAT$THETA) - 2)])) {
            gen.theta = settings$thetaGen
        }
        FitDATA <- list(XI = gen.xi, RP = rp, THETA = gen.theta, 
            A = A, AR = NA, B = B, C = C, xi = xi, xiError = xiError, 
            iError = iError, oError = oError, gain = gain, EZ = PSI$EZ, 
            EZZ = PSI$EZZ, That = THAT$THETA, Tmap = THAT$TMAP, 
            Tmaprot = NA, TRmap = NA, Theta = NA, Trot = NA, 
            settings = settings)
    }
    if (settings$empiricalse) {
        EmpSE <- GetEmpiricalSE(FitDATA, rp = rp)
        ThetaFix <- FixedParamTheta(FitDATA, rp = rp)
        if (settings$Adim > 1 & !is.na(AR)) {
            FitDATA <- list(XI = gen.xi, RP = rp, THETA = gen.theta, 
                A = A, AR = AR, B = B, C = C, xi = xi, xiError = xiError, 
                iError = iError, oError = oError, gain = gain, 
                EZ = PSI$EZ, EZZ = PSI$EZZ, That = THAT$THETA, 
                Tmap = THAT$TMAP, Tmaprot = TMAPROT, TRmap = THAT$TRMAP, 
                Theta = TROT[, 1:settings$Adim], Trot = TROT, 
                EmpSE = EmpSE, ThetaFix = ThetaFix, settings = settings)
        }
        else if (settings$Adim > 1 & is.na(AR)) {
            FitDATA <- list(XI = gen.xi, RP = rp, THETA = gen.theta, 
                A = A, AR = NA, B = B, C = C, xi = xi, xiError = xiError, 
                iError = iError, oError = oError, gain = gain, 
                EZ = PSI$EZ, EZZ = PSI$EZZ, That = THAT$THETA, 
                Tmap = THAT$TMAP, Tmaprot = NA, TRmap = NA, Theta = NA, 
                Trot = NA, EmpSE = EmpSE, ThetaFix = ThetaFix, 
                settings = settings)
        }
        else {
            FitDATA <- list(XI = gen.xi, RP = rp, THETA = gen.theta, 
                A = A, AR = NA, B = B, C = C, xi = xi, xiError = xiError, 
                iError = iError, oError = oError, gain = gain, 
                EZ = PSI$EZ, EZZ = PSI$EZZ, That = THAT$THETA, 
                Tmap = THAT$TMAP, Tmaprot = NA, TRmap = NA, Theta = NA, 
                Trot = NA, EmpSE = EmpSE, ThetaFix = ThetaFix, 
                settings = settings)
        }
    }
    if (grepl("\\.[Rr][Dd][Aa]", settings$estfile)) {
        filename = settings$estfile
    }
    else {
        filename = paste(settings$estfile, ".rda", sep = "")
    }
    if (tolower(settings$record) == "on") {
        MCMCDATA <- list(Aiter = Aiter, Biter = Biter, Citer = Citer, 
            LLiter = LLiter)
        save(FitDATA, MCMCDATA, settings, file = filename)
    }
    else {
        save(FitDATA, settings, file = filename)
    }
    return(FitDATA)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
