\name{GetEmpiricalSE}
\alias{GetEmpiricalSE}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
GetEmpiricalSE(FitDATA, rp, IT = settings$EmpIT, estgain = settings$estgain, thinA = settings$thinA, thinB = settings$thinB)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{FitDATA}{
%%     ~~Describe \code{FitDATA} here~~
}
  \item{rp}{
%%     ~~Describe \code{rp} here~~
}
  \item{IT}{
%%     ~~Describe \code{IT} here~~
}
  \item{estgain}{
%%     ~~Describe \code{estgain} here~~
}
  \item{thinA}{
%%     ~~Describe \code{thinA} here~~
}
  \item{thinB}{
%%     ~~Describe \code{thinB} here~~
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
function (FitDATA, rp, IT = settings$EmpIT, estgain = settings$estgain, 
    thinA = settings$thinA, thinB = settings$thinB) 
{
    cat("\n")
    print(paste("Starting", IT, "iterations of Empirical SEs, Thinning A:", 
        thinA, "; Thinning B:", thinB))
    settings = FitDATA$settings
    settings$fm <- "camilli"
    if (settings$Adim > 1) {
        if (!is.na(FitDATA$AR)) {
            A = as.matrix(FitDATA$AR$loadings)
            THat = as.matrix(FitDATA$Trot[, 1:settings$Adim])
        }
        else {
            A = as.matrix(FitDATA$A)
            THat = as.matrix(FitDATA$That[, 1:settings$Adim])
        }
        B = FitDATA$B
        C = FitDATA$C
        if (settings$fm \%in\% c("new", "pca")) 
            settings$fm <- "camilli"
    }
    else {
        A = FitDATA$A
        B = FitDATA$B
        C = FitDATA$C
        THat = FitDATA$That[, 1]
    }
    Alines <- A
    Blines <- B
    prior <- list(tmu = settings$tmu, tsigma = settings$tsigma)
    J <- ncol(rp)
    N <- nrow(rp)
    iDnew = rep(0, (settings$Adim + 1) * J)
    iGnew = mat.or.vec((settings$Adim + 1) * J, (settings$Adim + 
        1) * J)
    iJnew = mat.or.vec((settings$Adim + 1) * J, (settings$Adim + 
        1) * J)
    iDeltak0 <- rep(0, (settings$Adim + 1) * J)
    iDk0 <- mat.or.vec((settings$Adim + 1) * J, (settings$Adim + 
        1) * J)
    iGk0 <- mat.or.vec((settings$Adim + 1) * J, (settings$Adim + 
        1) * J)
    oDnew = rep(0, (settings$Adim + 1) * J)
    oGnew = mat.or.vec((settings$Adim + 1) * J, (settings$Adim + 
        1) * J)
    oJnew = mat.or.vec((settings$Adim + 1) * J, (settings$Adim + 
        1) * J)
    oDeltak0 <- rep(0, (settings$Adim + 1) * J)
    oDk0 <- mat.or.vec((settings$Adim + 1) * J, (settings$Adim + 
        1) * J)
    oGk0 <- mat.or.vec((settings$Adim + 1) * J, (settings$Adim + 
        1) * J)
    XIMean0 <- rep(1, J * (1 + settings$Adim))
    ifelse(settings$guess, W <- DrawW(aa = A, bb = B, cc = C, 
        tt = THat, rp = rp), W <- NA)
    Z <- SampZ(aa = A, bb = B, that = THat, rp = rp, w = W)
    ASEiter <- array(A, dim = c(J, settings$Adim, 1))
    BSEiter <- matrix(B, nrow = J, ncol = 1)
    ifelse(settings$guess, CSEiter <- matrix(C, nrow = J, ncol = 1), 
        CSEiter <- NA)
    ZSEiter <- matrix(Z, nrow = N, ncol = 1)
    TSEiter <- array(THat, dim = c(N, settings$Adim, 1))
    LLSEiter <- vector()
    for (it in 1:IT) {
        if (it\%\%10 == 0) 
            cat(".")
        if (it\%\%100 == 0) 
            cat(":")
        if (it\%\%1000 == 0) 
            cat(" ", it)
        if (settings$guess) {
            W <- DrawW(aa = A, bb = B, cc = C, tt = THat, rp = rp)
        }
        else {
            W <- NA
        }
        Z <- SampZ(aa = A, bb = B, that = THat, rp = rp, w = W)
        LL <- GIFAFullLL(A, B, Z, THat, prior = prior)
        if (settings$guess) {
            PSI <- GIFAEstimate(aa = A, bb = B, zz = Z, tt = THat, 
                settings = settings, w = W, rp = rp, EmpT = TRUE)
        }
        else {
            PSI <- GIFAEstimate(aa = A, bb = B, zz = Z, tt = THat, 
                settings = settings, EmpT = TRUE)
        }
        A <- PSI$A
        B <- PSI$B
        C <- PSI$C
        ASEiter <- abind(ASEiter, as.matrix(A), along = 3)
        BSEiter <- cbind(BSEiter, B)
        if (settings$guess) 
            CSEiter <- cbind(CSEiter, C)
        ZSEiter <- cbind(ZSEiter, Z)
        TSEiter <- abind(TSEiter, as.matrix(THat), along = 3)
        LLSEiter <- c(LLSEiter, LL)
        gain <- 1/it^estgain
        oJH <- GetErrorOgive(A = A, B = B, C = C, TH = THat, 
            Z = Z, RP = rp)
        oJacob <- oJH$Jacob
        oHess <- oJH$Hess
        oDk0 <- oDk0 + gain * (oHess - oDk0)
        oGk0 <- oGk0 + gain * (as.matrix(oJacob) \%*\% t(as.matrix(oJacob)) - 
            oGk0)
        oDeltak0 <- oDeltak0 + gain * (oJacob - oDeltak0)
        iJH <- GetErrorLogitApp(A = A, B = B, C = C, TH = THat, 
            RP = rp)
        iJacob <- iJH$Jacob
        iHess <- iJH$Hess
        iDk0 <- iDk0 + gain * (iHess - iDk0)
        iGk0 <- iGk0 + gain * (as.matrix(iJacob) \%*\% t(as.matrix(iJacob)) - 
            iGk0)
        iDeltak0 <- iDeltak0 + gain * (iJacob - iDeltak0)
        THat <- SampT(aa = A, bb = B, zz = Z, rp = rp, prior = prior)
    }
    MCthin <- MCthinA <- 1:floor(IT/thinA) * thinA
    MCthinB <- 1:floor(IT/thinB) * thinB
    if (settings$Adim > 1) {
        SEA <- as.matrix(apply(ASEiter[, , MCthinA], c(1, 2), 
            sd))
        MEA <- as.matrix(apply(ASEiter[, , MCthinA], c(1, 2), 
            mean))
        MCSA <- sqrt(as.matrix(apply(ASEiter[, , MCthinA], c(1, 
            2), function(x) (initseq(x)$var.pos))))
    }
    else {
        SEA <- as.vector(apply(as.matrix(ASEiter[, 1, MCthinA]), 
            1, sd))
        MEA <- as.vector(apply(as.matrix(ASEiter[, 1, MCthinA]), 
            1, mean))
        MCSA <- sqrt(as.vector(apply(ASEiter[, 1, MCthinA], 1, 
            function(x) (initseq(x)$var.pos))))
    }
    SEB <- as.vector(apply(BSEiter[, MCthinB], 1, sd))
    MEB <- as.vector(apply(BSEiter[, MCthinB], 1, mean))
    MCSB <- sqrt(as.vector(apply(BSEiter[, MCthinB], 1, function(x) (initseq(x)$var.pos))))
    ifelse(settings$guess, SEC <- as.vector(apply(CSEiter[, MCthin], 
        1, sd)), SEC <- NA)
    ifelse(settings$guess, MEC <- as.vector(apply(CSEiter[, MCthin], 
        1, mean)), MEC <- NA)
    ifelse(settings$guess, MCSC <- sqrt(as.vector(apply(CSEiter[, 
        MCthin], 1, function(x) (initseq(x)$var.pos)))), MCSC <- NA)
    if (settings$Adim > 1) {
        SET <- as.matrix(apply(TSEiter[, , MCthin], c(1, 2), 
            sd))
        MET <- as.matrix(apply(TSEiter[, , MCthin], c(1, 2), 
            mean))
        MCST <- sqrt(as.matrix(apply(TSEiter[, , MCthin], c(1, 
            2), function(x) (initseq(x)$var.pos))))
    }
    else {
        SET <- as.vector(apply(as.matrix(TSEiter[, 1, MCthin]), 
            1, sd))
        MET <- as.vector(apply(as.matrix(TSEiter[, 1, MCthin]), 
            1, mean))
        MCST <- sqrt(as.vector(apply(TSEiter[, 1, MCthin], 1, 
            function(x) (initseq(x)$var.pos))))
    }
    if (settings$plots) {
        par(mfrow = c(1, settings$Adim + 1), mar = c(3, 3, 1, 
            1))
        for (i in 1:settings$Adim) {
            plot(ASEiter[1, i, ], type = "n", xlim = c(0, length(ASEiter[1, 
                1, ])), ylim = range(ASEiter), main = paste("A", 
                i))
            for (j in 1:J) {
                lines(ASEiter[j, i, ], col = j)
            }
            abline(h = Alines[, i])
        }
        plot(BSEiter[1, ], type = "n", xlim = c(0, length(BSEiter[1, 
            ])), ylim = range(BSEiter), main = "B")
        for (j in 1:J) {
            lines(BSEiter[j, ], col = j)
        }
        abline(h = Blines)
    }
    oHk0 <- oDk0 + oGk0 - as.matrix(oDeltak0) \%*\% t(as.matrix(oDeltak0))
    iHk0 <- iDk0 + iGk0 - as.matrix(iDeltak0) \%*\% t(as.matrix(iDeltak0))
    VARLMIo <- diag(ginv((-1) * oHk0))
    VARLMIi <- diag(ginv((-1) * iHk0))
    VARLMIo <- matrix(VARLMIo, J, settings$Adim + 1, byrow = T)[, 
        c(2:(settings$Adim + 1), 1)]
    VARLMIi <- matrix(VARLMIi, J, settings$Adim + 1, byrow = T)[, 
        c(2:(settings$Adim + 1), 1)]
    return(list(SEA = SEA, MEA = MEA, MCSA = MCSA, SEB = SEB, 
        MEB = MEB, MCSB = MCSB, SEC = SEC, MEC = MEC, MCSC = MCSC, 
        SET = SET, MET = MET, MCST = MCST, VARLMIi = VARLMIi, 
        VARLMIo = VARLMIo))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
