\name{GIFAEstimate}
\alias{GIFAEstimate}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
GIFAEstimate(aa, bb, zz, tt, settings, gain = NA, w = NA, rp = NA, ez = NA, ezz = NA, EmpT = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{aa}{
%%     ~~Describe \code{aa} here~~
}
  \item{bb}{
%%     ~~Describe \code{bb} here~~
}
  \item{zz}{
%%     ~~Describe \code{zz} here~~
}
  \item{tt}{
%%     ~~Describe \code{tt} here~~
}
  \item{settings}{
%%     ~~Describe \code{settings} here~~
}
  \item{gain}{
%%     ~~Describe \code{gain} here~~
}
  \item{w}{
%%     ~~Describe \code{w} here~~
}
  \item{rp}{
%%     ~~Describe \code{rp} here~~
}
  \item{ez}{
%%     ~~Describe \code{ez} here~~
}
  \item{ezz}{
%%     ~~Describe \code{ezz} here~~
}
  \item{EmpT}{
%%     ~~Describe \code{EmpT} here~~
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
function (aa, bb, zz, tt, settings, gain = NA, w = NA, rp = NA, 
    ez = NA, ezz = NA, EmpT = FALSE) 
{
    if (is.na(gain)) 
        gain <- 1
    ttt <- Sys.time()
    C <- NA
    mT <- colMeans(as.matrix(tt))
    EZ <- colMeans(zz)
    EZZ <- t(zz) \%*\% zz/(nrow(zz) - 1)
    zz <- as.matrix(t(zz))
    aa <- as.matrix(aa)
    bb <- as.matrix(bb)
    N <- ncol(zz)
    J <- nrow(zz)
    ifelse(length(settings$tmu) == 1, L <- t(aa)/as.numeric(settings$tsigma + 
        t(aa) \%*\% aa), L <- ginv(settings$tsigma + t(aa) \%*\% 
        aa) \%*\% t(aa))
    if (tolower(settings$fm) == "camilli" | EmpT) {
        if (is.na(ez)) {
            S <- EZZ - as.matrix(EZ) \%*\% t(as.matrix(EZ))
            ifelse(ncol(aa) > 1, A <- (S \%*\% t(L)) \%*\% ginv(settings$tsigma - 
                L \%*\% aa + L \%*\% S \%*\% t(L)), A <- (S \%*\% t(L))/as.numeric(settings$tsigma - 
                L \%*\% aa + L \%*\% S \%*\% t(L)))
        }
        else {
            EZ <- ez + (EZ - ez) * gain
            EZZ <- ezz + (EZZ - ezz) * gain
            S <- EZZ - as.matrix(EZ) \%*\% t(as.matrix(EZ))
            ifelse(ncol(aa) > 1, A <- (S \%*\% t(L)) \%*\% ginv(settings$tsigma - 
                L \%*\% aa + L \%*\% S \%*\% t(L)), A <- (S \%*\% t(L))/as.numeric(settings$tsigma - 
                L \%*\% aa + L \%*\% S \%*\% t(L)))
        }
        B <- EZ * (-1)
    }
    else if (tolower(settings$fm) == "eigen") {
        if (is.na(ez)) {
            S <- EZZ - as.matrix(EZ) \%*\% t(as.matrix(EZ))
        }
        else {
            EZ <- ez + (EZ - ez) * gain
            EZZ <- ezz + (EZZ - ezz) * gain
            S <- EZZ - as.matrix(EZ) \%*\% t(as.matrix(EZ))
        }
        S = S - diag(J)
        out = eigen(S, symmetric = TRUE)
        Avec = out$values
        Aload = out$vectors
        if (settings$Adim > 1) {
            Aload <- as.matrix(apply(rbind(Aload[, 1:settings$Adim], 
                aa), 2, function(x) (if (mean(x[1:(length(x)/2)] * 
                x[length(x)/2 + 1:(length(x)/2)]) > 0 & sum(x[1:(length(x)/2)]) > 
                0) {
                return(x[1:(length(x)/2)])
            }
            else {
                return(-1 * x[1:(length(x)/2)])
            })), nrow(EZZ), settings$Adim)
            A = Aload[, 1:settings$Adim] \%*\% sqrt(diag(Avec[1:settings$Adim]))
        }
        else {
            ifelse(mean(Aload[, 1] * aa) > 0, Aload <- Aload[, 
                1], Aload <- (-1) * Aload[, 1])
            A = as.matrix(Aload * sqrt(Avec[1]))
        }
        B <- EZ * (-1)
    }
    else if (tolower(settings$fm) == "new") {
        if (is.na(ez)) {
            S <- EZZ - as.matrix(EZ) \%*\% t(as.matrix(EZ))
            out = princomp(covmat = S, scores = FALSE, cor = FALSE)
            Avec = out$sd^2
            if (settings$Adim > 1) {
                Aload <- as.matrix(apply(rbind(out$loadings[, 
                  1:settings$Adim], aa), 2, function(x) (if (mean(x[1:(length(x)/2)] * 
                  x[length(x)/2 + 1:(length(x)/2)]) > 0) {
                  return(x[1:(length(x)/2)])
                }
                else {
                  return(-1 * x[1:(length(x)/2)])
                })), nrow(EZZ), settings$Adim)
                A = Aload[, 1:settings$Adim] \%*\% sqrt(diag(Avec[1:settings$Adim]))
            }
            else {
                ifelse(mean(out$loadings[, 1] * aa) > 0, Aload <- out$loadings[, 
                  1], Aload <- (-1) * out$loadings[, 1])
                A = as.matrix(Aload * sqrt(Avec[1]))
            }
        }
        else {
            EZ <- ez + (EZ - ez) * gain
            EZZ <- ezz + (EZZ - ezz) * gain
            S <- EZZ - as.matrix(EZ) \%*\% t(as.matrix(EZ))
            out = princomp(covmat = S, scores = FALSE, cor = FALSE)
            Avec = out$sd^2
            if (settings$Adim > 1) {
                Aload <- as.matrix(apply(rbind(out$loadings[, 
                  1:settings$Adim], aa), 2, function(x) (if (mean(x[1:(length(x)/2)] * 
                  x[length(x)/2 + 1:(length(x)/2)]) > 0) {
                  return(x[1:(length(x)/2)])
                }
                else {
                  return(-1 * x[1:(length(x)/2)])
                })), nrow(EZZ), settings$Adim)
                A = Aload[, 1:settings$Adim] \%*\% sqrt(diag(Avec[1:settings$Adim]))
            }
            else {
                ifelse(mean(out$loadings[, 1] * aa) > 0, Aload <- out$loadings[, 
                  1], Aload <- (-1) * out$loadings[, 1])
                A = as.matrix(Aload * sqrt(Avec[1]))
            }
        }
        B <- EZ * (-1)
    }
    else if (tolower(settings$fm) == "old") {
        if (is.na(ez)) {
            S <- EZZ - as.matrix(EZ) \%*\% t(as.matrix(EZ))
            A <- as.matrix(fa(S, nfactors = ncol(as.matrix(aa)), 
                fm = "ml", covar = TRUE, rotate = "none", scores = "none")$loadings)
        }
        else {
            EZ <- ez + (EZ - ez) * gain
            EZZ <- ezz + (EZZ - ezz) * gain
            S <- EZZ - as.matrix(EZ) \%*\% t(as.matrix(EZ))
            A <- as.matrix(fa(S, nfactors = ncol(as.matrix(aa)), 
                fm = "ml", covar = TRUE, rotate = "none", scores = "none")$loadings)
        }
        B <- as.vector(A \%*\% as.matrix(mT) - EZ)
    }
    else if (tolower(settings$fm) == "pca") {
        if (is.na(ez)) {
            zb <- zz + bb \%*\% rep(1, ncol(zz))
            SZB <- zb \%*\% t(zb)/(nrow(zz) - 1)
            S <- EZZ - as.matrix(EZ) \%*\% t(as.matrix(EZ))
            out = princomp(covmat = S, scores = FALSE, cor = FALSE)
            Avec = out$sd^2
            Aload <- as.matrix(apply(out$loadings, 2, function(x) (if (mean((x > 
                0) + 0) > 0.5) {
                return(x)
            }
            else {
                return(-1 * x)
            })), nrow(EZZ), settings$Adim)
            A = Aload[, 1:settings$Adim] \%*\% sqrt(diag(Avec[1:settings$Adim] - 
                1))
        }
        else {
            EZ <- ez + (EZ - ez) * gain
            EZZ <- ezz + (EZZ - ezz) * gain
            S <- EZZ - as.matrix(EZ) \%*\% t(as.matrix(EZ))
            out = princomp(covmat = S, scores = FALSE, cor = FALSE)
            Avec = out$sd^2
            Aload <- as.matrix(apply(out$loadings, 2, function(x) (if (mean((x > 
                0) + 0) > 0.5) {
                return(x)
            }
            else {
                return(-1 * x)
            })), nrow(EZZ), settings$Adim)
            A = Aload[, 1:settings$Adim] \%*\% sqrt(diag(Avec[1:settings$Adim] - 
                1))
        }
        B <- EZ * (-1)
    }
    else {
        if (is.na(ez)) {
            S <- EZZ - as.matrix(EZ) \%*\% t(as.matrix(EZ))
            A <- as.matrix(fa(S, nfactors = ncol(aa), fm = tolower(settings$fm), 
                covar = TRUE, rotate = "none", scores = "none")$loadings)
        }
        else {
            EZ <- ez + (EZ - ez) * gain
            EZZ <- ezz + (EZZ - ezz) * gain
            S <- EZZ - as.matrix(EZ) \%*\% t(as.matrix(EZ))
            A <- as.matrix(fa(S, nfactors = ncol(aa), fm = tolower(settings$fm), 
                covar = TRUE, rotate = "none", scores = "none")$loadings)
        }
        B <- EZ * (-1)
    }
    if (!is.na(w) & !is.na(rp)) {
        C <- rep(0, ncol(w))
        s <- rep(0, ncol(w))
        t <- rep(0, ncol(w))
        for (j in 1:ncol(w)) {
            t[j] <- length(which(w[, j] == 0))
            s[j] <- sum(rp[which(w[, j] == 0), j])
            C[j] <- rbeta(1, s[j] + 1, t[j] - s[j] + 1)
            if (C[j] > 0.45) {
                C[j] <- runif(1, min = 0.05, max = 0.3)
            }
        }
    }
    return(list(A = A, B = B, C = C, EZ = EZ, EZZ = EZZ))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
