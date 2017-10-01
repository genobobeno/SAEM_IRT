\name{GetErrorOgive}
\alias{GetErrorOgive}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
GetErrorOgive(A, B, C = NA, TH, Z, RP, setting = settings)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{A}{
%%     ~~Describe \code{A} here~~
}
  \item{B}{
%%     ~~Describe \code{B} here~~
}
  \item{C}{
%%     ~~Describe \code{C} here~~
}
  \item{TH}{
%%     ~~Describe \code{TH} here~~
}
  \item{Z}{
%%     ~~Describe \code{Z} here~~
}
  \item{RP}{
%%     ~~Describe \code{RP} here~~
}
  \item{setting}{
%%     ~~Describe \code{setting} here~~
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
function (A, B, C = NA, TH, Z, RP, setting = settings) 
{
    d <- setting$Adim
    N <- nrow(RP)
    J <- ncol(RP)
    Phi <- dnorm(Z)
    Jacob <- NULL
    Hess <- matrix(0, (d + 1) * J, (d + 1) * J)
    if (d == 1) {
        THM <- as.matrix(TH) \%*\% matrix(rep(1, J), nrow = 1, 
            ncol = J)
    }
    else {
        THA <- array(0, dim = c(N, J, d))
        for (i in 1:d) {
            THA[, , i] <- as.matrix(TH[, i]) \%*\% matrix(1, nrow = 1, 
                ncol = J)
        }
    }
    P <- pnorm(Z)
    dP = P - P * P
    if (d == 1) {
        dA = THM * Phi * (RP - P)/dP
        dA[is.nan(dA)] = 0
        dA = colSums(dA)
        dB = Phi * (P - RP)/dP
        dB[is.nan(dB)] = 0
        dB = colSums(dB)
        ddA = THM * THM * Phi/dP * ((P - RP) * Z + Phi * (1 - 
            2 * P) * (P - RP)/dP - Phi)
        ddA[is.nan(ddA)] = 0
        ddA = colSums(ddA)
        ddB = Phi/dP * (Z * (P - RP) + Phi * (1 - 2 * P) * (P - 
            RP)/dP - Phi)
        ddB[is.nan(ddB)] = 0
        ddB = colSums(ddB)
        dAdB = THM * Phi/dP * ((RP - P) * Z + Phi * (RP - P) * 
            (1 - 2 * P)/dP + Phi)
        dAdB[is.nan(dAdB)] = 0
        dAdB = colSums(dAdB)
        Jacob = as.vector(t(cbind(dB, dA)))
        for (j in 1:J) {
            HJ = rbind(c(ddB[j], dAdB[j]), c(dAdB[j], ddA[j]))
            for (i in 1:nrow(HJ)) Hess[(j - 1) * 2 + i, (j - 
                1) * 2 + 1:2] <- HJ[i, ]
        }
    }
    else {
        dB = Phi * (P - RP)/dP
        dB[is.nan(dB)] = 0
        dB = colSums(dB)
        ddB = Phi/dP * (Z * (P - RP) + Phi * (1 - 2 * P) * (P - 
            RP)/dP - Phi)
        ddB[is.nan(ddB)] = 0
        ddB = colSums(ddB)
        dA = mat.or.vec(J, d)
        ddA = mat.or.vec(J, d)
        dAdB = mat.or.vec(J, d)
        dAdA = mat.or.vec(J, d * (d - 1)/2)
        for (i in 1:d) {
            dAd = THA[, , i] * Phi * (RP - P)/dP
            dAd[is.nan(dAd)] = 0
            dA[, i] = colSums(dAd)
            ddAd = THA[, , i] * THA[, , i] * Phi/dP * ((P - RP) * 
                Z + Phi * (1 - 2 * P) * (P - RP)/dP - Phi)
            ddAd[is.nan(ddAd)] = 0
            ddA[, i] = colSums(ddAd)
            dAdBd = THA[, , i] * Phi/dP * ((RP - P) * Z + Phi * 
                (RP - P) * (1 - 2 * P)/dP + Phi)
            dAdBd[is.nan(dAdBd)] = 0
            dAdB[, i] = colSums(dAdBd)
        }
        DA2 <- combn(1:d, m = 2)
        if (length(DA2) == 2) {
            dAdAd = THA[, , 1] * THA[, , 2] * Phi/dP * (Z * (P - 
                RP) + Phi * (P - RP) * (1 - 2 * P)/dP - Phi)
            dAdAd[is.nan(dAdAd)] = 0
            dAdA = colSums(dAdAd)
        }
        else {
            for (i in 1:ncol(DA2)) {
                dAdAd = THA[, , DA2[1, i]] * THA[, , DA2[2, i]] * 
                  Phi/dP * (Z * (P - RP) + Phi * (P - RP) * (1 - 
                  2 * P)/dP - Phi)
                dAdAd[is.nan(dAdAd)] = 0
                dAdA[, i] = colSums(dAdAd)
            }
        }
        Jacob = as.vector(t(cbind(dB, dA)))
        for (j in 1:J) {
            AJ = diag(ddA[j, ])
            if (length(DA2) == 2) {
                AJ[1, 2] = AJ[2, 1] = dAdA[j]
            }
            else {
                for (i in 1:ncol(DA2)) {
                  AJ[DA2[1, i], DA2[2, i]] <- dAdA[j, i]
                }
            }
            HJ = rbind(c(ddB[j], dAdB[j, ]), cbind(dAdB[j, ], 
                AJ))
            for (i in 1:nrow(HJ)) Hess[(j - 1) * (d + 1) + i, 
                (j - 1) * (d + 1) + 1:(d + 1)] <- HJ[i, ]
        }
    }
    return(list(Jacob = Jacob, Hess = Hess))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
