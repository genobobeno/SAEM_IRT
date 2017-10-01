\name{GetErrorLogitApp}
\alias{GetErrorLogitApp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
GetErrorLogitApp(A, B, C = NA, TH, RP, setting = settings)
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
function (A, B, C = NA, TH, RP, setting = settings) 
{
    d <- setting$Adim
    N <- nrow(RP)
    J <- ncol(RP)
    Jacob <- NULL
    Hess <- matrix(rep(0, ((d + 1) * J)^2), (d + 1) * J, (d + 
        1) * J)
    D = 1.7
    AT <- t(A \%*\% t(as.matrix(TH)))
    Bz <- t(B \%*\% t(as.matrix(rep(1, nrow(as.matrix(TH))))))
    P = 1/(1 + exp(-D * (AT - Bz)))
    Q = 1 - P
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
    if (d == 1) {
        dA = D * colSums((RP - P) * THM)
        dB = D * colSums(P - RP)
        ddA = -D^2 * colSums(THM * THM * P * Q)
        ddB = -D^2 * colSums(P * Q)
        dAdB = D^2 * colSums(P * Q * THM)
        Jacob = as.vector(t(cbind(dB, dA)))
        for (j in 1:J) {
            HJ = rbind(c(ddB[j], dAdB[j]), c(dAdB[j], ddA[j]))
            for (i in 1:nrow(HJ)) Hess[(j - 1) * 2 + i, (j - 
                1) * 2 + 1:2] <- HJ[i, ]
        }
    }
    else {
        dB = D * colSums(P - RP)
        ddB = -D^2 * colSums(P * Q)
        dA = mat.or.vec(J, d)
        ddA = mat.or.vec(J, d)
        dAdB = mat.or.vec(J, d)
        dAdA = mat.or.vec(J, d * (d - 1)/2)
        for (i in 1:d) {
            dA[, i] = D * colSums((RP - P) * THA[, , i])
            ddA[, i] = -D^2 * colSums(THA[, , i] * THA[, , i] * 
                P * Q)
            dAdB[, i] = D^2 * colSums(P * Q * THA[, , i])
        }
        DA2 <- combn(1:d, m = 2)
        if (length(DA2) == 2) {
            dAdA = -D^2 * colSums(THA[, , DA2[1]] * THA[, , DA2[2]] * 
                P * Q)
        }
        else {
            for (i in 1:ncol(DA2)) dAdA[, i] = -D^2 * colSums(THA[, 
                , DA2[1, i]] * THA[, , DA2[2, i]] * P * Q)
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
