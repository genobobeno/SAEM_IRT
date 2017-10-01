\name{GetJacobHessTheta}
\alias{GetJacobHessTheta}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
GetJacobHessTheta(A, B, Z, TH, setting = settings)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{A}{
%%     ~~Describe \code{A} here~~
}
  \item{B}{
%%     ~~Describe \code{B} here~~
}
  \item{Z}{
%%     ~~Describe \code{Z} here~~
}
  \item{TH}{
%%     ~~Describe \code{TH} here~~
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
function (A, B, Z, TH, setting = settings) 
{
    if (setting$Adim > 1) {
        ZZ <- t(Z)
        TH <- t(TH)
        K = nrow(ZZ)
        N = ncol(ZZ)
        Q = ncol(A)
        p <- (Q + 1) * K
        B <- matrix(B, K, 1)
        da = matrix(0, K, Q)
        Jac = matrix(0, K * (1 + Q), 1)
        T1 = 0
        T2 = matrix(0, Q, Q)
        ZT <- matrix(0, K, Q)
        dLdXi2 <- matrix(rep(0, (p)^2), p, p)
        dLdXi <- NULL
        ZBi <- ZZ + B \%*\% rep(1, N)
        ATi <- A \%*\% TH
        db <- matrix(rowSums(ATi - ZBi), K, 1)
        da <- (ZBi - ATi) \%*\% t(TH)
        dldab <- matrix(rowSums(TH), Q, 1)
        dda <- TH \%*\% t(TH)
        ddb <- N * (-1)
        for (jj in 1:K) {
            dlda <- da[jj, ]
            dldb <- db[jj]
            dLdXi <- c(dLdXi, dldb, dlda)
        }
        for (jj in 1:K) {
            HJ <- rbind(c(-N, dldab), cbind(dldab, -dda))
            for (ii in 1:nrow(HJ)) dLdXi2[(jj - 1) * (Q + 1) + 
                ii, (jj - 1) * (Q + 1) + 1:(Q + 1)] <- HJ[ii, 
                ]
        }
    }
    else {
        ZZ <- t(Z)
        K = nrow(ZZ)
        N = ncol(ZZ)
        TH <- matrix(TH, 1, N)
        B <- matrix(B, K, 1)
        A = matrix(A, K, 1)
        Jac = matrix(0, K * 2, 1)
        dLdXi2 <- matrix(rep(0, (2 * K)^2), 2 * K, 2 * K)
        dLdXi <- NULL
        ZBi <- ZZ + B \%*\% rep(1, N)
        ATi <- A \%*\% TH
        db <- matrix(rowSums(ATi - ZBi), K, 1)
        T1 <- matrix(1, K, 1) \%*\% TH
        da <- matrix(rowSums((ZBi - ATi) * T1, K, 1))
        dldab <- sum(TH)
        dda <- TH \%*\% t(TH)
        ddb <- N * (-1)
        for (jj in 1:K) {
            dlda <- da[jj]
            dldb <- db[jj]
            dLdXi <- c(dLdXi, dldb, dlda)
        }
        for (jj in 1:K) {
            HJ <- rbind(c(-N, dldab), c(dldab, -dda))
            for (ii in 1:nrow(HJ)) dLdXi2[(jj - 1) * 2 + ii, 
                (jj - 1) * 2 + 1:2] <- HJ[ii, ]
        }
    }
    Jac = dLdXi
    H = dLdXi2
    return(list(Jacob = Jac, Hess = H))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
