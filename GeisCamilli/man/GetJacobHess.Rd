\name{GetJacobHess}
\alias{GetJacobHess}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
GetJacobHess(A, B, Z, setting = settings)
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
function (A, B, Z, setting = settings) 
{
    ZZ <- t(Z)
    K = nrow(ZZ)
    N = ncol(ZZ)
    Q = setting$Adim
    B <- matrix(B, K, 1)
    db = matrix(0, K, 1)
    da = matrix(0, K, Q)
    Jac = matrix(0, K * (1 + Q), 1)
    T1 = 0
    T2 = matrix(0, Q, Q)
    ZT <- matrix(0, K, Q)
    dLdXi2 <- matrix(rep(0, ((Q + 1) * K)^2), (Q + 1) * K, (Q + 
        1) * K)
    dLdXi <- NULL
    ZBi <- ZZ + B \%*\% rep(1, N)
    ZB <- matrix(rowSums(ZBi), K, 1)
    IAA <- solve(diag(Q) + t(A) \%*\% A)
    Lamb <- IAA \%*\% t(A)
    for (ii in 1:N) {
        alpha <- Lamb \%*\% ZBi[, ii]
        T1 <- T1 + alpha
        T2 <- T2 + IAA + alpha \%*\% t(alpha)
        ZT <- ZT + matrix(ZBi[, ii], K, 1) \%*\% t(alpha)
    }
    db <- -ZB + A \%*\% T1
    da <- ZT - A \%*\% T2
    for (jj in 1:K) {
        dlda <- da[jj, ]
        dldb <- db[jj]
        dLdXi <- c(dLdXi, dldb, dlda)
    }
    Jac = dLdXi
    for (jj in 1:K) {
        HJ <- rbind(c(-N, T1), cbind(T1, -T2))
        for (ii in 1:nrow(HJ)) dLdXi2[(jj - 1) * (Q + 1) + ii, 
            (jj - 1) * (Q + 1) + 1:(Q + 1)] <- HJ[ii, ]
    }
    H = dLdXi2
    return(list(Jacob = Jac, Hess = H))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
