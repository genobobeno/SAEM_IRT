\name{GetJacobian}
\alias{GetJacobian}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
GetJacobian(A, B, Z)
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
function (A, B, Z) 
{
    ZZ <- t(Z)
    K = nrow(ZZ)
    N = ncol(ZZ)
    Q = ncol(A)
    B <- matrix(B, K, 1)
    db = matrix(0, K, 1)
    da = matrix(0, K, Q)
    sg = matrix(0, K, K)
    dLdXi <- vector()
    IAA <- ginv(diag(Q) + t(A) \%*\% A)
    Lamb <- IAA \%*\% t(A)
    ZBi <- ZZ + B \%*\% rep(1, N)
    ZB <- matrix(rowSums(ZBi), K, 1)
    db <- -ZB + A \%*\% Lamb \%*\% ZB
    for (i in 1:N) sg = sg + ZBi[, i] \%*\% t(ZBi[, i])
    da <- ((diag(K) - A \%*\% Lamb) \%*\% sg - N * diag(K)) \%*\% t(Lamb)
    for (j in 1:K) {
        dlda <- da[j, ]
        dldb <- db[j]
        dLdXi <- c(dLdXi, dldb, dlda)
    }
    return(dLdXi)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
