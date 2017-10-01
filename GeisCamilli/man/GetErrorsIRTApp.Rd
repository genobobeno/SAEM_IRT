\name{GetErrorsIRTApp}
\alias{GetErrorsIRTApp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
GetErrorsIRTApp(A, B, TH, RP, setting = settings)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{A}{
%%     ~~Describe \code{A} here~~
}
  \item{B}{
%%     ~~Describe \code{B} here~~
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
function (A, B, TH, RP, setting = settings) 
{
    q <- ncol(as.matrix(TH))
    k <- ncol(RP)
    Jacob <- NULL
    Hess <- matrix(rep(0, ((q + 1) * k)^2), (q + 1) * k, (q + 
        1) * k)
    D = 1.7
    for (i in 1:ncol(RP)) {
        sum1A = 0
        sum1b = 0
        sum2A = 0
        sum2b = 0
        sumAb = 0
        P <- 1/(1 + exp(-(D * (A[i] * TH - B[i]))))
        Q <- 1 - P
        u <- RP[, i]
        sum1A <- sum(D * (u - P) * TH)
        sum1b <- sum(-D * (u - P))
        sum2A <- sum(-D * D * P * Q * TH * TH)
        sum2b <- sum(-D * D * P * Q)
        sumAb <- sum(D * D * P * Q * TH)
        Jacob <- c(Jacob, sum1b, sum1A)
        HJ <- rbind(c(sum2b, sumAb), cbind(sumAb, sum2A))
        for (ii in 1:nrow(HJ)) Hess[(i - 1) * (q + 1) + ii, (i - 
            1) * (q + 1) + 1:(q + 1)] <- HJ[ii, ]
    }
    return(list(Jacob = Jacob, Hess = Hess))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
