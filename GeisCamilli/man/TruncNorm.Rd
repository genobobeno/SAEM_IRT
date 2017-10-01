\name{TruncNorm}
\alias{TruncNorm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
TruncNorm(Means, HighLow, prl)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Means}{
%%     ~~Describe \code{Means} here~~
}
  \item{HighLow}{
%%     ~~Describe \code{HighLow} here~~
}
  \item{prl}{
%%     ~~Describe \code{prl} here~~
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
function (Means, HighLow, prl) 
{
    N <- length(Means)
    if (prl) {
        if (HighLow == "high") {
            Q <- unlist(mclapply(Means, function(x) {
                if (x > 4) {
                  LB = 1e-07 + pnorm(0 - x)
                }
                else {
                  LB = pnorm(0 - x)
                }
                return(x + qnorm(runif(1, LB, 1) - 1e-08))
            }))
        }
        else {
            Q <- unlist(mclapply(Means, function(x) {
                if (x < (-4)) {
                  UB = pnorm(0 - x) - 1e-07
                }
                else {
                  UB = pnorm(0 - x)
                }
                return(x + qnorm(runif(1, 0, UB) + 1e-08))
            }))
        }
    }
    else {
        if (HighLow == "high") {
            Q <- unlist(lapply(Means, function(x) {
                if (x > 4) {
                  LB = 1e-07 + pnorm(0 - x)
                }
                else {
                  LB = pnorm(0 - x)
                }
                return(x + qnorm(runif(1, LB, 1) - 1e-08))
            }))
        }
        else {
            Q <- unlist(lapply(Means, function(x) {
                if (x < (-4)) {
                  UB = pnorm(0 - x) - 1e-07
                }
                else {
                  UB = pnorm(0 - x)
                }
                return(x + qnorm(runif(1, 0, UB) + 1e-08))
            }))
        }
    }
    return(Q)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
