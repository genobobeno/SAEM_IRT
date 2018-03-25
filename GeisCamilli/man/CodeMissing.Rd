\name{CodeMissing}
\alias{CodeMissing}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Code Missing Variables as 9
}
\description{
To simplify missing data, this function recodes all missing data as the value 9.
}
\usage{
CodeMissing(rp,settings)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{rp}{The response pattern}
  \item{settings}{The parameters/settings list for the package}
}
\details{
Response Data assumes numerical values of 0 thru K-1 where K is the number of categories in the polythomous data.
}
\value{
The matrix of the response data.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Eugene Geis
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
